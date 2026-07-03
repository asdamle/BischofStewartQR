// Randomized Bischof-Stewart column selection -- MEX backend.
//
// Mirrors matlab_rand/private/bsqr_rand_mfile.m in behaviour (not RNG). The
// kernel never maintains R11^{-1}R12 or column norms for every column; it tracks
// only the running squared inverse Frobenius norm f2 = ||R11^{-1}||_F^2 and
// samples candidate columns in blocks, brings them into the current reduced
// frame with the accumulated reflectors, and keeps those whose increment holds
// f2 under the per-step threshold.
//
// By default (opt.batched) the kernel runs IN-BLOCK: each sampled block is
// brought to the current frame once, then BSQR is run within it, selecting
// columns greedily while the per-step bound holds (the block's w-vectors and
// running norms are downdated incrementally, BLAS-2). One expensive apply yields
// many selections -- O(k^3) overall vs the single-select O(k^4). batched=false
// recovers the one-selection-per-block path. Both respect the same bound. The
// in-block residual-norm downdate carries the same Businger-Golub recompute
// safeguard as the deterministic kernel (opt.norm_recomp_tol, default sqrt(eps));
// the single-select and global-min paths recompute norms from scratch and need none.
//
// Performance design (see docs/RANDOMIZED_BSQR_PLAN.md):
//   * Block apply is BLAS-3 via the compact-WY form Q_nsel' = I - V T' V'
//     (one gemm + trmm + gemm), not a loop of nsel rank-1 updates -- this is
//     what makes the k<<n win materialize.
//   * Sampling touches only the columns actually tested: uniform uses a partial
//     Fisher-Yates with O(1) pivot removal; norm-weighted uses a Fenwick tree
//     so a step costs O(tested * log n) instead of an O(n log n) sort.
//
// The RNG is independent of MATLAB's, so pivot sequences need not match the
// m-file backend; both satisfy the same guarantees on orthonormal-row input.

#include "mex.h"
#include "blas.h"
#include "lapack.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <vector>

namespace {

enum class ThresholdMode { RunningMean, WorstcaseAllowance };
enum class Sampling { Uniform, NormWeighted };
enum class Pick { BestInBlock, First };

struct Options {
    mwSize k = 0;
    mwSize block_size = 0;   // 0 = auto (set to rand_default_block(k) once k is known)
    ThresholdMode threshold_mode = ThresholdMode::RunningMean;
    double slack = 1.0;
    double norm_recomp_tol = std::sqrt(std::numeric_limits<double>::epsilon());
    Sampling sampling = Sampling::NormWeighted;
    Pick pick = Pick::BestInBlock;
    bool has_seed = false;
    unsigned long long seed = 0;
    bool return_r12 = false;
    bool check_finite = true;
    bool batched = true;     // default kernel path: in-block BSQR (multiple picks/block)
};

void fail(const char *id, const char *msg) {
    mexErrMsgIdAndTxt(id, "%s", msg);
}

std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return s;
}

std::string get_string(const mxArray *a, const char *id) {
    const bool is_string_scalar = mxIsClass(a, "string") && mxGetNumberOfElements(a) == 1;
    if (!(mxIsChar(a) || is_string_scalar)) {
        mexErrMsgIdAndTxt(id, "Expected a string value.");
    }
    char *raw = mxArrayToString(a);
    if (!raw) {
        mexErrMsgIdAndTxt(id, "Could not parse string value.");
    }
    std::string out(raw);
    mxFree(raw);
    return out;
}

bool scalar_to_bool(const mxArray *a, const char *id) {
    if (!(mxIsLogicalScalar(a) || (mxIsNumeric(a) && mxIsScalar(a)))) {
        mexErrMsgIdAndTxt(id, "Expected a logical/numeric scalar.");
    }
    const double v = mxIsLogical(a) ? (mxIsLogicalScalarTrue(a) ? 1.0 : 0.0) : mxGetScalar(a);
    return v != 0.0;
}

double get_pos_scalar(const mxArray *a, const char *id) {
    if (!(mxIsNumeric(a) && mxIsScalar(a) && !mxIsComplex(a))) {
        mexErrMsgIdAndTxt(id, "Expected a real numeric scalar.");
    }
    return mxGetScalar(a);
}

void parse_options(const mxArray *A, int nrhs, const mxArray *prhs[], Options &opt) {
    const mwSize m = mxGetM(A);
    const mwSize n = mxGetN(A);
    const mwSize kmax = std::min(m, n);
    opt.k = kmax;

    if (((nrhs - 1) % 2) != 0) {
        fail("bsqr_rand:InvalidOptions", "Name-value options must be provided in pairs.");
    }

    for (int i = 1; i < nrhs; i += 2) {
        const std::string name = to_lower(get_string(prhs[i], "bsqr_rand:InvalidOptionName"));
        const mxArray *val = prhs[i + 1];

        if (name == "k") {
            const double kd = get_pos_scalar(val, "bsqr_rand:InvalidK");
            if (!std::isfinite(kd) || std::floor(kd) != kd || kd < 0.0 ||
                kd > static_cast<double>(kmax)) {
                fail("bsqr_rand:InvalidK", "k must be an integer in [0, min(size(A))].");
            }
            opt.k = static_cast<mwSize>(kd);
        } else if (name == "block_size") {
            const double bd = get_pos_scalar(val, "bsqr_rand:InvalidBlockSize");
            if (!std::isfinite(bd) || bd < 1.0) {
                fail("bsqr_rand:InvalidBlockSize", "block_size must be a positive integer.");
            }
            opt.block_size = static_cast<mwSize>(std::floor(bd + 0.5));
        } else if (name == "threshold_mode") {
            const std::string tm = to_lower(get_string(val, "bsqr_rand:InvalidThresholdMode"));
            if (tm == "running_mean") {
                opt.threshold_mode = ThresholdMode::RunningMean;
            } else if (tm == "worstcase_allowance") {
                opt.threshold_mode = ThresholdMode::WorstcaseAllowance;
            } else {
                fail("bsqr_rand:InvalidThresholdMode",
                     "threshold_mode must be \"running_mean\" or \"worstcase_allowance\".");
            }
        } else if (name == "slack") {
            const double sd = get_pos_scalar(val, "bsqr_rand:InvalidSlack");
            if (!std::isfinite(sd) || sd < 1.0) {
                fail("bsqr_rand:InvalidSlack", "slack must be a finite scalar >= 1.");
            }
            opt.slack = sd;
        } else if (name == "norm_recomp_tol") {
            if (!(mxIsNumeric(val) && mxIsScalar(val) && !mxIsComplex(val))) {
                fail("bsqr_rand:InvalidNormRecompTol", "norm_recomp_tol must be a real numeric scalar.");
            }
            const double t = mxGetScalar(val);
            if (!std::isfinite(t) || t < 0.0 || t > 1.0) {
                fail("bsqr_rand:InvalidNormRecompTol", "norm_recomp_tol must satisfy 0 <= value <= 1.");
            }
            opt.norm_recomp_tol = t;
        } else if (name == "sampling") {
            const std::string sm = to_lower(get_string(val, "bsqr_rand:InvalidSampling"));
            if (sm == "uniform") {
                opt.sampling = Sampling::Uniform;
            } else if (sm == "normweighted") {
                opt.sampling = Sampling::NormWeighted;
            } else {
                fail("bsqr_rand:InvalidSampling", "sampling must be \"uniform\" or \"normweighted\".");
            }
        } else if (name == "pick") {
            const std::string pk = to_lower(get_string(val, "bsqr_rand:InvalidPick"));
            if (pk == "best_in_block") {
                opt.pick = Pick::BestInBlock;
            } else if (pk == "first") {
                opt.pick = Pick::First;
            } else {
                fail("bsqr_rand:InvalidPick", "pick must be \"best_in_block\" or \"first\".");
            }
        } else if (name == "seed") {
            if (mxIsEmpty(val)) {
                opt.has_seed = false;
            } else {
                const double sd = get_pos_scalar(val, "bsqr_rand:InvalidSeed");
                if (!std::isfinite(sd) || sd < 0.0) {
                    fail("bsqr_rand:InvalidSeed", "seed must be a non-negative finite scalar.");
                }
                opt.has_seed = true;
                opt.seed = static_cast<unsigned long long>(sd);
            }
        } else if (name == "return_r12") {
            opt.return_r12 = scalar_to_bool(val, "bsqr_rand:InvalidReturnR12");
        } else if (name == "backend") {
            const std::string b = to_lower(get_string(val, "bsqr_rand:InvalidBackend"));
            if (!(b == "auto" || b == "mfile" || b == "mex")) {
                fail("bsqr_rand:InvalidBackend", "backend must be \"auto\", \"mfile\", or \"mex\".");
            }
        } else if (name == "check_finite") {
            opt.check_finite = scalar_to_bool(val, "bsqr_rand:InvalidCheckFinite");
        } else if (name == "batched") {
            opt.batched = scalar_to_bool(val, "bsqr_rand:InvalidBatched");
        } else {
            fail("bsqr_rand:UnknownOption", "Unknown bsqr_rand option.");
        }
    }
}

double threshold_value(ThresholdMode mode, double f2, mwSize nsel, mwSize k, mwSize n) {
    const double dn = static_cast<double>(n);
    const double di = static_cast<double>(nsel);
    const double dk = static_cast<double>(k);
    if (mode == ThresholdMode::RunningMean) {
        return (f2 + dn - 2.0 * di) / (dk - di);
    }
    const double Fhat_next = (di + 1.0) * (dn - di) / (dk - di);
    return Fhat_next - f2;
}

// Fenwick (BIT) for weighted sampling without replacement: sample an index with
// probability proportional to its current weight in O(log n), remove/restore a
// weight in O(log n). Used only for norm-weighted sampling.
struct Fenwick {
    std::vector<double> t;
    mwSize n = 0;
    int LOG = 0;
    double cur_total = 0.0;

    void init(const std::vector<double> &w) {
        n = static_cast<mwSize>(w.size());
        t.assign(n + 1, 0.0);
        LOG = 0;
        while ((static_cast<mwSize>(1) << (LOG + 1)) <= n) ++LOG;
        cur_total = 0.0;
        // O(n) bottom-up build: seed each node with its own leaf weight (1-indexed)
        // and fold it into its Fenwick parent in a single increasing pass. Children
        // k < i have already contributed t[k] into t[i] by the time i is processed.
        // (The naive alternative -- n separate add() calls -- is O(n log n).)
        for (mwSize i = 1; i <= n; ++i) {
            t[i] += w[i - 1];
            cur_total += w[i - 1];
            const mwSize parent = i + (i & (~i + 1));   // i + lowbit(i)
            if (parent <= n) t[parent] += t[i];
        }
    }
    void add(mwSize i, double d) {
        cur_total += d;
        for (mwSize j = i + 1; j <= n; j += j & (~j + 1)) t[j] += d;
    }
    // Smallest 0-based index whose cumulative weight exceeds x (0 <= x < total).
    mwSize find(double x) const {
        mwSize pos = 0;
        for (int pw = LOG; pw >= 0; --pw) {
            mwSize nx = pos + (static_cast<mwSize>(1) << pw);
            if (nx <= n && t[nx] <= x) {
                pos = nx;
                x -= t[nx];
            }
        }
        return pos;  // 0-based
    }
};

}  // namespace

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs < 1) {
        fail("bsqr_rand:NotEnoughInputs", "bsqr_rand requires at least one input matrix A.");
    }
    if (nlhs > 5) {
        fail("bsqr_rand:TooManyOutputs", "bsqr_rand supports at most 5 outputs.");
    }

    const mxArray *Ain = prhs[0];
    if (!(mxIsDouble(Ain) && !mxIsComplex(Ain) && mxGetNumberOfDimensions(Ain) == 2)) {
        fail("bsqr_rand:InvalidInput", "A must be a real double matrix.");
    }

    const mwSize m = mxGetM(Ain);
    const mwSize n = mxGetN(Ain);
    const double *A = mxGetPr(Ain);  // read-only; columns are gathered, never mutated

    Options opt;
    parse_options(Ain, nrhs, prhs, opt);

    if (nlhs >= 5 && !opt.return_r12) {
        fail("bsqr_rand:R12NotRequested", "The fifth output (R12) requires return_r12=true.");
    }

    if (opt.check_finite) {
        const size_t mn = static_cast<size_t>(m) * static_cast<size_t>(n);
        for (size_t t = 0; t < mn; ++t) {
            if (!std::isfinite(A[t])) {
                fail("bsqr_rand:NonFiniteInput", "A contains non-finite values.");
            }
        }
    }

    const mwSize k = opt.k;
    if (opt.block_size == 0) {   // auto block size (mirrors bsqr_rand_parse_options.m)
        if (opt.batched) {
            opt.block_size = std::max<mwSize>(k, 1);   // batched favours block = k
        } else {                                       // single-select: ceil(k/2) in [16,64]
            opt.block_size = std::min<mwSize>(64, std::max<mwSize>(16, (k + 1) / 2));
        }
    }
    const mwSize block = std::min(opt.block_size, std::max<mwSize>(n, 1));
    const bool weighted = (opt.sampling == Sampling::NormWeighted);

    // --- State -------------------------------------------------------------
    std::vector<double> V(static_cast<size_t>(m) * k, 0.0);    // reflector store (unit diag)
    std::vector<double> T(static_cast<size_t>(k) * k, 0.0);    // compact-WY factor, lda=k
    std::vector<double> R11(static_cast<size_t>(k) * k, 0.0);  // k x k, lda=k
    std::vector<double> tau(k, 0.0);
    std::vector<mwSize> selected(k, 0);
    std::vector<unsigned char> taken(n, 0);
    double f2 = 0.0;

    // Stats
    std::vector<double> st_f2(k, 0.0), st_crit(k, 0.0), st_thr(k, 0.0), st_Fhat(k, 0.0);
    std::vector<double> st_samples(k, 0.0), st_rounds(k, 0.0), st_fallback(k, 0.0);

    // Scratch
    std::vector<double> X(static_cast<size_t>(m) * std::max<mwSize>(block, 1), 0.0);
    std::vector<double> WB(static_cast<size_t>(std::max<mwSize>(k, 1)) * std::max<mwSize>(block, 1), 0.0);
    std::vector<double> Xtop(static_cast<size_t>(std::max<mwSize>(k, 1)) * std::max<mwSize>(block, 1), 0.0);
    std::vector<double> cbuf(std::max<mwSize>(block, 1), 0.0);
    std::vector<double> tcol(std::max<mwSize>(k, 1), 0.0);
    std::vector<double> xcol(m, 0.0);
    std::vector<mwSize> ids(std::max<mwSize>(block, 1), 0);

    // Uniform sampling pool (partial Fisher-Yates with O(1) removal).
    std::vector<mwSize> remaining;
    // Norm-weighted sampling (Fenwick over squared column norms).
    Fenwick bit;
    std::vector<double> g;          // original squared column norms
    std::vector<mwSize> drawn;      // columns drawn (and zeroed) within a step
    if (weighted) {
        g.assign(n, 0.0);
        ptrdiff_t mm = static_cast<ptrdiff_t>(m), inc1 = 1;
        for (mwSize j = 0; j < n; ++j) {
            const double nrm = dnrm2(&mm, const_cast<double *>(&A[static_cast<size_t>(j) * m]), &inc1);
            g[j] = nrm * nrm;
        }
        bit.init(g);
        drawn.reserve(n);
    } else {
        remaining.resize(n);
        std::iota(remaining.begin(), remaining.end(), static_cast<mwSize>(0));
    }

    std::mt19937_64 rng;
    if (opt.has_seed) {
        rng.seed(opt.seed);
    } else {
        std::random_device rd;
        rng.seed((static_cast<unsigned long long>(rd()) << 32) ^ rd());
    }
    std::uniform_real_distribution<double> unif01(0.0, 1.0);

    const double one = 1.0, zero = 0.0;
    ptrdiff_t inc1 = 1;
    const ptrdiff_t ldA = static_cast<ptrdiff_t>(m);
    const ptrdiff_t ldR = static_cast<ptrdiff_t>(k);

    mwSize rem_count = n;  // pool size at the start of the current step (single-select path)

    if (opt.batched) {
        // ===== Batched in-block BSQR =====
        // Each sampled block is brought to the current frame ONCE (one compact-WY
        // apply); we then run BSQR within it, selecting the in-block minimizer
        // greedily and incrementing f2 / re-deriving the threshold after EACH
        // selection, while the per-step bound holds. One expensive apply thus
        // yields many selections, amortizing the dominant cost. Active block
        // columns are kept contiguous in X[:,0..nact-1] (and the parallel
        // w-vector / running-norm buffers) so every in-block update is BLAS-2.
        const mwSize bmax = std::max<mwSize>(block, 1);
        std::vector<double> Wblk(static_cast<size_t>(std::max<mwSize>(k, 1)) * bmax, 0.0);  // in-block w-vecs
        std::vector<double> sblk(bmax, 0.0);     // running squared trailing norms
        std::vector<double> sblk_ref(bmax, 0.0); // last exact sblk (recompute reference)
        std::vector<double> wn2blk(bmax, 0.0);   // ||w||^2 per block column
        std::vector<double> dots(bmax, 0.0);     // w_a . w_pivot (pre-update)
        std::vector<double> betav(bmax, 0.0);    // alpha/diag per rest column
        std::vector<double> hw(bmax, 0.0);       // reflector-apply gemv scratch
        std::vector<mwSize> idblk(bmax, 0);      // original ids of the active block columns
        const ptrdiff_t ldW = static_cast<ptrdiff_t>(k);

        mwSize nsel = 0;
        mwSize rcount = n;       // uniform pool size (valid front of `remaining`)
        mwSize since_last = 0;   // columns sampled since the last selection (fallback trigger)

        while (nsel < k) {
            const mwSize rem = n - nsel;
            if (since_last >= rem) {
                // Global-min fallback. Random blocks keep missing the acceptable
                // columns (e.g. uniform sampling on concentrated-leverage input),
                // so scan ALL remaining columns in chunks for the true minimizer
                // and take it -- never force-accepting a poor column a block
                // happened to miss. One selection, then normal sampling resumes.
                const double theta = threshold_value(opt.threshold_mode, f2, nsel, k, n) * opt.slack;
                double gmin_c = std::numeric_limits<double>::infinity();
                mwSize gmin_id = n, scanned = 0, chunks = 0, j0 = 0;
                while (j0 < n) {
                    mwSize bc = 0;
                    while (j0 < n && bc < block) {
                        if (!taken[j0]) {
                            std::copy(&A[static_cast<size_t>(j0) * m], &A[static_cast<size_t>(j0) * m] + m,
                                      &X[static_cast<size_t>(bc) * m]);
                            ids[bc] = j0; ++bc;
                        }
                        ++j0;
                    }
                    if (bc == 0) break;
                    ++chunks; scanned += bc;
                    if (nsel > 0) {   // X <- Q_nsel' X (compact-WY), then Xtop <- R11^{-1} X(top)
                        const ptrdiff_t ns = static_cast<ptrdiff_t>(nsel), bcc = static_cast<ptrdiff_t>(bc), mm = ldA;
                        char tT = 'T', tN = 'N';
                        dgemm(&tT, &tN, &ns, &bcc, &mm, &one, V.data(), &ldA, X.data(), &ldA, &zero, WB.data(), &ns);
                        char side = 'L', uplo = 'U', tr = 'T', diag = 'N';
                        dtrmm(&side, &uplo, &tr, &diag, &ns, &bcc, &one, T.data(), &ldR, WB.data(), &ns);
                        const double neg1 = -1.0;
                        dgemm(&tN, &tN, &mm, &bcc, &ns, &neg1, V.data(), &ldA, WB.data(), &ns, &one, X.data(), &ldA);
                        for (mwSize t = 0; t < bc; ++t)
                            std::copy(&X[static_cast<size_t>(t) * m], &X[static_cast<size_t>(t) * m] + nsel,
                                      &Xtop[static_cast<size_t>(t) * nsel]);
                        char s2 = 'L', u2 = 'U', t2 = 'N', d2 = 'N';
                        ptrdiff_t nb = static_cast<ptrdiff_t>(nsel);
                        dtrsm(&s2, &u2, &t2, &d2, &ns, &bcc, &one, R11.data(), &ldR, Xtop.data(), &nb);
                    }
                    for (mwSize t = 0; t < bc; ++t) {
                        const double *xc = &X[static_cast<size_t>(t) * m];
                        double rho2 = 0.0;
                        for (mwSize r = nsel; r < m; ++r) rho2 = std::fma(xc[r], xc[r], rho2);
                        double wn2 = 0.0;
                        if (nsel > 0) {
                            const double *wc = &Xtop[static_cast<size_t>(t) * nsel];
                            for (mwSize r = 0; r < nsel; ++r) wn2 = std::fma(wc[r], wc[r], wn2);
                        }
                        const double c = (rho2 > 0.0) ? (1.0 + wn2) / rho2
                                                      : std::numeric_limits<double>::infinity();
                        if (c < gmin_c) { gmin_c = c; gmin_id = ids[t]; std::copy(xc, xc + m, xcol.data()); }
                    }
                }
                if (!std::isfinite(gmin_c)) {
                    fail("bsqr_rand:RankDeficient",
                         "All remaining columns have ~zero residual; the input appears "
                         "rank-deficient for the requested k. Reduce k to the numerical rank.");
                }
                // Reduce the global minimizer (its reduced column is already in xcol).
                for (mwSize r = 0; r < nsel; ++r) R11[r + static_cast<size_t>(nsel) * k] = xcol[r];
                double tau_i = 0.0;
                ptrdiff_t len = static_cast<ptrdiff_t>(m - nsel);
                if (len > 1) dlarfg(&len, &xcol[nsel], &xcol[nsel + 1], &inc1, &tau_i);
                R11[nsel + static_cast<size_t>(nsel) * k] = xcol[nsel];
                V[nsel + static_cast<size_t>(nsel) * m] = 1.0;
                for (mwSize r = nsel + 1; r < m; ++r) V[r + static_cast<size_t>(nsel) * m] = xcol[r];
                tau[nsel] = tau_i;
                if (nsel > 0) {
                    const ptrdiff_t ns = static_cast<ptrdiff_t>(nsel), mm = ldA;
                    char tT = 'T';
                    dgemv(&tT, &mm, &ns, &one, V.data(), &ldA, &V[static_cast<size_t>(nsel) * m], &inc1,
                          &zero, tcol.data(), &inc1);
                    const double negtau = -tau_i;
                    for (mwSize r = 0; r < nsel; ++r) tcol[r] *= negtau;
                    char uplo = 'U', tN = 'N', diag = 'N';
                    dtrmv(&uplo, &tN, &diag, &ns, T.data(), &ldR, tcol.data(), &inc1);
                    for (mwSize r = 0; r < nsel; ++r) T[r + static_cast<size_t>(nsel) * k] = tcol[r];
                }
                T[nsel + static_cast<size_t>(nsel) * k] = tau_i;
                f2 += gmin_c;
                selected[nsel] = gmin_id;
                taken[gmin_id] = 1;
                if (weighted) {
                    bit.add(gmin_id, -g[gmin_id]);
                } else {
                    mwSize w = 0;
                    for (mwSize jj = 0; jj < rcount; ++jj)
                        if (!taken[remaining[jj]]) remaining[w++] = remaining[jj];
                    rcount = w;
                }
                st_f2[nsel] = f2;
                st_crit[nsel] = gmin_c;
                st_thr[nsel] = theta;
                st_Fhat[nsel] = (static_cast<double>(nsel) + 1.0) *
                    (static_cast<double>(n) - nsel) / (static_cast<double>(k) - nsel);
                st_fallback[nsel] = (gmin_c > theta) ? 1.0 : 0.0;
                st_samples[nsel] = static_cast<double>(scanned);
                st_rounds[nsel] = static_cast<double>(chunks);
                since_last = 0;
                ++nsel;
                continue;
            }
            const mwSize bcount = std::min(block, rem);

            // --- draw bcount distinct columns; gather into X (m x bcount) ---
            if (weighted) {
                drawn.clear();
                for (mwSize t = 0; t < bcount; ++t) {
                    mwSize id;
                    const double total = bit.cur_total;
                    if (total <= 0.0) {
                        id = n;
                        for (mwSize j = 0; j < n; ++j) { if (!taken[j]) { id = j; break; } }
                    } else {
                        double x = unif01(rng) * total;
                        if (x >= total) x = std::nextafter(total, 0.0);
                        id = bit.find(x);
                    }
                    bit.add(id, -g[id]);
                    drawn.push_back(id);
                    idblk[t] = id;
                    const double *src = &A[static_cast<size_t>(id) * m];
                    std::copy(src, src + m, &X[static_cast<size_t>(t) * m]);
                }
            } else {
                for (mwSize t = 0; t < bcount; ++t) {
                    std::uniform_int_distribution<mwSize> d(t, rcount - 1);
                    std::swap(remaining[t], remaining[d(rng)]);
                    const mwSize id = remaining[t];
                    idblk[t] = id;
                    const double *src = &A[static_cast<size_t>(id) * m];
                    std::copy(src, src + m, &X[static_cast<size_t>(t) * m]);
                }
            }
            since_last += bcount;

            // --- bring the block to the current frame: X <- Q_nsel' X (compact-WY) ---
            if (nsel > 0) {
                const ptrdiff_t ns = static_cast<ptrdiff_t>(nsel), bc = static_cast<ptrdiff_t>(bcount), mm = ldA;
                char tT = 'T', tN = 'N';
                dgemm(&tT, &tN, &ns, &bc, &mm, &one, V.data(), &ldA, X.data(), &ldA, &zero, WB.data(), &ns);
                char side = 'L', uplo = 'U', tr = 'T', diag = 'N';
                dtrmm(&side, &uplo, &tr, &diag, &ns, &bc, &one, T.data(), &ldR, WB.data(), &ns);
                const double neg1 = -1.0;
                dgemm(&tN, &tN, &mm, &bc, &ns, &neg1, V.data(), &ldA, WB.data(), &ns, &one, X.data(), &ldA);
            }

            // --- initialise in-block state: w-vectors (one trsm) and running norms ---
            if (nsel > 0) {
                for (mwSize t = 0; t < bcount; ++t) {
                    std::copy(&X[static_cast<size_t>(t) * m], &X[static_cast<size_t>(t) * m] + nsel,
                              &Wblk[static_cast<size_t>(t) * k]);
                }
                char side = 'L', uplo = 'U', tr = 'N', diag = 'N';
                ptrdiff_t ns = static_cast<ptrdiff_t>(nsel), bc = static_cast<ptrdiff_t>(bcount);
                dtrsm(&side, &uplo, &tr, &diag, &ns, &bc, &one, R11.data(), &ldR, Wblk.data(), &ldW);
            }
            for (mwSize t = 0; t < bcount; ++t) {
                const double *xc = &X[static_cast<size_t>(t) * m];
                double s = 0.0;
                for (mwSize r = nsel; r < m; ++r) s = std::fma(xc[r], xc[r], s);
                sblk[t] = s;
                sblk_ref[t] = s;   // exact value seeds the recompute reference
                const double *wc = &Wblk[static_cast<size_t>(t) * k];
                double w = 0.0;
                for (mwSize r = 0; r < nsel; ++r) w = std::fma(wc[r], wc[r], w);
                wn2blk[t] = w;
            }

            mwSize nact = bcount;
            const mwSize block_start = nsel;

            // --- in-block greedy BSQR ---
            while (nsel < k && nact > 0) {
                mwSize pbest = 0;
                double cbest = std::numeric_limits<double>::infinity();
                for (mwSize a = 0; a < nact; ++a) {
                    const double c = (sblk[a] > 0.0) ? (1.0 + wn2blk[a]) / sblk[a]
                                                     : std::numeric_limits<double>::infinity();
                    if (c < cbest) { cbest = c; pbest = a; }
                }
                const double theta = threshold_value(opt.threshold_mode, f2, nsel, k, n) * opt.slack;
                if (!std::isfinite(cbest) || cbest > theta) break;   // none acceptable -> resample

                // Swap the pivot to the end of the active region (rest stays contiguous).
                const mwSize nr = nact - 1;
                if (pbest != nr) {
                    std::swap_ranges(&X[static_cast<size_t>(pbest) * m],
                                     &X[static_cast<size_t>(pbest) * m] + m, &X[static_cast<size_t>(nr) * m]);
                    std::swap_ranges(&Wblk[static_cast<size_t>(pbest) * k],
                                     &Wblk[static_cast<size_t>(pbest) * k] + k, &Wblk[static_cast<size_t>(nr) * k]);
                    std::swap(sblk[pbest], sblk[nr]);
                    std::swap(sblk_ref[pbest], sblk_ref[nr]);
                    std::swap(wn2blk[pbest], wn2blk[nr]);
                    std::swap(idblk[pbest], idblk[nr]);
                }

                // Build the reflector from the pivot column (now at position nr).
                std::copy(&X[static_cast<size_t>(nr) * m], &X[static_cast<size_t>(nr) * m] + m, xcol.data());
                for (mwSize r = 0; r < nsel; ++r) R11[r + static_cast<size_t>(nsel) * k] = xcol[r];
                double tau_i = 0.0;
                ptrdiff_t len = static_cast<ptrdiff_t>(m - nsel);
                if (len > 1) dlarfg(&len, &xcol[nsel], &xcol[nsel + 1], &inc1, &tau_i);
                R11[nsel + static_cast<size_t>(nsel) * k] = xcol[nsel];
                V[nsel + static_cast<size_t>(nsel) * m] = 1.0;
                for (mwSize r = nsel + 1; r < m; ++r) V[r + static_cast<size_t>(nsel) * m] = xcol[r];
                tau[nsel] = tau_i;

                // Compact-WY: new column of T (dlarft forward recurrence) for the next apply.
                // The reflector v is zero above row nsel, so V'v only needs rows nsel..m-1.
                if (nsel > 0) {
                    const ptrdiff_t ns = static_cast<ptrdiff_t>(nsel), mlen = static_cast<ptrdiff_t>(m - nsel);
                    char tT = 'T';
                    dgemv(&tT, &mlen, &ns, &one, &V[static_cast<size_t>(nsel)], &ldA,
                          &V[static_cast<size_t>(nsel) + static_cast<size_t>(nsel) * m], &inc1,
                          &zero, tcol.data(), &inc1);
                    const double negtau = -tau_i;
                    for (mwSize r = 0; r < nsel; ++r) tcol[r] *= negtau;
                    char uplo = 'U', tN = 'N', diag = 'N';
                    dtrmv(&uplo, &tN, &diag, &ns, T.data(), &ldR, tcol.data(), &inc1);
                    for (mwSize r = 0; r < nsel; ++r) T[r + static_cast<size_t>(nsel) * k] = tcol[r];
                }
                T[nsel + static_cast<size_t>(nsel) * k] = tau_i;

                const double beta = xcol[nsel];
                const double invdiag = (beta != 0.0) ? 1.0 / beta : 0.0;
                const double wcoeff = wn2blk[nr] + 1.0;   // pivot ||w||^2 + 1

                // Apply the new reflector to the rest of the block, then downdate their
                // w-vectors / running norms incrementally (BLAS-2 over nr contiguous cols).
                // The reflector touches only rows nsel..m-1 (zero above), so restrict there.
                if (nr > 0) {
                    const ptrdiff_t mlen = static_cast<ptrdiff_t>(m - nsel), nrr = static_cast<ptrdiff_t>(nr);
                    double *vtail = &V[static_cast<size_t>(nsel) + static_cast<size_t>(nsel) * m]; // [1; tail]
                    double *Xtail = &X[static_cast<size_t>(nsel)];   // row nsel; columns stride m
                    char tT = 'T';
                    dgemv(&tT, &mlen, &nrr, &one, Xtail, &ldA, vtail, &inc1, &zero, hw.data(), &inc1);
                    const double negtau = -tau_i;
                    dger(&mlen, &nrr, &negtau, vtail, &inc1, hw.data(), &inc1, Xtail, &ldA);
                    for (mwSize a = 0; a < nr; ++a)
                        betav[a] = Xtail[static_cast<size_t>(a) * m] * invdiag;   // alpha/diag
                    if (nsel > 0) {
                        double *wpiv = &Wblk[static_cast<size_t>(nr) * k];           // pivot w-vector
                        const ptrdiff_t ns = static_cast<ptrdiff_t>(nsel);
                        dgemv(&tT, &ns, &nrr, &one, Wblk.data(), &ldW, wpiv, &inc1, &zero, dots.data(), &inc1);
                        const double neg1 = -1.0;   // W(:,0:nr-1) <- W - wpivot beta_vec'
                        dger(&ns, &nrr, &neg1, wpiv, &inc1, betav.data(), &inc1, Wblk.data(), &ldW);
                    } else {
                        for (mwSize a = 0; a < nr; ++a) dots[a] = 0.0;
                    }
                    for (mwSize a = 0; a < nr; ++a) {
                        const double bv = betav[a];
                        Wblk[nsel + static_cast<size_t>(a) * k] = bv;            // new w-row
                        wn2blk[a] += -2.0 * bv * dots[a] + bv * bv * wcoeff;     // ||w||^2 downdate
                        const double alpha = X[nsel + static_cast<size_t>(a) * m];
                        sblk[a] = std::max(sblk[a] - alpha * alpha, 0.0);        // running-norm downdate

                        // Businger-Golub safeguard, mirroring the deterministic kernel
                        // (matlab/private/bsqr_mfile.m): once the running residual norm has
                        // decayed past norm_recomp_tol * (last exact value), recompute it
                        // exactly from the just-updated block column. Inline rather than
                        // deferred (as in the deterministic *panel*) because the in-block
                        // apply above has already updated X[:,a]; nsel has not yet been
                        // incremented, so the post-elimination residual is rows nsel+1..m-1.
                        // Residual norm only -- wn2blk is refreshed exactly at every block
                        // boundary, so it is not separately guarded.
                        if (sblk[a] <= sblk_ref[a] * opt.norm_recomp_tol) {
                            double s_exact = 0.0;
                            const double *xa = &X[static_cast<size_t>(a) * m];
                            for (mwSize r = nsel + 1; r < m; ++r) s_exact = std::fma(xa[r], xa[r], s_exact);
                            sblk[a] = sblk_ref[a] = s_exact;
                        }
                    }
                }

                // Commit the selection.
                f2 += cbest;
                selected[nsel] = idblk[nr];
                taken[idblk[nr]] = 1;
                since_last = 0;

                st_f2[nsel] = f2;
                st_crit[nsel] = cbest;
                st_thr[nsel] = theta;
                st_Fhat[nsel] = (static_cast<double>(nsel) + 1.0) *
                    (static_cast<double>(n) - nsel) / (static_cast<double>(k) - nsel);
                st_fallback[nsel] = 0.0;
                if (nsel == block_start) {        // attribute the block's apply to its first pick
                    st_samples[nsel] = static_cast<double>(bcount);
                    st_rounds[nsel] = 1.0;
                }

                nact = nr;   // pivot removed
                ++nsel;
            }

            // --- return the unselected drawn columns to the pool ---
            if (weighted) {
                for (mwSize id : drawn) {
                    if (!taken[id]) bit.add(id, g[id]);
                }
            } else {
                // Only the drawn columns (remaining[0..bcount-1]) can have been
                // taken this block, so swap-pop the taken ones out of the pool in
                // O(bcount) rather than rescanning the whole O(rcount) pool (which
                // would make a failed sampling block O(n) and the uniform-on-needle
                // cost O(n^2/k) instead of linear). Refill each hole from the tail
                // and re-examine it -- the tail can itself be a taken drawn column.
                mwSize t = 0;
                while (t < bcount && t < rcount) {
                    if (taken[remaining[t]]) {
                        remaining[t] = remaining[rcount - 1];
                        --rcount;
                    } else {
                        ++t;
                    }
                }
            }
        }
    } else {
    for (mwSize nsel = 0; nsel < k; ++nsel) {
        const double theta = threshold_value(opt.threshold_mode, f2, nsel, k, n) * opt.slack;
        const double Fhat_next = (static_cast<double>(nsel) + 1.0) *
            (static_cast<double>(n) - nsel) / (static_cast<double>(k) - nsel);
        st_thr[nsel] = theta;
        st_Fhat[nsel] = Fhat_next;

        if (weighted) drawn.clear();

        double best_c = std::numeric_limits<double>::infinity();
        mwSize best_id = 0, best_pos = 0;     // best_pos: position in `remaining` (uniform only)
        bool accepted = false;
        mwSize accept_id = 0, accept_c_local = 0, accept_pos = 0;
        double accept_c = std::numeric_limits<double>::infinity();
        mwSize tested = 0, rounds = 0;

        for (mwSize pos = 0; pos < rem_count && !accepted; pos += block) {
            const mwSize bcount = std::min(block, rem_count - pos);

            // Draw bcount candidates and gather them into X (m x bcount).
            for (mwSize t = 0; t < bcount; ++t) {
                mwSize id;
                if (!weighted) {
                    const mwSize p = pos + t;  // partial Fisher-Yates: freeze position p
                    std::uniform_int_distribution<mwSize> d(p, rem_count - 1);
                    std::swap(remaining[p], remaining[d(rng)]);
                    id = remaining[p];
                } else {
                    double total = bit.cur_total;
                    if (total <= 0.0) {  // degenerate: all remaining weights ~0
                        id = n;          // sentinel -> linear fallback below
                        for (mwSize j = 0; j < n; ++j) {
                            if (!taken[j]) { id = j; break; }
                        }
                    } else {
                        double x = unif01(rng) * total;
                        if (x >= total) x = std::nextafter(total, 0.0);
                        id = bit.find(x);
                    }
                    bit.add(id, -g[id]);   // remove for without-replacement draws
                    drawn.push_back(id);
                }
                ids[t] = id;
                const double *src = &A[static_cast<size_t>(id) * m];
                std::copy(src, src + m, &X[static_cast<size_t>(t) * m]);
            }

            // Apply the accumulated reflectors as a block: X <- Q_nsel' X
            //   = (I - V T' V') X  (compact-WY, BLAS-3).
            if (nsel > 0) {
                const ptrdiff_t ns = static_cast<ptrdiff_t>(nsel);
                const ptrdiff_t bc = static_cast<ptrdiff_t>(bcount);
                const ptrdiff_t mm = ldA;
                char tT = 'T', tN = 'N';
                dgemm(&tT, &tN, &ns, &bc, &mm, &one, V.data(), &ldA, X.data(), &ldA,
                      &zero, WB.data(), &ns);
                char side = 'L', uplo = 'U', tr = 'T', diag = 'N';
                dtrmm(&side, &uplo, &tr, &diag, &ns, &bc, &one, T.data(), &ldR, WB.data(), &ns);
                const double neg1 = -1.0;
                dgemm(&tN, &tN, &mm, &bc, &ns, &neg1, V.data(), &ldA, WB.data(), &ns,
                      &one, X.data(), &ldA);
            }

            // wn2 = column norms of R11^{-1} * (top nsel rows of X).
            if (nsel > 0) {
                for (mwSize t = 0; t < bcount; ++t) {
                    std::copy(&X[static_cast<size_t>(t) * m], &X[static_cast<size_t>(t) * m] + nsel,
                              &Xtop[static_cast<size_t>(t) * nsel]);
                }
                char side = 'L', uplo = 'U', tr = 'N', diag = 'N';
                ptrdiff_t ns = static_cast<ptrdiff_t>(nsel), bc = static_cast<ptrdiff_t>(bcount);
                ptrdiff_t ldb = static_cast<ptrdiff_t>(nsel);
                dtrsm(&side, &uplo, &tr, &diag, &ns, &bc, &one, R11.data(), &ldR, Xtop.data(), &ldb);
            }

            mwSize blk_best = 0;
            double blk_best_c = std::numeric_limits<double>::infinity();
            for (mwSize t = 0; t < bcount; ++t) {
                double rho2 = 0.0;
                const double *xc = &X[static_cast<size_t>(t) * m];
                for (mwSize r = nsel; r < m; ++r) rho2 = std::fma(xc[r], xc[r], rho2);
                double wn2 = 0.0;
                if (nsel > 0) {
                    const double *wc = &Xtop[static_cast<size_t>(t) * nsel];
                    for (mwSize r = 0; r < nsel; ++r) wn2 = std::fma(wc[r], wc[r], wn2);
                }
                const double c = (rho2 > 0.0) ? (1.0 + wn2) / rho2
                                              : std::numeric_limits<double>::infinity();
                cbuf[t] = c;
                if (c < blk_best_c) { blk_best_c = c; blk_best = t; }
            }

            tested += bcount;
            ++rounds;
            if (blk_best_c < best_c) {
                best_c = blk_best_c;
                best_id = ids[blk_best];
                best_pos = pos + blk_best;
            }

            if (opt.pick == Pick::First) {
                for (mwSize t = 0; t < bcount; ++t) {
                    if (cbuf[t] <= theta) {
                        accept_id = ids[t]; accept_c = cbuf[t];
                        accept_c_local = t; accept_pos = pos + t; accepted = true;
                        break;
                    }
                }
            } else if (blk_best_c <= theta) {
                accept_id = ids[blk_best]; accept_c = blk_best_c;
                accept_c_local = blk_best; accept_pos = pos + blk_best; accepted = true;
            }
        }

        bool fallback = false;
        if (!accepted) {
            accept_id = best_id; accept_c = best_c; accept_pos = best_pos; fallback = true;
        }

        // Rank guard: a non-finite increment means every remaining column had a
        // ~zero residual (rho^2 = 0) -- the input is rank-deficient for this k.
        // Fail loudly rather than propagate Inf into f2/R11 (the bound cannot be
        // maintained past the numerical rank; targets the full-rank GKS setting).
        if (!std::isfinite(accept_c)) {
            fail("bsqr_rand:RankDeficient",
                 "All remaining columns have ~zero residual; the input appears "
                 "rank-deficient for the requested k. Reduce k to the numerical rank.");
        }

        // Reduce the accepted column. If it was accepted from the last block its
        // reduced form is still in X (col accept_c_local); otherwise re-apply.
        if (accepted) {
            const double *xsrc = &X[static_cast<size_t>(accept_c_local) * m];
            std::copy(xsrc, xsrc + m, xcol.data());
        } else {
            std::copy(&A[static_cast<size_t>(accept_id) * m],
                      &A[static_cast<size_t>(accept_id) * m] + m, xcol.data());
            if (nsel > 0) {  // single-column compact-WY apply
                const ptrdiff_t ns = static_cast<ptrdiff_t>(nsel), mm = ldA;
                char tT = 'T', tN = 'N';
                dgemv(&tT, &mm, &ns, &one, V.data(), &ldA, xcol.data(), &inc1, &zero, tcol.data(), &inc1);
                char uplo = 'U', diag = 'N';
                dtrmv(&uplo, &tT, &diag, &ns, T.data(), &ldR, tcol.data(), &inc1);
                const double neg1 = -1.0;
                dgemv(&tN, &mm, &ns, &neg1, V.data(), &ldA, tcol.data(), &inc1, &one, xcol.data(), &inc1);
            }
        }

        // R11 column (top nsel) and the new reflector / diagonal.
        for (mwSize r = 0; r < nsel; ++r) R11[r + static_cast<size_t>(nsel) * k] = xcol[r];
        double tau_i = 0.0;
        ptrdiff_t len = static_cast<ptrdiff_t>(m - nsel);
        if (len > 1) dlarfg(&len, &xcol[nsel], &xcol[nsel + 1], &inc1, &tau_i);
        R11[nsel + static_cast<size_t>(nsel) * k] = xcol[nsel];
        V[nsel + static_cast<size_t>(nsel) * m] = 1.0;
        for (mwSize r = nsel + 1; r < m; ++r) V[r + static_cast<size_t>(nsel) * m] = xcol[r];
        tau[nsel] = tau_i;

        // Incremental compact-WY update: new column of T (dlarft, forward).
        if (nsel > 0) {
            const ptrdiff_t ns = static_cast<ptrdiff_t>(nsel), mm = ldA;
            char tT = 'T';
            dgemv(&tT, &mm, &ns, &one, V.data(), &ldA, &V[static_cast<size_t>(nsel) * m], &inc1,
                  &zero, tcol.data(), &inc1);
            const double negtau = -tau_i;
            for (mwSize r = 0; r < nsel; ++r) tcol[r] *= negtau;
            char uplo = 'U', tN = 'N', diag = 'N';
            dtrmv(&uplo, &tN, &diag, &ns, T.data(), &ldR, tcol.data(), &inc1);
            for (mwSize r = 0; r < nsel; ++r) T[r + static_cast<size_t>(nsel) * k] = tcol[r];
        }
        T[nsel + static_cast<size_t>(nsel) * k] = tau_i;

        f2 += accept_c;
        selected[nsel] = accept_id;
        taken[accept_id] = 1;

        // Remove the pivot from the pool; restore the other drawn weights.
        if (!weighted) {
            remaining[accept_pos] = remaining[rem_count - 1];  // O(1) swap-pop
        } else {
            for (mwSize id : drawn) {
                if (id != accept_id) bit.add(id, g[id]);  // restore non-selected draws
            }
            // accept_id stays removed (taken).
        }
        --rem_count;

        st_f2[nsel] = f2;
        st_crit[nsel] = accept_c;
        st_samples[nsel] = static_cast<double>(tested);
        st_rounds[nsel] = static_cast<double>(rounds);
        st_fallback[nsel] = fallback ? 1.0 : 0.0;
    }
    }

    // --- Outputs -----------------------------------------------------------
    if (nlhs == 0) return;

    mxArray *pArr = mxCreateDoubleMatrix(1, n, mxREAL);
    double *pp = mxGetPr(pArr);
    for (mwSize j = 0; j < k; ++j) pp[j] = static_cast<double>(selected[j] + 1);
    mwSize w = k;
    for (mwSize j = 0; j < n; ++j) {
        if (!taken[j]) pp[w++] = static_cast<double>(j + 1);
    }
    plhs[0] = pArr;
    if (nlhs == 1) return;

    // Q (economy, m x k), formed lazily -- only when this output is requested.
    // V is already in geqrf packed layout (reflector tails below the diagonal;
    // dorgqr ignores the diagonal and above), so accumulate straight into the
    // output buffer with one dorgqr call, O(m k^2).
    mxArray *Qm = mxCreateDoubleMatrix(m, k, mxREAL);
    if (k > 0) {
        double *q = mxGetPr(Qm);
        std::copy(V.begin(), V.end(), q);
        ptrdiff_t mm = static_cast<ptrdiff_t>(m), kk = static_cast<ptrdiff_t>(k), info = 0;
        double wq = 0.0;
        ptrdiff_t lwork = -1;
        dorgqr(&mm, &kk, &kk, q, &mm, tau.data(), &wq, &lwork, &info);  // workspace query
        lwork = std::max<ptrdiff_t>(static_cast<ptrdiff_t>(wq), kk);
        std::vector<double> work(static_cast<size_t>(lwork));
        dorgqr(&mm, &kk, &kk, q, &mm, tau.data(), work.data(), &lwork, &info);
        if (info != 0) {
            fail("bsqr_rand:LapackError", "dorgqr failed to form Q.");
        }
    }
    plhs[1] = Qm;
    if (nlhs == 2) return;

    mxArray *R11m = mxCreateDoubleMatrix(k, k, mxREAL);
    double *rptr = mxGetPr(R11m);
    for (mwSize j = 0; j < k; ++j) {
        for (mwSize i = 0; i < k; ++i) {
            rptr[i + static_cast<size_t>(j) * k] = (i <= j) ? R11[i + static_cast<size_t>(j) * k] : 0.0;
        }
    }
    plhs[2] = R11m;
    if (nlhs == 3) return;

    const char *sfields[] = {"f2", "crit", "threshold", "Fhat", "samples_tested",
                             "rounds", "fallback", "frob_inv", "osinsky_bound", "total_tested",
                             "blocks_sampled"};
    mxArray *stats = mxCreateStructMatrix(1, 1, 11, sfields);
    auto set_row = [&](const char *name, const std::vector<double> &v) {
        mxArray *a = mxCreateDoubleMatrix(1, k, mxREAL);
        std::copy(v.begin(), v.end(), mxGetPr(a));
        mxSetField(stats, 0, name, a);
    };
    set_row("f2", st_f2);
    set_row("crit", st_crit);
    set_row("threshold", st_thr);
    set_row("Fhat", st_Fhat);
    set_row("samples_tested", st_samples);
    set_row("rounds", st_rounds);
    set_row("fallback", st_fallback);
    double total_tested = std::accumulate(st_samples.begin(), st_samples.end(), 0.0);
    mxSetField(stats, 0, "frob_inv", mxCreateDoubleScalar(std::sqrt(std::max(f2, 0.0))));
    mxSetField(stats, 0, "osinsky_bound",
               mxCreateDoubleScalar(std::sqrt(static_cast<double>(k) * (n - k + 1))));
    mxSetField(stats, 0, "total_tested", mxCreateDoubleScalar(total_tested));
    double blocks_sampled = std::accumulate(st_rounds.begin(), st_rounds.end(), 0.0);
    mxSetField(stats, 0, "blocks_sampled", mxCreateDoubleScalar(blocks_sampled));
    plhs[3] = stats;
    if (nlhs == 4) return;

    // R12 (optional): Q(:,1:k)' * A(:,leftover), top k rows, via compact-WY.
    mwSize n12 = (k < n) ? (n - k) : 0;
    mxArray *R12m = mxCreateDoubleMatrix(k, n12, mxREAL);
    if (opt.return_r12 && k > 0 && n12 > 0) {
        double *r12 = mxGetPr(R12m);
        std::vector<double> Xr(static_cast<size_t>(m) * n12);
        mwSize col = 0;
        for (mwSize j = 0; j < n; ++j) {
            if (!taken[j]) {
                const double *src = &A[static_cast<size_t>(j) * m];
                std::copy(src, src + m, &Xr[static_cast<size_t>(col++) * m]);
            }
        }
        std::vector<double> WBr(static_cast<size_t>(k) * n12);
        const ptrdiff_t kk = static_cast<ptrdiff_t>(k), nn = static_cast<ptrdiff_t>(n12), mm = ldA;
        char tT = 'T', tN = 'N';
        dgemm(&tT, &tN, &kk, &nn, &mm, &one, V.data(), &ldA, Xr.data(), &ldA, &zero, WBr.data(), &kk);
        char side = 'L', uplo = 'U', tr = 'T', diag = 'N';
        dtrmm(&side, &uplo, &tr, &diag, &kk, &nn, &one, T.data(), &ldR, WBr.data(), &kk);
        const double neg1 = -1.0;
        dgemm(&tN, &tN, &mm, &nn, &kk, &neg1, V.data(), &ldA, WBr.data(), &kk, &one, Xr.data(), &ldA);
        for (mwSize j = 0; j < n12; ++j) {
            std::copy(&Xr[static_cast<size_t>(j) * m], &Xr[static_cast<size_t>(j) * m] + k,
                      &r12[static_cast<size_t>(j) * k]);
        }
    }
    plhs[4] = R12m;
}
