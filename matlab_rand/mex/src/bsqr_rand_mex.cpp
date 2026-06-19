// Randomized Bischof-Stewart column selection -- MEX backend.
//
// Mirrors matlab_rand/private/bsqr_rand_mfile.m. The kernel never maintains
// R11^{-1}R12 or column norms for every column; it tracks only the running
// squared inverse Frobenius norm f2 = ||R11^{-1}||_F^2 and, per step, samples
// candidate columns in blocks, brings each into the current reduced frame with
// the accumulated reflectors (BLAS-2/3), and accepts the first/best one whose
// increment keeps f2 under the per-step threshold.
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
    mwSize block_size = 16;
    ThresholdMode threshold_mode = ThresholdMode::RunningMean;
    double slack = 1.0;
    Sampling sampling = Sampling::Uniform;
    Pick pick = Pick::BestInBlock;
    bool has_seed = false;
    unsigned long long seed = 0;
    bool return_r12 = false;
    bool check_finite = true;
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
            // Accepted for API compatibility; ignored inside the MEX.
            const std::string b = to_lower(get_string(val, "bsqr_rand:InvalidBackend"));
            if (!(b == "auto" || b == "mfile" || b == "mex")) {
                fail("bsqr_rand:InvalidBackend", "backend must be \"auto\", \"mfile\", or \"mex\".");
            }
        } else if (name == "check_finite") {
            opt.check_finite = scalar_to_bool(val, "bsqr_rand:InvalidCheckFinite");
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
    const mwSize block = std::min(opt.block_size, std::max<mwSize>(n, 1));

    // --- State -------------------------------------------------------------
    std::vector<double> V(static_cast<size_t>(m) * k, 0.0);       // reflector store (unit diag)
    std::vector<double> R11(static_cast<size_t>(k) * k, 0.0);     // k x k, lda = k
    std::vector<double> tau(k, 0.0);
    std::vector<mwSize> selected(k, 0);
    std::vector<mwSize> remaining(n);
    std::iota(remaining.begin(), remaining.end(), static_cast<mwSize>(0));
    double f2 = 0.0;

    // --- Stats -------------------------------------------------------------
    std::vector<double> st_f2(k, 0.0), st_crit(k, 0.0), st_thr(k, 0.0), st_Fhat(k, 0.0);
    std::vector<double> st_samples(k, 0.0), st_rounds(k, 0.0), st_fallback(k, 0.0);

    // --- Scratch -----------------------------------------------------------
    std::vector<double> X(static_cast<size_t>(m) * std::max<mwSize>(block, 1), 0.0);
    std::vector<double> Xtop(static_cast<size_t>(std::max<mwSize>(k, 1)) * std::max<mwSize>(block, 1), 0.0);
    std::vector<double> wbuf(std::max<mwSize>(block, 1), 0.0);
    std::vector<double> cbuf(std::max<mwSize>(block, 1), 0.0);
    std::vector<double> xcol(m, 0.0);
    std::vector<mwSize> visit(n, 0);
    std::vector<double> keys;
    std::vector<mwSize> g_order;
    std::vector<double> g;  // starting squared column norms (normweighted only)

    if (opt.sampling == Sampling::NormWeighted) {
        g.assign(n, 0.0);
        ptrdiff_t mm = static_cast<ptrdiff_t>(m);
        ptrdiff_t inc1 = 1;
        for (mwSize j = 0; j < n; ++j) {
            const double nrm = dnrm2(&mm, const_cast<double *>(&A[static_cast<size_t>(j) * m]), &inc1);
            g[j] = nrm * nrm;
        }
    }

    std::mt19937_64 rng;
    if (opt.has_seed) {
        rng.seed(opt.seed);
    } else {
        std::random_device rd;
        rng.seed((static_cast<unsigned long long>(rd()) << 32) ^ rd());
    }
    std::uniform_real_distribution<double> unif01(std::numeric_limits<double>::min(), 1.0);

    const double one = 1.0, zero = 0.0;
    ptrdiff_t inc1 = 1;
    const ptrdiff_t ldA = static_cast<ptrdiff_t>(m);
    const ptrdiff_t ldR = static_cast<ptrdiff_t>(k);

    for (mwSize nsel = 0; nsel < k; ++nsel) {
        const mwSize rem_count = n - nsel;
        const double theta =
            threshold_value(opt.threshold_mode, f2, nsel, k, n) * opt.slack;
        const double Fhat_next =
            (static_cast<double>(nsel) + 1.0) * (static_cast<double>(n) - nsel) /
            (static_cast<double>(k) - nsel);
        st_thr[nsel] = theta;
        st_Fhat[nsel] = Fhat_next;

        // Visiting order over the remaining columns.
        for (mwSize t = 0; t < rem_count; ++t) {
            visit[t] = remaining[t];
        }
        if (opt.sampling == Sampling::Uniform) {
            for (mwSize i = rem_count; i > 1; --i) {  // Fisher-Yates
                std::uniform_int_distribution<mwSize> d(0, i - 1);
                std::swap(visit[i - 1], visit[d(rng)]);
            }
        } else {
            keys.assign(rem_count, 0.0);
            g_order.assign(rem_count, 0);
            for (mwSize t = 0; t < rem_count; ++t) {
                const double w = std::max(g[remaining[t]], std::numeric_limits<double>::min());
                keys[t] = -std::log(unif01(rng)) / w;  // Efraimidis-Spirakis
                g_order[t] = t;
            }
            std::sort(g_order.begin(), g_order.end(),
                      [&keys](mwSize a, mwSize b) { return keys[a] < keys[b]; });
            for (mwSize t = 0; t < rem_count; ++t) {
                visit[t] = remaining[g_order[t]];
            }
        }

        double best_c = std::numeric_limits<double>::infinity();
        mwSize best_id = visit[0];
        bool accepted = false;
        mwSize accept_id = 0;
        double accept_c = std::numeric_limits<double>::infinity();
        mwSize accept_local = 0;   // column index of the accepted pivot within the last block
        mwSize tested = 0, rounds = 0;

        for (mwSize pos = 0; pos < rem_count && !accepted; pos += block) {
            const mwSize bcount = std::min(block, rem_count - pos);

            // Gather candidate columns into X (m x bcount).
            for (mwSize t = 0; t < bcount; ++t) {
                const double *src = &A[static_cast<size_t>(visit[pos + t]) * m];
                std::copy(src, src + m, &X[static_cast<size_t>(t) * m]);
            }

            // Apply the nsel accumulated reflectors: X <- Q_nsel * X.
            for (mwSize s = 0; s < nsel; ++s) {
                if (tau[s] == 0.0) {
                    continue;
                }
                const char trT = 'T';
                ptrdiff_t mm = ldA, nn = static_cast<ptrdiff_t>(bcount);
                dgemv(&trT, &mm, &nn, &one, X.data(), &ldA, &V[static_cast<size_t>(s) * m],
                      &inc1, &zero, wbuf.data(), &inc1);
                const double negtau = -tau[s];
                dger(&mm, &nn, &negtau, &V[static_cast<size_t>(s) * m], &inc1, wbuf.data(),
                     &inc1, X.data(), &ldA);
            }

            // wn2 = column norms of R11^{-1} * (top nsel rows of X).
            if (nsel > 0) {
                for (mwSize t = 0; t < bcount; ++t) {
                    std::copy(&X[static_cast<size_t>(t) * m],
                              &X[static_cast<size_t>(t) * m] + nsel,
                              &Xtop[static_cast<size_t>(t) * nsel]);
                }
                const char side = 'L', uplo = 'U', tr = 'N', diag = 'N';
                ptrdiff_t mm = static_cast<ptrdiff_t>(nsel), nn = static_cast<ptrdiff_t>(bcount);
                ptrdiff_t ldb = static_cast<ptrdiff_t>(nsel);
                dtrsm(&side, &uplo, &tr, &diag, &mm, &nn, &one, R11.data(), &ldR, Xtop.data(), &ldb);
            }

            mwSize blk_best = 0;
            double blk_best_c = std::numeric_limits<double>::infinity();
            for (mwSize t = 0; t < bcount; ++t) {
                double rho2 = 0.0;
                const double *xc = &X[static_cast<size_t>(t) * m];
                for (mwSize r = nsel; r < m; ++r) {
                    rho2 = std::fma(xc[r], xc[r], rho2);
                }
                double wn2 = 0.0;
                if (nsel > 0) {
                    const double *wc = &Xtop[static_cast<size_t>(t) * nsel];
                    for (mwSize r = 0; r < nsel; ++r) {
                        wn2 = std::fma(wc[r], wc[r], wn2);
                    }
                }
                const double c = (rho2 > 0.0) ? (1.0 + wn2) / rho2
                                              : std::numeric_limits<double>::infinity();
                cbuf[t] = c;
                if (c < blk_best_c) {
                    blk_best_c = c;
                    blk_best = t;
                }
            }

            tested += bcount;
            ++rounds;
            if (blk_best_c < best_c) {
                best_c = blk_best_c;
                best_id = visit[pos + blk_best];
            }

            if (opt.pick == Pick::First) {
                for (mwSize t = 0; t < bcount; ++t) {
                    if (cbuf[t] <= theta) {
                        accept_id = visit[pos + t];
                        accept_c = cbuf[t];
                        accept_local = t;
                        accepted = true;
                        break;
                    }
                }
            } else if (blk_best_c <= theta) {
                accept_id = visit[pos + blk_best];
                accept_c = blk_best_c;
                accept_local = blk_best;
                accepted = true;
            }
        }

        bool fallback = false;
        if (!accepted) {
            accept_id = best_id;
            accept_c = best_c;
            fallback = true;
        }

        // Reduce the accepted column and append its reflector / R11 column. When
        // the pivot was accepted from a sampled block its reduced form already
        // sits in X (the reflectors were applied there) -- reuse it. Only the
        // exhaustive-fallback case (X overwritten) needs a fresh apply.
        if (accepted) {
            const double *xsrc = &X[static_cast<size_t>(accept_local) * m];
            std::copy(xsrc, xsrc + m, xcol.data());
        } else {
            std::copy(&A[static_cast<size_t>(accept_id) * m],
                      &A[static_cast<size_t>(accept_id) * m] + m, xcol.data());
            for (mwSize s = 0; s < nsel; ++s) {
                if (tau[s] == 0.0) {
                    continue;
                }
                ptrdiff_t mm = ldA;
                const double d = ddot(&mm, &V[static_cast<size_t>(s) * m], &inc1, xcol.data(), &inc1);
                const double coef = -tau[s] * d;
                daxpy(&mm, &coef, &V[static_cast<size_t>(s) * m], &inc1, xcol.data(), &inc1);
            }
        }
        for (mwSize r = 0; r < nsel; ++r) {
            R11[r + static_cast<size_t>(nsel) * k] = xcol[r];
        }
        double tau_i = 0.0;
        ptrdiff_t len = static_cast<ptrdiff_t>(m - nsel);
        if (len > 1) {
            dlarfg(&len, &xcol[nsel], &xcol[nsel + 1], &inc1, &tau_i);
        }
        const double beta_i = xcol[nsel];
        R11[nsel + static_cast<size_t>(nsel) * k] = beta_i;
        V[nsel + static_cast<size_t>(nsel) * m] = 1.0;
        for (mwSize r = nsel + 1; r < m; ++r) {
            V[r + static_cast<size_t>(nsel) * m] = xcol[r];
        }
        tau[nsel] = tau_i;

        f2 += accept_c;
        selected[nsel] = accept_id;
        // Remove accept_id from remaining (order irrelevant; reshuffled each step).
        for (mwSize t = 0; t < rem_count; ++t) {
            if (remaining[t] == accept_id) {
                remaining[t] = remaining[rem_count - 1];
                remaining.pop_back();
                break;
            }
        }

        st_f2[nsel] = f2;
        st_crit[nsel] = accept_c;
        st_samples[nsel] = static_cast<double>(tested);
        st_rounds[nsel] = static_cast<double>(rounds);
        st_fallback[nsel] = fallback ? 1.0 : 0.0;
    }

    // --- Outputs -----------------------------------------------------------
    if (nlhs == 0) {
        return;
    }

    // p (1 x n): selected first, then the leftover columns.
    mxArray *pArr = mxCreateDoubleMatrix(1, n, mxREAL);
    double *pp = mxGetPr(pArr);
    for (mwSize j = 0; j < k; ++j) {
        pp[j] = static_cast<double>(selected[j] + 1);
    }
    for (mwSize j = 0; j < remaining.size(); ++j) {
        pp[k + j] = static_cast<double>(remaining[j] + 1);
    }
    plhs[0] = pArr;
    if (nlhs == 1) {
        return;
    }

    // reflectors struct.
    const char *rfields[] = {"V", "tau", "m", "k"};
    mxArray *refl = mxCreateStructMatrix(1, 1, 4, rfields);
    mxArray *Vm = mxCreateDoubleMatrix(m, k, mxREAL);
    std::copy(V.begin(), V.end(), mxGetPr(Vm));
    mxSetField(refl, 0, "V", Vm);
    mxArray *taum = mxCreateDoubleMatrix(k, 1, mxREAL);
    std::copy(tau.begin(), tau.end(), mxGetPr(taum));
    mxSetField(refl, 0, "tau", taum);
    mxSetField(refl, 0, "m", mxCreateDoubleScalar(static_cast<double>(m)));
    mxSetField(refl, 0, "k", mxCreateDoubleScalar(static_cast<double>(k)));
    plhs[1] = refl;
    if (nlhs == 2) {
        return;
    }

    mxArray *R11m = mxCreateDoubleMatrix(k, k, mxREAL);
    double *rptr = mxGetPr(R11m);
    for (mwSize j = 0; j < k; ++j) {
        for (mwSize i = 0; i < k; ++i) {
            rptr[i + static_cast<size_t>(j) * k] = (i <= j) ? R11[i + static_cast<size_t>(j) * k] : 0.0;
        }
    }
    plhs[2] = R11m;
    if (nlhs == 3) {
        return;
    }

    // stats struct.
    const char *sfields[] = {"f2", "crit", "threshold", "Fhat", "samples_tested",
                             "rounds", "fallback", "frob_inv", "osinsky_bound", "total_tested"};
    mxArray *stats = mxCreateStructMatrix(1, 1, 10, sfields);
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
    plhs[3] = stats;
    if (nlhs == 4) {
        return;
    }

    // R12 (optional): Q(:,1:k)' * A(:,remaining), top k rows.
    mwSize n12 = (k < n) ? (n - k) : 0;
    mxArray *R12m = mxCreateDoubleMatrix(k, n12, mxREAL);
    if (opt.return_r12 && k > 0 && n12 > 0) {
        double *r12 = mxGetPr(R12m);
        std::vector<double> Xr(static_cast<size_t>(m) * n12);
        for (mwSize j = 0; j < n12; ++j) {
            const double *src = &A[static_cast<size_t>(remaining[j]) * m];
            std::copy(src, src + m, &Xr[static_cast<size_t>(j) * m]);
        }
        for (mwSize s = 0; s < k; ++s) {
            if (tau[s] == 0.0) {
                continue;
            }
            const char trT = 'T';
            ptrdiff_t mm = ldA, nn = static_cast<ptrdiff_t>(n12);
            std::vector<double> w(n12);
            dgemv(&trT, &mm, &nn, &one, Xr.data(), &ldA, &V[static_cast<size_t>(s) * m], &inc1,
                  &zero, w.data(), &inc1);
            const double negtau = -tau[s];
            dger(&mm, &nn, &negtau, &V[static_cast<size_t>(s) * m], &inc1, w.data(), &inc1,
                 Xr.data(), &ldA);
        }
        for (mwSize j = 0; j < n12; ++j) {
            std::copy(&Xr[static_cast<size_t>(j) * m], &Xr[static_cast<size_t>(j) * m] + k,
                      &r12[static_cast<size_t>(j) * k]);
        }
    }
    plhs[4] = R12m;
}
