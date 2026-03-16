#include "mex.h"
#include "blas.h"
#include "lapack.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace {

struct Options {
    mwSize k = 0;
    bool return_rinv_r12 = false;
    bool pivot_vector = false;  // false => matrix
    double norm_recomp_tol = std::sqrt(std::numeric_limits<double>::epsilon());
    bool check_finite = true;
};

struct Workspace {
    std::vector<double> A;
    std::vector<double> tau;
    std::vector<double> W;
    std::vector<double> wnorm2;
    std::vector<double> s;
    std::vector<double> s_ref;
    std::vector<mwSize> p;
    std::vector<double> beta_vec;
    std::vector<double> dots;
    std::vector<double> hh_work;
    std::vector<double> q_work;
};

Workspace &workspace() {
    static Workspace ws;
    return ws;
}

inline double &Aat(std::vector<double> &A, mwSize m, mwSize r, mwSize c) {
    return A[r + c * m];
}

inline const double &AatConst(const std::vector<double> &A, mwSize m, mwSize r, mwSize c) {
    return A[r + c * m];
}

inline double &Wat(std::vector<double> &W, mwSize k, mwSize r, mwSize c) {
    return W[r + c * k];
}

inline const double &WatConst(const std::vector<double> &W, mwSize k, mwSize r, mwSize c) {
    return W[r + c * k];
}

void fail(const char *id, const char *msg) {
    mexErrMsgIdAndTxt(id, "%s", msg);
}

bool scalar_to_bool(const mxArray *a, const char *id) {
    if (!(mxIsLogicalScalar(a) || (mxIsNumeric(a) && mxIsScalar(a)))) {
        mexErrMsgIdAndTxt(id, "Expected a logical/numeric scalar.");
    }
    const double v = mxIsLogical(a) ? (mxIsLogicalScalarTrue(a) ? 1.0 : 0.0) : mxGetScalar(a);
    return v != 0.0;
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

void parse_options(const mxArray *A, int nrhs, const mxArray *prhs[], Options &opt) {
    const mwSize m = mxGetM(A);
    const mwSize n = mxGetN(A);
    const mwSize kmax = std::min(m, n);
    opt.k = kmax;

    if (((nrhs - 1) % 2) != 0) {
        fail("bsqr:InvalidOptions", "Name-value options must be provided in pairs.");
    }

    for (int i = 1; i < nrhs; i += 2) {
        const std::string name = to_lower(get_string(prhs[i], "bsqr:InvalidOptionName"));
        const mxArray *val = prhs[i + 1];

        if (name == "k") {
            if (!(mxIsNumeric(val) && mxIsScalar(val) && !mxIsComplex(val))) {
                fail("bsqr:InvalidK", "k must be a real numeric scalar.");
            }
            const double kd = mxGetScalar(val);
            if (!std::isfinite(kd) || std::floor(kd) != kd) {
                fail("bsqr:InvalidK", "k must be an integer in [0, min(size(A))].");
            }
            if (kd < 0.0 || kd > static_cast<double>(kmax)) {
                fail("bsqr:InvalidK", "k must satisfy 0 <= k <= min(size(A)).");
            }
            opt.k = static_cast<mwSize>(kd);
        } else if (name == "return_rinv_r12") {
            opt.return_rinv_r12 = scalar_to_bool(val, "bsqr:InvalidReturnRinvR12");
        } else if (name == "pivot_format") {
            const std::string pf = to_lower(get_string(val, "bsqr:InvalidPivotFormat"));
            if (pf == "vector") {
                opt.pivot_vector = true;
            } else if (pf == "matrix") {
                opt.pivot_vector = false;
            } else {
                fail("bsqr:InvalidPivotFormat", "pivot_format must be \"matrix\" or \"vector\".");
            }
        } else if (name == "backend") {
            // Accepted for API compatibility; ignored inside MEX.
            const std::string b = to_lower(get_string(val, "bsqr:InvalidBackend"));
            if (!(b == "auto" || b == "mfile" || b == "mex")) {
                fail("bsqr:InvalidBackend", "backend must be \"auto\", \"mfile\", or \"mex\".");
            }
        } else if (name == "norm_recomp_tol") {
            if (!(mxIsNumeric(val) && mxIsScalar(val) && !mxIsComplex(val))) {
                fail("bsqr:InvalidNormRecompTol", "norm_recomp_tol must be a real numeric scalar.");
            }
            const double t = mxGetScalar(val);
            if (!std::isfinite(t) || t < 0.0 || t > 1.0) {
                fail("bsqr:InvalidNormRecompTol", "norm_recomp_tol must satisfy 0 <= value <= 1.");
            }
            opt.norm_recomp_tol = t;
        } else if (name == "check_finite") {
            opt.check_finite = scalar_to_bool(val, "bsqr:InvalidCheckFinite");
        } else {
            fail("bsqr:UnknownOption", "Unknown bsqr option.");
        }
    }
}

void householder_column(std::vector<double> &A, mwSize m, mwSize i, double &tau, double &beta) {
    const ptrdiff_t n = static_cast<ptrdiff_t>(m - i);
    if (n <= 1) {
        tau = 0.0;
        beta = AatConst(A, m, i, i);
        return;
    }
    ptrdiff_t inc1 = 1;
    double *alpha = &Aat(A, m, i, i);
    double *x = &Aat(A, m, i + 1, i);
    // Use LAPACK's stable reflector constructor.
    dlarfg(&n, alpha, x, &inc1, &tau);
    beta = *alpha;
}

void apply_householder_left(
    std::vector<double> &A,
    mwSize m,
    mwSize n,
    mwSize i,
    double tau,
    std::vector<double> &workbuf
) {
    const ptrdiff_t rows = static_cast<ptrdiff_t>(m - i);
    const ptrdiff_t cols_total = static_cast<ptrdiff_t>(n - i - 1);
    if (tau == 0.0 || cols_total <= 0) {
        return;
    }

    const char trans = 'T';
    ptrdiff_t inc1 = 1;
    ptrdiff_t ldc = static_cast<ptrdiff_t>(m);
    // Use the full trailing width in one BLAS call. Small fixed blocks add
    // call overhead and can suppress vendor BLAS heuristics on larger sizes.
    const ptrdiff_t block_cols = cols_total;
    const size_t need = static_cast<size_t>(std::max<ptrdiff_t>(1, block_cols));
    if (workbuf.size() < need) {
        workbuf.resize(need);
    }

    double *v = &Aat(A, m, i, i);
    const double one = 1.0;
    const double zero = 0.0;
    const double alpha = -tau;

    for (ptrdiff_t off = 0; off < cols_total; off += block_cols) {
        const ptrdiff_t cols = std::min(block_cols, cols_total - off);
        double *C = &Aat(A, m, i, i + 1 + static_cast<mwSize>(off));
        // work := C' * v
        dgemv(&trans, &rows, &cols, &one, C, &ldc, v, &inc1, &zero, workbuf.data(), &inc1);
        // C := C - tau * v * work'
        dger(&rows, &cols, &alpha, v, &inc1, workbuf.data(), &inc1, C, &ldc);
    }
}

void build_q(
    const std::vector<double> &A,
    mwSize m,
    mwSize k,
    const std::vector<double> &tau,
    double *Q,
    std::vector<double> &q_work
) {
    std::fill(Q, Q + static_cast<std::ptrdiff_t>(m) * static_cast<std::ptrdiff_t>(k), 0.0);
    for (mwSize c = 0; c < k; ++c) {
        Q[c + c * m] = 1.0;
    }

    if (k == 0) {
        return;
    }

    const char side = 'L';
    const char trans = 'N';
    const ptrdiff_t mm = static_cast<ptrdiff_t>(m);
    const ptrdiff_t nn = static_cast<ptrdiff_t>(k);
    const ptrdiff_t kk = static_cast<ptrdiff_t>(k);
    const ptrdiff_t lda = static_cast<ptrdiff_t>(m);
    const ptrdiff_t ldc = static_cast<ptrdiff_t>(m);

    ptrdiff_t info = 0;
    ptrdiff_t lwork = -1;
    double work_query = 0.0;
    dormqr(&side, &trans, &mm, &nn, &kk, A.data(), &lda, tau.data(), Q, &ldc, &work_query, &lwork, &info);
    if (info != 0) {
        fail("bsqr:QBuildFailed", "LAPACK dormqr workspace query failed.");
    }

    lwork = std::max<ptrdiff_t>(1, static_cast<ptrdiff_t>(work_query));
    if (q_work.size() < static_cast<size_t>(lwork)) {
        q_work.resize(static_cast<size_t>(lwork));
    }
    dormqr(&side, &trans, &mm, &nn, &kk, A.data(), &lda, tau.data(), Q, &ldc, q_work.data(), &lwork, &info);
    if (info != 0) {
        fail("bsqr:QBuildFailed", "LAPACK dormqr failed while forming Q.");
    }
}

mxArray *make_R(const std::vector<double> &A, mwSize m, mwSize n, mwSize k) {
    mxArray *R = mxCreateDoubleMatrix(k, n, mxREAL);
    double *rptr = mxGetPr(R);
    for (mwSize j = 0; j < n; ++j) {
        for (mwSize i = 0; i < k; ++i) {
            rptr[i + j * k] = (i <= j) ? AatConst(A, m, i, j) : 0.0;
        }
    }
    return R;
}

mxArray *make_Q(
    const std::vector<double> &A,
    mwSize m,
    mwSize k,
    const std::vector<double> &tau,
    std::vector<double> &q_work
) {
    mxArray *Qm = mxCreateDoubleMatrix(m, k, mxREAL);
    double *qptr = mxGetPr(Qm);
    build_q(A, m, k, tau, qptr, q_work);
    return Qm;
}

mxArray *make_pivot_vector(const std::vector<mwSize> &p) {
    const mwSize n = p.size();
    mxArray *pv = mxCreateDoubleMatrix(1, n, mxREAL);
    double *ptr = mxGetPr(pv);
    for (mwSize j = 0; j < n; ++j) {
        ptr[j] = static_cast<double>(p[j] + 1);
    }
    return pv;
}

mxArray *make_pivot_matrix(const std::vector<mwSize> &p) {
    const mwSize n = p.size();
    mxArray *E = mxCreateDoubleMatrix(n, n, mxREAL);
    double *eptr = mxGetPr(E);
    std::fill(eptr, eptr + static_cast<size_t>(n) * static_cast<size_t>(n), 0.0);
    for (mwSize j = 0; j < n; ++j) {
        eptr[p[j] + j * n] = 1.0;
    }
    return E;
}

mxArray *make_rinv(const std::vector<double> &W, mwSize k, mwSize n, bool want) {
    if (!want) {
        return mxCreateDoubleMatrix(0, 0, mxREAL);
    }
    if (k == 0) {
        return mxCreateDoubleMatrix(0, n, mxREAL);
    }
    if (k >= n) {
        return mxCreateDoubleMatrix(k, 0, mxREAL);
    }
    const mwSize n12 = n - k;
    mxArray *out = mxCreateDoubleMatrix(k, n12, mxREAL);
    double *ptr = mxGetPr(out);
    for (mwSize j = 0; j < n12; ++j) {
        for (mwSize i = 0; i < k; ++i) {
            ptr[i + j * k] = WatConst(W, k, i, k + j);
        }
    }
    return out;
}

mwSize select_pivot_column(
    const std::vector<double> &s,
    const std::vector<double> &wnorm2,
    mwSize i,
    mwSize n
) {
    mwSize best_j = i;
    double best_c = std::numeric_limits<double>::infinity();
    for (mwSize j = i; j < n; ++j) {
        const double sj = s[j];
        // Bischof-Stewart pivot score: minimize (1 + ||w_j||^2) / ||a_j^(i)||^2.
        const double cj = (sj > 0.0) ? (1.0 + wnorm2[j]) / sj : std::numeric_limits<double>::infinity();
        if (cj < best_c) {
            best_c = cj;
            best_j = j;
        }
    }
    return best_j;
}

void swap_pivot_state(
    std::vector<double> &A,
    std::vector<double> &W,
    std::vector<double> &s,
    std::vector<double> &s_ref,
    std::vector<double> &wnorm2,
    std::vector<mwSize> &p,
    mwSize m,
    mwSize k,
    mwSize i,
    mwSize best_j
) {
    if (best_j == i) {
        return;
    }

    ptrdiff_t inc1 = 1;
    const ptrdiff_t mptr = static_cast<ptrdiff_t>(m);
    dswap(&mptr, &Aat(A, m, 0, i), &inc1, &Aat(A, m, 0, best_j), &inc1);
    if (i > 0) {
        const ptrdiff_t iptr = static_cast<ptrdiff_t>(i);
        dswap(&iptr, &Wat(W, k, 0, i), &inc1, &Wat(W, k, 0, best_j), &inc1);
    }

    std::swap(s[i], s[best_j]);
    std::swap(s_ref[i], s_ref[best_j]);
    std::swap(wnorm2[i], wnorm2[best_j]);
    std::swap(p[i], p[best_j]);
}

void update_trailing_state(
    std::vector<double> &A,
    std::vector<double> &W,
    std::vector<double> &wnorm2,
    std::vector<double> &s,
    std::vector<double> &s_ref,
    std::vector<double> &beta_vec,
    std::vector<double> &dots,
    mwSize m,
    mwSize n,
    mwSize k,
    mwSize i,
    double beta_i,
    double norm_recomp_tol
) {
    const mwSize nrem = n - i - 1;
    if (nrem == 0) {
        return;
    }

    const double invdiag = (beta_i != 0.0) ? (1.0 / beta_i) : 0.0;
    for (mwSize t = 0; t < nrem; ++t) {
        const mwSize j = i + 1 + t;
        const double aij = AatConst(A, m, i, j);
        beta_vec[t] = aij * invdiag;
        Wat(W, k, i, j) = beta_vec[t];
    }

    ptrdiff_t inc1 = 1;
    if (i > 0) {
        const ptrdiff_t mm = static_cast<ptrdiff_t>(i);
        const ptrdiff_t nn = static_cast<ptrdiff_t>(nrem);
        const ptrdiff_t lda = static_cast<ptrdiff_t>(k);
        const double one = 1.0;
        const double zero = 0.0;
        const double neg1 = -1.0;
        char trans = 'T';
        double *Wprefix = &Wat(W, k, 0, i + 1);
        double *wpivot = &Wat(W, k, 0, i);
        dgemv(&trans, &mm, &nn, &one, Wprefix, &lda, wpivot, &inc1, &zero, dots.data(), &inc1);
        dger(&mm, &nn, &neg1, wpivot, &inc1, beta_vec.data(), &inc1, Wprefix, &lda);
    }

    const double wcoeff = wnorm2[i] + 1.0;
    for (mwSize t = 0; t < nrem; ++t) {
        const mwSize j = i + 1 + t;
        const double b = beta_vec[t];

        double wn = 0.0;
        if (i > 0) {
            wn = wnorm2[j] - 2.0 * b * dots[t] + b * b * wcoeff;
        } else {
            wn = wnorm2[j] + b * b;
        }
        wnorm2[j] = std::max(wn, 0.0);

        const double old_s = s[j];
        if (!(old_s > 0.0)) {
            s[j] = 0.0;
            continue;
        }

        const double alpha = AatConst(A, m, i, j);
        const double sj = std::max(old_s - alpha * alpha, 0.0);
        s[j] = sj;

        // LAPACK-style guard: refresh exact norm after sufficient decay.
        if (sj <= s_ref[j] * norm_recomp_tol) {
            double exact = 0.0;
            if (i + 1 < m) {
                const ptrdiff_t len = static_cast<ptrdiff_t>(m - i - 1);
                // For short tails, inline accumulation avoids BLAS call overhead.
                constexpr ptrdiff_t kRecomputeLoopThreshold = 256;
                if (len <= kRecomputeLoopThreshold) {
                    double acc = 0.0;
                    for (mwSize r = i + 1; r < m; ++r) {
                        const double v = AatConst(A, m, r, j);
                        acc = std::fma(v, v, acc);
                    }
                    exact = acc;
                } else {
                    const double nrm = dnrm2(&len, &Aat(A, m, i + 1, j), &inc1);
                    exact = nrm * nrm;
                }
            }
            s_ref[j] = exact;
            s[j] = exact;
        }
    }
}

}  // namespace

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs < 1) {
        fail("bsqr:NotEnoughInputs", "bsqr requires at least one input matrix A.");
    }
    if (nlhs > 4) {
        fail("bsqr:TooManyOutputs", "bsqr supports at most 4 outputs.");
    }

    const mxArray *Ain = prhs[0];
    if (!(mxIsDouble(Ain) && !mxIsComplex(Ain) && mxGetNumberOfDimensions(Ain) == 2)) {
        fail("bsqr:InvalidInput", "A must be a real double matrix.");
    }

    const mwSize m = mxGetM(Ain);
    const mwSize n = mxGetN(Ain);

    Options opt;
    parse_options(Ain, nrhs, prhs, opt);

    Workspace &ws = workspace();
    const double *aptr = mxGetPr(Ain);
    const size_t mn = static_cast<size_t>(m) * static_cast<size_t>(n);
    ws.A.resize(mn);
    std::copy(aptr, aptr + mn, ws.A.data());
    std::vector<double> &A = ws.A;

    if (opt.check_finite) {
        for (double v : A) {
            if (!std::isfinite(v)) {
                fail("bsqr:NonFiniteInput", "A contains non-finite values.");
            }
        }
    }

    const mwSize k = opt.k;
    ws.tau.assign(k, 0.0);
    ws.W.assign(static_cast<size_t>(k) * static_cast<size_t>(n), 0.0);
    ws.wnorm2.assign(n, 0.0);
    ws.s.assign(n, 0.0);
    ws.s_ref.assign(n, 0.0);
    ws.p.resize(n);

    std::vector<double> &tau = ws.tau;
    std::vector<double> &W = ws.W;
    std::vector<double> &wnorm2 = ws.wnorm2;
    std::vector<double> &s = ws.s;
    std::vector<double> &s_ref = ws.s_ref;
    std::vector<mwSize> &p = ws.p;

    ptrdiff_t inc1 = 1;
    for (mwSize j = 0; j < n; ++j) {
        p[j] = j;
        const ptrdiff_t len = static_cast<ptrdiff_t>(m);
        const double nj = dnrm2(&len, &Aat(A, m, 0, j), &inc1);
        const double sj = nj * nj;
        s[j] = sj;
        s_ref[j] = sj;
    }

    ws.beta_vec.resize(static_cast<size_t>(n));
    ws.dots.resize(static_cast<size_t>(n));
    ws.hh_work.resize(static_cast<size_t>(n));
    std::vector<double> &beta_vec = ws.beta_vec;
    std::vector<double> &dots = ws.dots;
    std::vector<double> &hh_work = ws.hh_work;

    for (mwSize i = 0; i < k; ++i) {
        const mwSize best_j = select_pivot_column(s, wnorm2, i, n);
        swap_pivot_state(A, W, s, s_ref, wnorm2, p, m, k, i, best_j);

        double tau_i = 0.0, beta_i = 0.0;
        householder_column(A, m, i, tau_i, beta_i);
        tau[i] = tau_i;

        if (tau_i != 0.0 && i + 1 < n) {
            Aat(A, m, i, i) = 1.0;
            apply_householder_left(A, m, n, i, tau_i, hh_work);
        }

        Aat(A, m, i, i) = beta_i;
        s[i] = beta_i * beta_i;
        s_ref[i] = s[i];

        update_trailing_state(A, W, wnorm2, s, s_ref, beta_vec, dots, m, n, k, i, beta_i, opt.norm_recomp_tol);
    }

    if (nlhs == 0) {
        return;
    }

    mxArray *R = make_R(A, m, n, k);
    if (nlhs == 1) {
        plhs[0] = R;
        return;
    }

    mxArray *Q = make_Q(A, m, k, tau, ws.q_work);
    plhs[0] = Q;
    plhs[1] = R;

    if (nlhs >= 3) {
        plhs[2] = opt.pivot_vector ? make_pivot_vector(p) : make_pivot_matrix(p);
    }

    if (nlhs >= 4) {
        plhs[3] = make_rinv(W, k, n, opt.return_rinv_r12);
    }
}
