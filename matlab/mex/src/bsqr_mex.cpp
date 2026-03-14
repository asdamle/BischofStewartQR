#include "mex.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Options {
    mwSize k = 0;
    bool k_set = false;
    bool return_rinv_r12 = false;
    bool pivot_vector = false;  // false => matrix
    double norm_recomp_tol = std::sqrt(std::numeric_limits<double>::epsilon());
    bool check_finite = true;
};

inline double &Aat(std::vector<double> &A, mwSize m, mwSize r, mwSize c) {
    return A[r + c * m];
}

inline double AatConst(const std::vector<double> &A, mwSize m, mwSize r, mwSize c) {
    return A[r + c * m];
}

inline double &Wat(std::vector<double> &W, mwSize k, mwSize r, mwSize c) {
    return W[r + c * k];
}

inline double WatConst(const std::vector<double> &W, mwSize k, mwSize r, mwSize c) {
    return W[r + c * k];
}

void fail(const char *id, const char *msg) {
    mexErrMsgIdAndTxt(id, "%s", msg);
}

bool scalar_to_bool(const mxArray *a, const char *id) {
    if (!(mxIsLogicalScalar(a) || (mxIsNumeric(a) && mxIsScalar(a)))) {
        mexErrMsgIdAndTxt(id, "Expected a logical/numeric scalar.");
    }
    double v = mxIsLogical(a) ? (mxIsLogicalScalarTrue(a) ? 1.0 : 0.0) : mxGetScalar(a);
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
        std::string name = to_lower(get_string(prhs[i], "bsqr:InvalidOptionName"));
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
            opt.k_set = true;
        } else if (name == "return_rinv_r12") {
            opt.return_rinv_r12 = scalar_to_bool(val, "bsqr:InvalidReturnRinvR12");
        } else if (name == "pivot_format") {
            std::string pf = to_lower(get_string(val, "bsqr:InvalidPivotFormat"));
            if (pf == "vector") {
                opt.pivot_vector = true;
            } else if (pf == "matrix") {
                opt.pivot_vector = false;
            } else {
                fail("bsqr:InvalidPivotFormat", "pivot_format must be \"matrix\" or \"vector\".");
            }
        } else if (name == "backend") {
            // Accepted for API compatibility; ignored inside MEX.
            std::string b = to_lower(get_string(val, "bsqr:InvalidBackend"));
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
    const mwSize n = m - i;
    const double alpha = AatConst(A, m, i, i);
    if (n == 1) {
        tau = 0.0;
        beta = alpha;
        return;
    }

    double xnorm2 = 0.0;
    for (mwSize r = i + 1; r < m; ++r) {
        const double v = AatConst(A, m, r, i);
        xnorm2 += v * v;
    }
    const double xnorm = std::sqrt(xnorm2);
    if (xnorm == 0.0) {
        tau = 0.0;
        beta = alpha;
        return;
    }

    const double sgn = (alpha >= 0.0) ? 1.0 : -1.0;
    beta = -sgn * std::hypot(alpha, xnorm);
    tau = (beta - alpha) / beta;
    const double scale = 1.0 / (alpha - beta);
    for (mwSize r = i + 1; r < m; ++r) {
        Aat(A, m, r, i) *= scale;
    }
    Aat(A, m, i, i) = beta;
}

void apply_householder_left(std::vector<double> &A, mwSize m, mwSize n, mwSize i, double tau) {
    if (tau == 0.0 || i + 1 >= n) {
        return;
    }

    for (mwSize j = i + 1; j < n; ++j) {
        double dot = AatConst(A, m, i, j);
        for (mwSize r = i + 1; r < m; ++r) {
            dot += AatConst(A, m, r, i) * AatConst(A, m, r, j);
        }
        const double coeff = tau * dot;
        Aat(A, m, i, j) -= coeff;
        for (mwSize r = i + 1; r < m; ++r) {
            Aat(A, m, r, j) -= coeff * AatConst(A, m, r, i);
        }
    }
}

void build_q(const std::vector<double> &A, mwSize m, mwSize k, const std::vector<double> &tau, std::vector<double> &Q) {
    Q.assign(m * k, 0.0);
    for (mwSize c = 0; c < k; ++c) {
        Q[c + c * m] = 1.0;
    }

    for (mwSignedIndex ii = static_cast<mwSignedIndex>(k) - 1; ii >= 0; --ii) {
        const mwSize i = static_cast<mwSize>(ii);
        const double tau_i = tau[i];
        if (tau_i == 0.0) {
            continue;
        }
        for (mwSize col = 0; col < k; ++col) {
            double dot = Q[i + col * m];
            for (mwSize r = i + 1; r < m; ++r) {
                dot += AatConst(A, m, r, i) * Q[r + col * m];
            }
            const double coeff = tau_i * dot;
            Q[i + col * m] -= coeff;
            for (mwSize r = i + 1; r < m; ++r) {
                Q[r + col * m] -= coeff * AatConst(A, m, r, i);
            }
        }
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

mxArray *make_Q(const std::vector<double> &A, mwSize m, mwSize k, const std::vector<double> &tau) {
    mxArray *Qm = mxCreateDoubleMatrix(m, k, mxREAL);
    double *qptr = mxGetPr(Qm);
    std::vector<double> Q;
    build_q(A, m, k, tau, Q);
    std::copy(Q.begin(), Q.end(), qptr);
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
    std::fill(eptr, eptr + n * n, 0.0);
    for (mwSize j = 0; j < n; ++j) {
        const mwSize r = p[j];
        eptr[r + j * n] = 1.0;
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

    const double *aptr = mxGetPr(Ain);
    std::vector<double> A(aptr, aptr + static_cast<size_t>(m) * static_cast<size_t>(n));

    if (opt.check_finite) {
        for (double v : A) {
            if (!std::isfinite(v)) {
                fail("bsqr:NonFiniteInput", "A contains non-finite values.");
            }
        }
    }

    const mwSize k = opt.k;
    std::vector<double> tau(k, 0.0);
    std::vector<double> W(k * n, 0.0);
    std::vector<double> wnorm2(n, 0.0), s(n, 0.0), s_ref(n, 0.0);
    std::vector<mwSize> p(n, 0);

    for (mwSize j = 0; j < n; ++j) {
        p[j] = j;
        double sj = 0.0;
        for (mwSize r = 0; r < m; ++r) {
            const double v = AatConst(A, m, r, j);
            sj += v * v;
        }
        s[j] = sj;
        s_ref[j] = sj;
    }

    std::vector<double> beta_vec;
    std::vector<double> dots;

    for (mwSize i = 0; i < k; ++i) {
        mwSize best_j = i;
        double best_c = std::numeric_limits<double>::infinity();
        for (mwSize j = i; j < n; ++j) {
            const double sj = s[j];
            const double cj = (sj > 0.0) ? (1.0 + wnorm2[j]) / sj : std::numeric_limits<double>::infinity();
            if (cj < best_c) {
                best_c = cj;
                best_j = j;
            }
        }

        if (best_j != i) {
            for (mwSize r = 0; r < m; ++r) {
                std::swap(Aat(A, m, r, i), Aat(A, m, r, best_j));
            }
            for (mwSize r = 0; r < i; ++r) {
                std::swap(Wat(W, k, r, i), Wat(W, k, r, best_j));
            }
            std::swap(s[i], s[best_j]);
            std::swap(s_ref[i], s_ref[best_j]);
            std::swap(wnorm2[i], wnorm2[best_j]);
            std::swap(p[i], p[best_j]);
        }

        double tau_i = 0.0, beta_i = 0.0;
        householder_column(A, m, i, tau_i, beta_i);
        tau[i] = tau_i;

        if (tau_i != 0.0 && i + 1 < n) {
            Aat(A, m, i, i) = 1.0;
            apply_householder_left(A, m, n, i, tau_i);
        }

        Aat(A, m, i, i) = beta_i;
        s[i] = beta_i * beta_i;
        s_ref[i] = s[i];

        const mwSize nrem = n - i - 1;
        if (nrem == 0) {
            continue;
        }

        beta_vec.assign(nrem, 0.0);
        const double invdiag = (beta_i != 0.0) ? (1.0 / beta_i) : 0.0;
        for (mwSize t = 0; t < nrem; ++t) {
            const mwSize j = i + 1 + t;
            const double aij = AatConst(A, m, i, j);
            beta_vec[t] = aij * invdiag;
            Wat(W, k, i, j) = beta_vec[t];
        }

        dots.assign(nrem, 0.0);
        if (i > 0) {
            for (mwSize t = 0; t < nrem; ++t) {
                const mwSize j = i + 1 + t;
                double d = 0.0;
                for (mwSize r = 0; r < i; ++r) {
                    d += WatConst(W, k, r, j) * WatConst(W, k, r, i);
                }
                dots[t] = d;
            }
            for (mwSize t = 0; t < nrem; ++t) {
                const mwSize j = i + 1 + t;
                for (mwSize r = 0; r < i; ++r) {
                    Wat(W, k, r, j) -= WatConst(W, k, r, i) * beta_vec[t];
                }
            }
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
            if (wn < 0.0) {
                wn = 0.0;
            }
            wnorm2[j] = wn;

            const double old_s = s[j];
            if (!(old_s > 0.0)) {
                s[j] = 0.0;
                continue;
            }

            const double alpha = AatConst(A, m, i, j);
            double sj = old_s - alpha * alpha;
            if (sj < 0.0) {
                sj = 0.0;
            }
            s[j] = sj;

            if (sj <= s_ref[j] * opt.norm_recomp_tol) {
                double exact = 0.0;
                if (i + 1 < m) {
                    for (mwSize r = i + 1; r < m; ++r) {
                        const double v = AatConst(A, m, r, j);
                        exact += v * v;
                    }
                }
                s_ref[j] = exact;
                s[j] = exact;
            }
        }
    }

    const mxArray *Q = nullptr;
    const mxArray *R = make_R(A, m, n, k);
    if (nlhs >= 2) {
        Q = make_Q(A, m, k, tau);
    }

    if (nlhs == 1) {
        plhs[0] = const_cast<mxArray *>(R);
        return;
    }

    plhs[0] = const_cast<mxArray *>(Q);
    plhs[1] = const_cast<mxArray *>(R);

    if (nlhs >= 3) {
        plhs[2] = opt.pivot_vector ? make_pivot_vector(p) : make_pivot_matrix(p);
    }

    if (nlhs >= 4) {
        plhs[3] = make_rinv(W, k, n, opt.return_rinv_r12);
    }
}
