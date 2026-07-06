function [Q, R, p, rinv_r12, trace] = bsqr_mfile(A, opts)
%BSQR_MFILE Pure MATLAB BSQR backend.
%   The fifth output is the per-step validation trace (docs/VALIDATION.md
%   V4): trace.crit(i) is the criterion value of the selected pivot and
%   trace.nrecomp(i) the number of norm-recompute events at step i.

Awork = double(A);
[m, n] = size(Awork);
k = opts.k;

p = 1:n;
tau = zeros(k, 1);
W = zeros(k, n);
wnorm2 = zeros(1, n);
[s, s_ref] = init_column_norm_state(Awork, n);
trace = struct('crit', zeros(1, k), 'nrecomp', zeros(1, k));

for i = 1:k
    [best_j, best_c] = select_pivot_column(i, n, s, wnorm2);
    trace.crit(i) = best_c;
    if best_j ~= i
        [Awork, W, s, s_ref, wnorm2, p] = swap_pivot_state(Awork, W, s, s_ref, wnorm2, p, i, best_j);
    end

    [tau_i, beta_i, col] = householder_column(Awork(i:m, i));
    tau(i) = tau_i;
    Awork(i:m, i) = col;

    if tau_i ~= 0 && i < n
        % Apply H_i from the left with v(1)=1 staging in-place.
        Awork(i, i) = 1;
        Awork(i:m, i+1:n) = apply_householder_left(Awork(i:m, i+1:n), Awork(i:m, i), tau_i);
    end

    Awork(i, i) = beta_i;
    s(i) = beta_i * beta_i;
    s_ref(i) = s(i);

    nrem = n - i;
    if nrem == 0
        continue;
    end

    [W, wnorm2, s, s_ref, nrecomp] = update_trailing_state( ...
        Awork, W, wnorm2, s, s_ref, i, m, n, beta_i, opts.norm_recomp_tol);
    trace.nrecomp(i) = nrecomp;
end

if k == 0
    Q = zeros(m, 0);
    R = zeros(0, n);
else
    Q = build_q_from_factors(Awork, tau, k, m);
    R = triu(Awork(1:k, :));
end

if opts.return_rinv_r12
    if k == 0
        rinv_r12 = zeros(0, n);
    elseif k < n
        rinv_r12 = W(1:k, k+1:n);
    else
        rinv_r12 = zeros(k, 0);
    end
else
    rinv_r12 = [];
end

end

function [s, s_ref] = init_column_norm_state(Awork, n)
s = zeros(1, n);
s_ref = zeros(1, n);
for j = 1:n
    nj = norm(Awork(:, j));
    sj = nj * nj;
    s(j) = sj;
    s_ref(j) = sj;
end
end

function [best_j, best_c] = select_pivot_column(i, n, s, wnorm2)
best_j = i;
best_c = inf;
for j = i:n
    sj = s(j);
    if sj > 0
        % Bischof-Stewart criterion: minimize (1 + ||w_j||^2) / ||a_j^(i)||^2.
        cj = (1 + wnorm2(j)) / sj;
    else
        cj = inf;
    end
    if cj < best_c
        best_c = cj;
        best_j = j;
    end
end
end

function [Awork, W, s, s_ref, wnorm2, p] = swap_pivot_state(Awork, W, s, s_ref, wnorm2, p, i, best_j)
Awork(:, [i, best_j]) = Awork(:, [best_j, i]);
if i > 1
    W(1:i-1, [i, best_j]) = W(1:i-1, [best_j, i]);
end
s([i, best_j]) = s([best_j, i]);
s_ref([i, best_j]) = s_ref([best_j, i]);
wnorm2([i, best_j]) = wnorm2([best_j, i]);
p([i, best_j]) = p([best_j, i]);
end

function [W, wnorm2, s, s_ref, nrecomp] = update_trailing_state(Awork, W, wnorm2, s, s_ref, i, m, n, beta_i, norm_recomp_tol)
nrem = n - i;
nrecomp = 0;
alpha = Awork(i, i+1:n);
if beta_i ~= 0
    invdiag = 1 / beta_i;
else
    invdiag = 0;
end
beta_vec = alpha .* invdiag;
W(i, i+1:n) = beta_vec;

if i > 1
    Wprefix = W(1:i-1, i+1:n);
    wpivot = W(1:i-1, i);
    dots = Wprefix.' * wpivot;
    W(1:i-1, i+1:n) = Wprefix - wpivot * beta_vec;
else
    dots = zeros(nrem, 1);
end

wcoeff = wnorm2(i) + 1;
for t = 1:nrem
    j = i + t;
    b = beta_vec(t);

    if i > 1
        d = dots(t);
        wn = wnorm2(j) - 2 * b * d + b * b * wcoeff;
    else
        wn = wnorm2(j) + b * b;
    end
    wnorm2(j) = max(wn, 0);

    old_s = s(j);
    if ~(old_s > 0)
        s(j) = 0;
        continue;
    end

    sj = max(old_s - alpha(t) * alpha(t), 0);
    s(j) = sj;

    % Refresh exact tail norm only after sufficient decay in the running estimate.
    if sj <= s_ref(j) * norm_recomp_tol
        if i < m
            tail = Awork(i+1:m, j);
            sj_exact = dot(tail, tail);
        else
            sj_exact = 0;
        end
        s_ref(j) = sj_exact;
        s(j) = sj_exact;
        nrecomp = nrecomp + 1;
    end
end
end

function [tau, beta, col] = householder_column(col)
% Reflector zeroing col below its first entry: H = I - tau*v*v' with
% v = [1; col(2:end)] (unit-diagonal packed form). On return col(1) = beta
% (the new diagonal) and col(2:end) is the scaled reflector tail; tau = 0
% encodes the identity reflector (length-1 or zero-tail column, col unchanged).
n = numel(col);
alpha = col(1);
if n == 1
    tau = 0;
    beta = alpha;
    return;
end

xnorm = norm(col(2:end));
if xnorm == 0
    tau = 0;
    beta = alpha;
    return;
end

if alpha >= 0
    sgn = 1;
else
    sgn = -1;
end
beta = -sgn * hypot(alpha, xnorm);
tau = (beta - alpha) / beta;
scale = 1 / (alpha - beta);
col(2:end) = col(2:end) * scale;
col(1) = beta;
end

function B = apply_householder_left(B, v, tau)
% Applies H = I - tau*v*v' from the left. v must be the unit-diagonal packed
% form (v(1) = 1, tail below); tau = 0 encodes the identity reflector.
if tau == 0 || isempty(B)
    return;
end
B = B - (tau * v) * (v.' * B);
end

function Q = build_q_from_factors(Afact, tau, ksteps, m)
% Economy Q = H_1*...*H_ksteps applied to eye(m, ksteps), accumulated
% backward (i = ksteps:-1:1) so each reflector is applied once; tau(i) = 0
% reflectors (identity) are skipped.
Q = eye(m, ksteps);
for i = ksteps:-1:1
    tau_i = tau(i);
    if tau_i == 0
        continue;
    end
    v = [1; Afact(i+1:m, i)];
    block = Q(i:m, :);
    Q(i:m, :) = block - (tau_i * v) * (v.' * block);
end
end
