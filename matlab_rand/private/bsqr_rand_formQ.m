function Q = bsqr_rand_formQ(reflectors, ncols)
%BSQR_RAND_FORMQ Materialize Q from the m-file kernel's reflector store.
%   Internal helper: bsqr_rand forms its lazy Q output from the reflector
%   struct returned by bsqr_rand_mfile (the MEX backend uses LAPACK dorgqr
%   instead). Q = BSQR_RAND_FORMQ(REFLECTORS) returns the economy m-by-k
%   factor; BSQR_RAND_FORMQ(REFLECTORS, NCOLS) its leading NCOLS columns.
%   Q satisfies Q' * A(:,P(1:k)) = R11 for the bsqr_rand outputs.

V = reflectors.V;
tau = reflectors.tau;
[m, k] = size(V);
if nargin < 2 || isempty(ncols)
    ncols = k;
end

Q = eye(m, ncols);
for t = k:-1:1
    tt = tau(t);
    if tt ~= 0
        v = V(:, t);          % zeros above row t, 1 at row t, tail below
        Q = Q - (tt * v) * (v' * Q);
    end
end
end
