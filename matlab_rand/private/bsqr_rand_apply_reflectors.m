function Xred = bsqr_rand_apply_reflectors(X, V, tau, nsel, m)
%BSQR_RAND_APPLY_REFLECTORS Apply the first NSEL accumulated reflectors to X.
%   Computes Q_nsel * X where Q_nsel = H_nsel ... H_1 (the transforms already
%   applied to the selected columns), bringing the original candidate columns
%   X (m-by-b) into the current partially-reduced frame. Rows 1:nsel of the
%   result are the R12 part (r), rows nsel+1:m are the unreduced tail.
%
%   V is the m-by-k reflector store in unit-diagonal packed form: column t has
%   zeros in rows 1:t-1, a 1 in row t, and the reflector tail below. This is a
%   sequence of BLAS-2 rank-1 updates on the whole block (one per reflector).

Xred = X;
for t = 1:nsel
    tt = tau(t);
    if tt ~= 0
        v = V(:, t);
        Xred = Xred - (tt * v) * (v' * Xred);
    end
end
end
