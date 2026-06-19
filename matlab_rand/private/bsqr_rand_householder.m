function [tau, beta, vtail] = bsqr_rand_householder(col)
%BSQR_RAND_HOUSEHOLDER Elementary reflector zeroing col(2:end).
%   Same convention as matlab/private/bsqr_mfile.m: H = I - tau v v' with
%   v = [1; vtail], H*col = beta*e1. Returns the tail v(2:end) explicitly so
%   the caller can pack it into the reflector store.

nlen = numel(col);
alpha = col(1);
vtail = zeros(max(nlen - 1, 0), 1);
if nlen <= 1
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
vtail = col(2:end) / (alpha - beta);
end
