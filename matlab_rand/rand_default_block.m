function b = rand_default_block(k)
%RAND_DEFAULT_BLOCK Default single-select sampling block size (batched uses k).
%   b = ceil(k/2), clamped to [16, 64].
%
%   Block size trades compute for *realized* selection quality: with
%   pick='best_in_block' the accepted pivot is the minimum-criterion column over
%   the block, so a larger block minimizes over more candidates and drives the
%   per-step increment (hence the realized ||R11^{-1}||_F) toward the
%   deterministic greedy -- at higher cost. The guaranteed per-step bound is
%   unchanged (it is enforced by the acceptance threshold at any block size); only
%   the realized conditioning improves. See docs/RANDOMIZED_BSQR_PLAN.md.
%
%   The C++ MEX backend mirrors this formula; keep the two in sync.

b = min(64, max(16, ceil(k / 2)));
end
