function [p, reflectors, R11, stats, R12] = bsqr_rand_mfile(A, opts)
%BSQR_RAND_MFILE Randomized Bischof-Stewart column selection (reference m-file).
%
%   Selects K columns of A without maintaining R11^{-1}R12 / column norms for
%   every column at every step. Instead it tracks only the running squared
%   inverse Frobenius norm f2 = ||R11^{-1}||_F^2 and samples candidate columns
%   in blocks, brought into the current reduced frame by applying the
%   accumulated reflectors, keeping those whose increment holds f2 under the
%   per-step threshold (see bsqr_rand_threshold and docs/RANDOMIZED_BSQR_PLAN.md).
%
%   Two paths (opts.batched, default true):
%     batched=true  - in-block BSQR: each sampled block is reduced once, then
%                     BSQR runs within it, taking as many columns as the bound
%                     allows (f2/threshold updated after EACH in-block pick)
%                     before resampling. Amortizes the per-block apply.
%     batched=false - single-select: one accepted column per sampled block.
%   If a whole pass over the remaining columns finds none acceptable (only near
%   rounding ties) it falls back to the global minimizer, guaranteed to exist.
%
%   Returns the column subset (p(1:k)), the accumulated reflectors, R11, an
%   instrumentation struct, and -- only when requested -- R12.

Awork = double(A);
[m, n] = size(Awork);
k = opts.k;

if ~isempty(opts.seed)
    rng(opts.seed);
end

V = zeros(m, k);          % reflector store, unit-diagonal packed form
R11 = zeros(k, k);        % upper-triangular factor of the selected columns
tau = zeros(k, 1);
selected = zeros(1, k);
remaining = 1:n;          % original indices still available
f2 = 0.0;                 % running ||R11^{-1}||_F^2

stats = local_init_stats(k, n);

needweights = strcmp(opts.sampling, 'normweighted');
g = [];
if needweights
    g = sum(Awork.^2, 1);          % original squared column norms (O(m*n))
end

b = opts.block_size;

if opts.batched
    % ===== Batched in-block BSQR (default kernel path) =====
    % Each sampled block is brought to the current frame ONCE (one reflector
    % apply); we then run BSQR within it, selecting the in-block minimizer
    % greedily -- incrementing f2 and re-deriving the threshold after EACH
    % selection -- while the per-step bound is met. One expensive apply yields
    % multiple selections (amortizing the dominant cost, like rejection_rpqr).
    nsel = 0;
    since_last_select = 0;          % columns sampled since the last selection
    while nsel < k
        rem_count = n - nsel;
        force_all = since_last_select >= rem_count;   % progress safety net
        if force_all
            ids = remaining; bcount = rem_count;
        else
            bcount = min(b, rem_count);
            ids = sample_block(remaining, bcount, needweights, g);
        end
        Xred = bsqr_rand_apply_reflectors(Awork(:, ids), V, tau, nsel, m);
        since_last_select = since_last_select + bcount;
        cand = true(1, bcount);
        nsel_at_block_start = nsel;
        first = true;
        while nsel < k && any(cand)
            % rho2 / wn2 are recomputed exactly from Xred each in-block step (no
            % running downdate), so opts.norm_recomp_tol -- the MEX's safeguard knob
            % -- has nothing to guard here and is intentionally unused.
            ci = find(cand);
            Xc = Xred(:, ci);
            rho2 = sum(Xc(nsel + 1:m, :).^2, 1);
            if nsel > 0
                Wb = R11(1:nsel, 1:nsel) \ Xc(1:nsel, :);
                wn2 = sum(Wb.^2, 1);
            else
                wn2 = zeros(1, numel(ci));
            end
            cvals = (1 + wn2) ./ rho2;
            theta = bsqr_rand_threshold(f2, nsel, k, n, opts.threshold_mode) * opts.slack;
            [bmin, brel] = min(cvals);
            if ~isfinite(bmin)
                if force_all && first
                    error('bsqr_rand:RankDeficient', ...
                        ['At step %d all remaining columns have ~zero residual; the ', ...
                         'input appears rank-deficient for k=%d.'], nsel + 1, k);
                end
                break;
            end
            forced = false;
            if bmin > theta
                if force_all && first
                    forced = true;     % rounding-tie fallback: take the global min
                else
                    break;             % no in-block column meets the bound -> resample
                end
            end
            first = false;

            step = nsel + 1;
            jblk = ci(brel);
            xsel = Xred(:, jblk);
            R11(1:nsel, step) = xsel(1:nsel);
            [tau_i, beta_i, vtail] = bsqr_rand_householder(xsel(nsel + 1:m));
            R11(step, step) = beta_i;
            V(step, step) = 1;
            if m > step; V(step + 1:m, step) = vtail; end
            tau(step) = tau_i;

            f2 = f2 + bmin;            % <-- f_i incremented per selection, even in-block
            selected(step) = ids(jblk);
            remaining(remaining == ids(jblk)) = [];
            cand(jblk) = false;
            since_last_select = 0;

            rest = find(cand);          % apply the new reflector to the rest of the block
            if ~isempty(rest) && tau_i ~= 0
                v = [zeros(step - 1, 1); 1; vtail];
                Xred(:, rest) = Xred(:, rest) - (tau_i * v) * (v' * Xred(:, rest));
            end

            stats.f2(step) = f2;
            stats.crit(step) = bmin;
            stats.threshold(step) = theta;
            stats.Fhat(step) = (nsel + 1) * (n - nsel) / (k - nsel);
            stats.fallback(step) = forced;
            if nsel == nsel_at_block_start
                stats.samples_tested(step) = bcount;   % attribute the apply to the 1st pick
                stats.rounds(step) = 1;
            end
            nsel = nsel + 1;
        end
    end
else
for nsel = 0:k-1
    step = nsel + 1;
    rem_count = n - nsel;

    theta = bsqr_rand_threshold(f2, nsel, k, n, opts.threshold_mode) * opts.slack;
    Fhat_next = (nsel + 1) * (n - nsel) / (k - nsel);
    stats.threshold(step) = theta;
    stats.Fhat(step) = Fhat_next;

    % (weighted) random visiting order over the remaining columns
    if needweights
        wts = g(remaining);
        keys = -log(rand(1, rem_count)) ./ max(wts, realmin);   % Efraimidis-Spirakis
        [~, ord] = sort(keys);
    else
        ord = randperm(rem_count);
    end
    visit = remaining(ord);

    best_c = inf; best_id = visit(1);
    accept_id = 0; accept_c = inf;
    tested = 0; rounds = 0;

    pos = 1;
    while pos <= rem_count
        bcount = min(b, rem_count - pos + 1);
        ids = visit(pos:pos + bcount - 1);
        Xred = bsqr_rand_apply_reflectors(Awork(:, ids), V, tau, nsel, m);

        rho2 = sum(Xred(nsel + 1:m, :).^2, 1);          % 1 x bcount
        if nsel > 0
            Wb = R11(1:nsel, 1:nsel) \ Xred(1:nsel, :); % nsel x bcount, triangular solve
            wn2 = sum(Wb.^2, 1);
        else
            wn2 = zeros(1, bcount);
        end
        cvals = (1 + wn2) ./ rho2;                      % Inf where rho2 == 0
        rounds = rounds + 1;
        tested = tested + bcount;

        [bmin, bidx] = min(cvals);
        if bmin < best_c
            best_c = bmin; best_id = ids(bidx);
        end

        if strcmp(opts.pick, 'first')
            rel = find(cvals <= theta, 1, 'first');
            if ~isempty(rel)
                accept_id = ids(rel); accept_c = cvals(rel);
                break;
            end
        else  % best_in_block
            if bmin <= theta
                accept_id = ids(bidx); accept_c = bmin;
                break;
            end
        end
        pos = pos + bcount;
    end

    fallback = false;
    if accept_id == 0
        accept_id = best_id; accept_c = best_c;
        fallback = true;
    end

    % Rank guard: a non-finite increment means every remaining column had a
    % ~zero residual (rho^2 = 0), i.e. the input is rank-deficient for this k.
    % The bound cannot be maintained past the numerical rank; fail loudly rather
    % than propagate Inf into f2/R11. (Randomized BSQR targets the full-rank GKS
    % setting; the deterministic kernel's rank_stop has no analogue here.)
    if ~isfinite(accept_c)
        error('bsqr_rand:RankDeficient', ...
            ['At step %d all remaining columns have ~zero residual; the input ', ...
             'appears rank-deficient for k=%d. Reduce k to the numerical rank.'], step, k);
    end

    % Reduce the accepted column and append its reflector / R11 column.
    xred = bsqr_rand_apply_reflectors(Awork(:, accept_id), V, tau, nsel, m);
    R11(1:nsel, step) = xred(1:nsel);
    [tau_i, beta_i, vtail] = bsqr_rand_householder(xred(nsel + 1:m));
    R11(step, step) = beta_i;
    V(step, step) = 1;
    if m > step
        V(step + 1:m, step) = vtail;
    end
    tau(step) = tau_i;

    f2 = f2 + accept_c;
    selected(step) = accept_id;
    remaining(remaining == accept_id) = [];

    stats.f2(step) = f2;
    stats.crit(step) = accept_c;
    stats.samples_tested(step) = tested;
    stats.rounds(step) = rounds;
    stats.fallback(step) = fallback;
end
end

p = [selected, remaining];
R11 = triu(R11);
reflectors = struct('V', V, 'tau', tau, 'm', m, 'k', k);

stats.frob_inv = sqrt(max(f2, 0));
stats.osinsky_bound = sqrt(k * (n - k + 1));
stats.total_tested = sum(stats.samples_tested);
stats.blocks_sampled = sum(stats.rounds);

if nargout >= 5
    if opts.return_r12 && k > 0 && k < n
        Xrem = bsqr_rand_apply_reflectors(Awork(:, remaining), V, tau, k, m);
        R12 = Xrem(1:k, :);
    elseif k < n
        R12 = zeros(k, n - k);
    else
        R12 = zeros(k, 0);
    end
end
end

function stats = local_init_stats(k, n)
stats = struct();
stats.f2 = zeros(1, k);
stats.crit = zeros(1, k);
stats.threshold = zeros(1, k);
stats.Fhat = zeros(1, k);
stats.samples_tested = zeros(1, k);
stats.rounds = zeros(1, k);
stats.fallback = false(1, k);
stats.frob_inv = 0;
stats.osinsky_bound = sqrt(k * (n - k + 1));
stats.total_tested = 0;
end

function ids = sample_block(remaining, bcount, needweights, g)
% Draw bcount distinct columns from `remaining` (weighted by starting squared
% column norms when needweights, else uniform).
rc = numel(remaining);
if needweights
    keys = -log(rand(1, rc)) ./ max(g(remaining), realmin);   % Efraimidis-Spirakis
    [~, ord] = sort(keys);
else
    ord = randperm(rc);
end
ids = remaining(ord(1:bcount));
end
