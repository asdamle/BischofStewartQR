function val = bsqr_bench_getenv_default(key, default)
%BSQR_BENCH_GETENV_DEFAULT Return environment value or fallback default.

val = getenv(key);
if isempty(val)
    val = default;
end
end
