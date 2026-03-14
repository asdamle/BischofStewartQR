function tf = bsqr_bench_parse_bool(v, error_id)
%BSQR_BENCH_PARSE_BOOL Parse logical-like values.

if nargin < 2
    error_id = 'bsqr_benchmark:InvalidBool';
end

if islogical(v)
    tf = logical(v);
    return;
end
if isnumeric(v)
    tf = logical(v ~= 0);
    return;
end

s = lower(strtrim(string(v)));
if s == "1" || s == "true" || s == "yes" || s == "on"
    tf = true;
elseif s == "0" || s == "false" || s == "no" || s == "off" || s == ""
    tf = false;
else
    error(error_id, 'Could not parse boolean value: %s', string(v));
end
end
