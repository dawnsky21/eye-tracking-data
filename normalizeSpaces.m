function s = normalizeSpaces(s)
    s = string(s);
    s = replace(s, char(160), " ");     % NBSP -> space
    s = regexprep(s, "\s+", " ");       % collapse whitespace
    s = strtrim(s);
end

