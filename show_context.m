function show_context(ascFile, lineNums, win)
if nargin < 3, win = 20; end
txt = fileread(ascFile);
lines = regexp(txt, '\r\n|\n', 'split');

for k = 1:numel(lineNums)
    i = lineNums(k);
    a = max(1, i-win);
    b = min(numel(lines), i+win);
    fprintf('\n===== CONTEXT around line %d (showing %d..%d) =====\n', i, a, b);
    for j = a:b
        fprintf('[%d] %s\n', j, lines{j});
    end
end
end
