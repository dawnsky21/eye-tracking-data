function diag_nonsample_numeric_lines(ascFile, maxKeep)
% Robust version: handles variable token counts per line.

if nargin < 2, maxKeep = 5000; end
if nargin < 1 || isempty(ascFile)
    error('ascFile을 넘겨줘야 합니다. 예: diag_nonsample_numeric_lines(ascFile, 8000)');
end

txt = fileread(ascFile);
lines = regexp(txt, '\r\n|\n|\r', 'split');
n = numel(lines);

s = string(lines(:));
sTrim = strip(s);

% 1) "숫자로 시작" 후보
isNumStart = ~cellfun('isempty', regexp(cellstr(sTrim), '^[+-]?\d', 'once'));
candIdx = find(isNumStart);
cand = sTrim(candIdx);

% 2) 라인별 토큰화(길이 달라도 안전)
Kc = numel(cand);
firstTok = strings(Kc,1);
nTok = zeros(Kc,1);

for i = 1:Kc
    tks = regexp(cand(i), '\S+', 'match');  % 공백/탭 포함 어떤 whitespace든 OK
    nTok(i) = numel(tks);
    if ~isempty(tks)
        firstTok(i) = string(tks{1});
    else
        firstTok(i) = "";
    end
end

% 3) "정상 sample 후보" (time 정수 + 최소 4토큰)
isIntTime = ~cellfun('isempty', regexp(cellstr(firstTok), '^\d+$', 'once'));
isLikelySample = isIntTime & (nTok >= 4);

% 4) Non-sample numeric-start
nonMask = ~isLikelySample;
nonIdx = candIdx(nonMask);
nonLines = sTrim(nonIdx);

fprintf('Total lines: %d\n', n);
fprintf('Numeric-start candidates: %d\n', numel(candIdx));
fprintf('Likely-sample among numeric-start: %d\n', sum(isLikelySample));
fprintf('Non-sample numeric-start: %d\n', numel(nonIdx));

% --- 샘플링 ---
K = min(maxKeep, numel(nonIdx));
keepIdx = nonIdx(1:K);
keepLines = sTrim(keepIdx);

% 5) 패턴 분류
type = strings(K,1);
for i = 1:K
    L = keepLines(i);

    if contains(L, "e-") || contains(L, "e+") || ~isempty(regexp(L, '^\d+(\.\d+)?e[+-]?\d+', 'once'))
        type(i) = "scientific_notation";
    elseif contains(L, ",")
        type(i) = "has_comma";
    elseif contains(L, "→") || contains(L, char(8594))
        type(i) = "has_arrow";
    elseif contains(L, char(9)) % tab
        type(i) = "has_tab";
    elseif ~isempty(regexp(L, '^\d+(\.\d+)?$', 'once'))
        type(i) = "single_number_only";
    else
        type(i) = "other";
    end
end

% 6) 요약표
ut = unique(type);
cnt = zeros(numel(ut),1);
for j = 1:numel(ut)
    cnt(j) = sum(type == ut(j));
end
T = table(ut, cnt, 'VariableNames', {'type','count'});
T = sortrows(T, 'count','descend');
disp(T);

% 7) 타입별 대표 예시(라인 번호 포함)
topTypes = T.type(1:min(5,height(T)));
for j = 1:numel(topTypes)
    tt = topTypes(j);
    fprintf('\n--- TYPE: %s (show up to 5) ---\n', tt);
    ex = find(type == tt, 5, 'first');
    for k2 = 1:numel(ex)
        ii = ex(k2);
        fprintf('[line %d] %s\n', keepIdx(ii), keepLines(ii));
    end
end

% 8) 어디에 몰려있는지(파일 라인 인덱스 분포)
figure;
histogram(nonIdx, 50);
xlabel('Line index in ASC'); ylabel('Count');
title('Where non-sample numeric lines occur (by file line index)');

end

