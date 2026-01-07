%% Eye-tracking 문장 읽기 실험 분석 파이프라인 개요
% 이 파일은 EDF → ASC → MATLAB 분석까지 전체 흐름을 정리한 주석 블록이다.
% 실제 코드는 이 구조를 기준으로 함수/스크립트로 나누어 구현한다.

% 0. 참가자·trial 퀄리티 체크 (요약 지표 계산 후 적용)
%   - 참가자 레벨:
%       - comprehension 질문 정답률이 너무 낮은 참가자 제외 (예: < 75% 등 기준 설정)
%       - 전체 trial 중 tracking 실패(눈 미검출 비율, 비정상적으로 짧은 읽기 등)가
%         일정 비율 이상이면 해당 참가자 제외 후보로 표시
%   - trial 레벨:
%       - 문장 길이에 비해 fixation 수가 비정상적으로 적거나 많은 trial 제외
%       - 문장 시작점/끝점 근처에 시선이 전혀 없는 trial은 “제대로 읽지 않은” trial로 제외
%   - 실제 구현 순서:
%       - 6. word × trial 지표 계산 → summary table 생성 후
%         여기에서 참가자/trial 퀄리티 기준을 적용하는 것이 현실적이다.

% 1. EDF → ASC 변환

% 폴더 및 ASC 파일 경로 설정
baseDir = pwd;
subjDir = '20251202_121351_sub01';
ascFileName = '01.asc';
ascFile = fullfile(baseDir, subjDir, ascFileName);

disp(ascFile)   % 실제 경로 확인용
if ~isfile(ascFile)
    error('ASC 파일을 찾을 수 없습니다: %s', ascFile);
end

% --- mainscript 결과(ROI + design 정보) 불러오기 ---
elFile = fullfile(baseDir, subjDir, 'EL_demo.mat');
S = load(elFile);   % 안에 있는 변수 전부 읽기

if ~isfield(S, 'results') || ~isfield(S, 'dp')
    error('EL_demo.mat 안에 results 또는 dp 변수가 없습니다.');
end
results = S.results;
dp      = S.dp;

% --- Main 시트 읽기 (designRow 매핑 + target용) ---
xlsxPath = fullfile(baseDir, 'Experimental stimulus_실험용_수정본7.xlsx');
if ~isfile(xlsxPath)
    error('엑셀 파일을 찾을 수 없습니다: %s', xlsxPath);
end
Tmain = readtable(xlsxPath, 'Sheet','Main', 'TextType','string');

% --- designRow 확보: (1) Final.designRow 있으면 그대로, 없으면 문장 매칭으로 복구 ---
if isfield(S, 'Final') && isfield(S.Final, 'designRow')
    Final     = S.Final;
    designRow = double(Final.designRow);
else
    warning('Final/designRow 없음 → Pilot3 호환 모드: 문장 내용으로 designRow를 복구합니다.');

    nTrials = numel(results.words);
    designRow = zeros(nTrials,1);

    % 1) 엑셀 문장 → 단어 리스트 전처리
    nRows = height(Tmain);
    sentWords = cell(nRows,1);
    for r = 1:nRows
        s = string(Tmain.sentence(r));
        s = normalizeSpaces(s);                 
        if exist('cleanSentence','file')
            s = cleanSentence(s);
        end
        s = normalizeSpaces(s);                 
        ws = split(s);                          % strtrim 불필요
        sentWords{r} = string(ws(:));
    end

    % 2) 각 trial의 results.words{t}와 매칭되는 row 찾기
    for t = 1:nTrials
        wTrial = string(results.words{t}(:));
        if isempty(wTrial)
            continue;
        end

        matchIdx = [];
        for r = 1:nRows
            wRow = sentWords{r};
            if numel(wRow) ~= numel(wTrial)
                continue;
            end
            if all(wRow == wTrial)
                matchIdx(end+1) = r; %#ok<AGROW>
            end
        end

        if numel(matchIdx) == 1
            designRow(t) = matchIdx;
        elseif isempty(matchIdx)
            warning('Trial %d: 엑셀 Main에서 매칭 문장을 찾지 못했습니다.', t);
        else
            warning('Trial %d: 엑셀 Main에서 매칭 row가 여러 개(%d개)입니다. 첫 번째만 사용.', ...
                    t, numel(matchIdx));
            designRow(t) = matchIdx(1);
        end
    end
end

% === (A) 확인할 trial 지정 ===
trialsToInspect = [105 215];

% === (B) 엑셀 key 미리 만들어두기 (Tmain.sentence 기준) ===
excelKey = strings(height(Tmain),1);
for r = 1:height(Tmain)
    excelKey(r) = makeSentenceKey(normalizeSpaces(Tmain.sentence(r)));  
end

% === (C) trial별로 원문/키 출력 + (있으면) 가장 가까운 후보 몇 개 보여주기 ===
for t = trialsToInspect
    fprintf('\n============================\n');
    fprintf('[INSPECT] Trial %d\n', t);

    % 1) trial 원문(단어열) 출력
    if t <= numel(results.words) && ~isempty(results.words{t})
        wTrial = string(results.words{t}(:));
        disp("wTrial (words) = ");
        disp(wTrial');
        trialSentence = strjoin(wTrial, " ");
    else
        disp("wTrial is empty or trial index out of range.");
        trialSentence = "";
    end

    trialSentence = normalizeSpaces(trialSentence);  

    % 2) trial key 출력
    tKey = makeSentenceKey(trialSentence);
    fprintf('trialSentence = %s\n', trialSentence);
    fprintf('trialKey      = %s\n', tKey);

    % 3) 엑셀에서 exact key 매칭되는 row 있는지
    hit = find(excelKey == tKey, 1);
    if ~isempty(hit)
        fprintf('[MATCH] Excel row = %d\n', hit);
        fprintf('excelSentence(raw) = %s\n', string(Tmain.sentence(hit)));
        continue;
    end

    fprintf('[NO EXACT MATCH] Showing a few candidates by similarity...\n');

    % 4) 간단 유사도: 공통 문자 비율(대충)로 상위 5개 후보
    %    (정교한 edit distance가 필요하면 추가해줄게)
    scores = zeros(height(Tmain),1);
    for r = 1:height(Tmain)
        s = excelKey(r);
        if strlength(s)==0 || strlength(tKey)==0
            scores(r) = 0;
        else
            % 공통 토큰 기반(공백 분할 key) 점수
            a = split(tKey, " ");
            b = split(s, " ");
            scores(r) = numel(intersect(a,b)) / max(numel(a),1);
        end
    end

    [~, idx] = sort(scores, 'descend');
    topK = idx(1:min(5,numel(idx)));

    fprintf('Top candidates (row | score | sentence):\n');
    for k = 1:numel(topK)
        r = topK(k);
        fprintf('  %4d | %.3f | %s\n', r, scores(r), string(Tmain.sentence(r)));
    end
end

% --- designRow 범위 체크 + targetIdxPerTrial 만들기 ---
bad = (designRow < 1) | (designRow > height(Tmain));
if any(bad)
    warning('일부 designRow가 Main 시트 범위를 벗어났습니다(비매칭 trial일 수 있음).');
end

% 1) 엑셀에서 target index 벡터 가져오기
%    네가 보여준 엑셀 헤더가
%    sentence, freq, valence, is_catch, target_word, target_idx
%    였으니까, 우선 target_idx를 찾도록 짜자.
varNameCandidates = {'target_idx','targetIndex','target_idx_','targetPos'};  % 열 이름이 다를 경우 대비
nameIdx = find(ismember(Tmain.Properties.VariableNames, varNameCandidates), 1);

if isempty(nameIdx)
    error('Main 시트에서 target index 변수(target_idx 등)를 찾지 못했습니다.');
end

targetVarName = Tmain.Properties.VariableNames{nameIdx};
tmp = double(Tmain.(targetVarName));   % 엑셀에서 target index 벡터

% 2) NaN은 catch trial로 보고 0으로
tmp(isnan(tmp)) = 0;   % catch trial → 0

% 3) trial별 targetIdxPerTrial 채우기
targetIdxPerTrial = zeros(size(designRow));
valid = designRow >= 1 & designRow <= height(Tmain);
targetIdxPerTrial(valid) = tmp(designRow(valid));

% results.wordRects{t} : 각 trial의 [nWords×4] (L T R B) ROI
% results.words{t}     : 각 trial의 단어 string 리스트(있다면)

% 2. ASC → MATLAB struct (sample + event)
%   → parseAscToStruct 내부에서 sample, fix/sacc/blink, MSG, START/END까지 다 파싱해서
%     subj.sample / subj.event.*에 정리함.

% ---- 화면 해상도(고정) ----
screenRect = [0 0 1920 1080];

subj = parseAscToStruct(ascFile);

% ---- FIX/OK: displayRect status ----
if ~isfield(subj,"displayRect") || any(isnan(subj.displayRect))
    subj.displayRect = screenRect;
    fprintf("[FIX] subj.displayRect missing → forced to [%d %d %d %d]\n", screenRect);
else
    fprintf("[OK ] subj.displayRect already set = [%d %d %d %d]\n", subj.displayRect);
end

dr = subj.displayRect;
w = dr(3)-dr(1)+1; 
h = dr(4)-dr(2)+1;
assert(w==1920 && h==1080, "Display size mismatch: %dx%d", w, h);

% === MSG 토큰 빈도 + TRIAL 관련 MSG 샘플(진단) ===
msgText = string({subj.event.msg.text});

tok = regexprep(msgText, "\s+.*$", "");  % 첫 토큰만
[ut,~,ic] = unique(tok);
cnt = accumarray(ic,1);
[~,ord] = sort(cnt,'descend');

k = min(30, numel(ord));
if k == 0
    disp(table(string.empty(0,1), double.empty(0,1), ...
        'VariableNames', {'msgToken','count'}));
else
    tokTop = ut(ord(1:k));  tokTop = tokTop(:);   % column 강제
    cntTop = cnt(ord(1:k)); cntTop = cntTop(:);   % column 강제
    disp(table(tokTop, cntTop, 'VariableNames', {'msgToken','count'}));
end

ix = contains(msgText,"TRIAL","IgnoreCase",true);
disp(msgText(find(ix, 30, 'first')));

% 2.1 샘플 validity 플래그
subj = addSampleValidity(subj, screenRect, 50);

% 2.2 MSG 기반 trial 경계 정의 및 trial struct 만들기
subj = addTrialsFromMsg(subj);
subj = addReadWindowFromMsg(subj);

% [ANCHOR] Trial split (practice/main)  =====
ids = string({subj.trial.id})';
ids = replace(ids, char(160), " ");
ids = regexprep(ids, "\s+", " ");
ids = strtrim(ids);

isPractice = startsWith(ids, "PRACTICE_TRIALID", "IgnoreCase", true);
isMain     = startsWith(ids, "TRIALID",          "IgnoreCase", true);

% ===== CANONICAL (do not overwrite later) =====
mainTrialIdx = find(isMain);        % subj.trial 기준 MAIN 인덱스 (예: 11..234)
practiceTrialIdx = find(isPractice);

% (선택) 즉시 확인
assert(~isempty(mainTrialIdx), "No MAIN trials found.");

fprintf("nTrials total=%d | practice=%d | main=%d\n", numel(ids), sum(isPractice), sum(isMain));
fprintf("ids sample: [%s] | [%s] | [%s]\n", ids(1), ids(min(10,end)), ids(min(11,end)));

assert(any(isPractice) && any(isMain), ...
    "Trial split failed: practice=%d main=%d", sum(isPractice), sum(isMain));

% [CHECK-1] first/last MAIN id
disp(ids(mainTrialIdx([1 numel(mainTrialIdx)])));

% ===== (3) MSG 토큰 빈도/샘플 출력: practice vs main 분리 진단 =====
msgText = string({subj.event.msg.text})';

% (B) 각 msg가 어느 trial에 속하는지 태그
msgTrial = nan(numel(subj.event.msg),1);
for t=1:numel(subj.trial)
    mi = subj.trial(t).msgIdx(:);
    mi = mi(mi>=1 & mi<=numel(subj.event.msg));
    msgTrial(mi) = t;
end

% (C) "MSG 123..." 같은 래퍼 제거(있을 때만)
tmp = regexprep(msgText, "^(MSG→\d+\s+|MSG\s+\d+\s+)", "");
tok = regexprep(tmp, "\s+.*$", "");    % 첫 토큰

% (D) 토큰 빈도 Top30 (전체)
[ut,~,ic] = unique(tok);
cnt = accumarray(ic,1);
[~,ord] = sort(cnt,'descend');
k = min(30, numel(ord));
tokTop = ut(ord(1:k));  tokTop = tokTop(:);
cntTop = cnt(ord(1:k)); cntTop = cntTop(:);
disp(table(tokTop, cntTop, 'VariableNames',{'msgToken','count'}));

% (E) 핵심 토큰별 카운트(Practice/Main 분리) + 샘플
need = ["TRIALID ","PRACTICE_TRIALID ", ...
        "TRIAL_RESULT","PRACTICE_TRIAL_RESULT", ...
        "SENTENCE_ONSET_RIGHTABC","PRACTICE_SENTENCE_ONSET", ...
        "PROMPT_ONSET","PRACTICE_PROMPT_ONSET","RESPONSE","PRACTICE_RESPONSE"];

% --- 안전 마스크: msgTrial이 유효한 trial 인덱스인 경우만 ---
validTrial = isfinite(msgTrial) & msgTrial>=1 & msgTrial<=numel(isPractice);

% msg가 practice/main trial에 속하는지(크기: numMsg×1)
msgIsPractice = false(size(msgTrial));
msgIsMain     = false(size(msgTrial));
msgIsPractice(validTrial) = isPractice(msgTrial(validTrial));
msgIsMain(validTrial)     = isMain(msgTrial(validTrial));

for s = need
    ix = contains(tmp, s, "IgnoreCase", true);

    ip = ix & msgIsPractice;   % ✅ 이제 인덱싱 없음 (안전)
    im = ix & msgIsMain;

    fprintf('\n[%s] PRACTICE n=%d | MAIN n=%d\n', s, sum(ip), sum(im));

    if any(ip), disp(tmp(find(ip,3,'first'))); end
    if any(im), disp(tmp(find(im,3,'first'))); end
end

% (F) main trial에서 PROMPT_ONSET 누락 trial 찾기 (핵심 QC)
missingPrompt = false(numel(mainTrialIdx),1);
for i=1:numel(mainTrialIdx)
    t = mainTrialIdx(i);
    mi = subj.trial(t).msgIdx(:);
    mi = mi(mi>=1 & mi<=numel(subj.event.msg));
    if isempty(mi), missingPrompt(i)=true; continue; end
    txt = tmp(mi);
    missingPrompt(i) = ~any(contains(txt,"PROMPT_ONSET","IgnoreCase",true));
end
fprintf('\n[CHECK] main trials missing PROMPT_ONSET: %d/%d\n', sum(missingPrompt), numel(missingPrompt));
if any(missingPrompt)
    disp(ids(mainTrialIdx(missingPrompt)));
end

fprintf("nTrials total=%d | practice=%d | main=%d\n", numel(ids), sum(isPractice), sum(isMain));

%% 3. drift 기준점(xTrue, yTrue) 자동 계산 (main-only: results.wordRects{1..224} 기반)

assert(isfield(results,'wordRects') && ~isempty(results.wordRects), ...
       "results.wordRects missing/empty");

nMain = numel(results.wordRects);   % main-only 개수

xT = nan(nMain,1);
yT = nan(nMain,1);

for k = 1:nMain
    rects = results.wordRects{k};
    if isempty(rects) || size(rects,2) < 4
        continue;
    end

    r1 = rects(1,:);               % 첫 단어 ROI
    if any(~isfinite(r1)), continue; end

    xT(k) = mean(r1([1 3]));
    yT(k) = mean(r1([2 4]));
end

xTrue = mean(xT,'omitnan');
yTrue = mean(yT,'omitnan');

fprintf('Drift 기준점(xTrue, yTrue) = (%.2f, %.2f) | finite=%d/%d\n', ...
    xTrue, yTrue, sum(isfinite(xT)&isfinite(yT)), nMain);

r = results.wordRects{1}(1,:);
fprintf('[CHECK] word1 ROI example = [%.1f %.1f %.1f %.1f]\n', r);

fprintf('[CHECK] drift applied? gxCorr exists=%d | fix.xCorr exists=%d\n', ...
    isfield(subj.sample,'gxCorr'), isfield(subj.event.fix,'xCorr'));

% (선택) outlier 빠르게 점검 (기준: |dx|>60 or |dy|>40)isPractice = contains(ids,"PRACTICE","IgnoreCase",true);

dx = xT - xTrue;  dy = yT - yTrue;
out = find(abs(dx)>60 | abs(dy)>40);
fprintf('[CHECK] outliers(|dx|>60 or |dy|>40): %d\n', numel(out));
if ~isempty(out)
    disp(table(out, xT(out), yT(out), dx(out), dy(out), ...
        'VariableNames',{'k','x_t','y_t','dx','dy'}));
end

% driftCfg에 반영
driftCfg = struct();
driftCfg.xTrue         = xTrue;
driftCfg.yTrue         = yTrue;
driftCfg.refEvent      = 'SENTENCE_ONSET';
driftCfg.minFixDur     = 60;
driftCfg.maxAbsOffsetX = 60;
driftCfg.maxAbsOffsetY = 40;

subj = applyDriftCorrection(subj, driftCfg);
fprintf('[AFTER applyDrift] gxCorr exists=%d | fix.xCorr exists=%d\n', ...
    isfield(subj.sample,'gxCorr'), isfield(subj.event.fix,'xCorr'));

% 이제 subj.sample.gxCorr / gyCorr, subj.event.fix.xCorr / yCorr,
% subj.trial(t).driftOffset / isBadDrift 등이 채워져 있음.
%
% ---- sanity check (시각 확인) ----
% 몇 개 trial을 골라 raw vs corrected time–gx 플롯을 비교해본다.
mainTrials = mainTrialIdx;
pick = [1 5 10];
pick = pick(pick <= numel(mainTrials));
trialsToCheck = mainTrials(pick);

for ii = 1:numel(trialsToCheck)
    tIdx = trialsToCheck(ii);
    if tIdx > numel(subj.trial)
        continue;
    end

    tr = subj.trial(tIdx);
    sampIdx = tr.sampleIdx;

    tRel = subj.sample.time(sampIdx) - tr.startTime;   % trial 시작 기준 상대 시간
    gxRaw = subj.sample.gx(sampIdx);
    gxCorr = subj.sample.gxCorr(sampIdx);
    gyRaw = subj.sample.gy(sampIdx);
    gyCorr = subj.sample.gyCorr(sampIdx);

    figure;
    subplot(2,1,1);
    plot(tRel, gxRaw, '.-');
    hold on;
    yline(driftCfg.xTrue, 'k--');  % 이상적인 x 위치
    xlabel('Time from trial start (ms)');
    ylabel('gx (raw)');
    title(sprintf('Trial %d - Raw gx', tIdx));

    subplot(2,1,2);
    plot(tRel, gxCorr, '.-');
    hold on;
    yline(driftCfg.xTrue, 'k--');
    xlabel('Time from trial start (ms)');
    ylabel('gx (corrected)');
    title(sprintf('Trial %d - Corrected gx', tIdx));
end

% 추가로, trial별 drift 오프셋 분포를 한 번에 보는 플롯 (요약용)
dx = arrayfun(@(tr) tr.driftOffset(1), subj.trial);
dy = arrayfun(@(tr) tr.driftOffset(2), subj.trial);
badDriftTrials = find([subj.trial.isBadDrift]);

figure;
subplot(1,2,1);
stem(dx, 'filled');
hold on;
yline(driftCfg.maxAbsOffsetX, 'r--');
yline(-driftCfg.maxAbsOffsetX, 'r--');
xlabel('Trial');
ylabel('dx (px)');
title('Drift offset X per trial');

subplot(1,2,2);
stem(dy, 'filled');
hold on;
yline(driftCfg.maxAbsOffsetY, 'r--');
yline(-driftCfg.maxAbsOffsetY, 'r--');
xlabel('Trial');
ylabel('dy (px)');
title('Drift offset Y per trial');

% 4. blink / out-of-range / line 기반 전처리
lineYRange = [500 620];
subj = cleanFixations(subj, [0 0 1920 1080], lineYRange, 50, 70);

%% ========================================================================
% [ANCHOR PATCH-CORE v1] OFFICIAL WORD-LEVEL CORE (single source of truth)
% Place this block immediately AFTER:
% Delete legacy blocks that overlap with this:
%   - "% 5. fixation → word ROI 매핑 (main-only 224 → all trials ... padding=10 ...)"
%   - Any old "makeWordTrialTable + onset-anchor landing" sections
%   - Any old target/spillover code assuming wordTbl.trial == results index
%% ========================================================================

%% ===== [P0] REQUIREMENTS =====
assert(exist('subj','var')==1 && isstruct(subj), 'subj missing');
assert(exist('results','var')==1 && isstruct(results), 'results missing');
assert(exist('targetIdxPerTrial','var')==1, 'targetIdxPerTrial missing');
assert(isfield(results,'wordRects') && ~isempty(results.wordRects), 'results.wordRects missing/empty');

%% ===== [P1] MAIN TRIAL INDEX (canonical) =====
ids = string({subj.trial.id})';
ids = replace(ids, char(160), " ");
ids = regexprep(ids, "\s+", " ");
ids = strtrim(ids);

isMain = startsWith(ids,"TRIALID","IgnoreCase",true);
mainTrialIdx = find(isMain);
assert(~isempty(mainTrialIdx), "No MAIN trials found.");

% results.wordRects is main-only (nMain)
nMain = numel(results.wordRects);
assert(numel(mainTrialIdx)==nMain, "mainTrialIdx(%d) != numel(results.wordRects)(%d)", numel(mainTrialIdx), nMain);

% subj(trial) -> results(main) index mapping
subj2res = nan(numel(subj.trial),1);
subj2res(mainTrialIdx) = 1:nMain;

% results(main) -> subj(trial) mapping
res2subj = nan(nMain,1);
res2subj(:) = mainTrialIdx(:);

fprintf("[P1] nTrials=%d | main=%d | nMain(wordRects)=%d\n", numel(subj.trial), numel(mainTrialIdx), nMain);

%% ===== [P2] ROI MAPPING (OFFICIAL): padding=30 =====
paddingPx = 30;
fprintf("[P2] paddingPx=%d (forced)\n", paddingPx);
assert(paddingPx==30);

% all-trials wordRects cell (main only filled)
wordRectsCell = cell(numel(subj.trial),1);
wordRectsCell(mainTrialIdx) = results.wordRects;

% reset + map
if isfield(subj,'event') && isfield(subj.event,'fix') && ~isempty(subj.event.fix)
    for k=1:numel(subj.event.fix), subj.event.fix(k).word = 0; end
end
subj = mapFixationsToWordROIs(subj, wordRectsCell, paddingPx);

% ===== [P2b] DURATION CLEAN (OFFICIAL): must run AFTER ROI mapping =====
shortThresh = 60;
longThresh  = 1200;
subj = cleanFixationDurations(subj, shortThresh, longThresh);

% quick hit rate on MAIN fixes
mainFixIdx = [];
for t = mainTrialIdx(:)'
    fi = subj.trial(t).fixIdx(:);
    fi = fi(fi>=1 & fi<=numel(subj.event.fix));
    mainFixIdx = [mainFixIdx; fi]; %#ok<AGROW>
end

mainFixIdx = unique(mainFixIdx);
wFix = [subj.event.fix.word]';

fprintf("[P2 ROI] pad=%d | hit(main fixIdx)=%.3f | nFix=%d\n", ...
    paddingPx, mean(wFix(mainFixIdx)>0), numel(mainFixIdx));

% durClean 기준 (공식)
mainFixIdx_dc = [];
for t = mainTrialIdx(:)'
    if isfield(subj.trial,'fixIdxDurClean') && ~isempty(subj.trial(t).fixIdxDurClean)
        fi = subj.trial(t).fixIdxDurClean(:);
    else
        fi = subj.trial(t).fixIdx(:);
    end
    fi = fi(fi>=1 & fi<=numel(subj.event.fix));
    mainFixIdx_dc = [mainFixIdx_dc; fi]; %#ok<AGROW>
end
mainFixIdx_dc = unique(mainFixIdx_dc);

fprintf("[P2 ROI] pad=%d | hit(main durClean)=%.3f | nFix(dc)=%d\n", ...
    paddingPx, mean(wFix(mainFixIdx_dc)>0), numel(mainFixIdx_dc));

fprintf("[P2b durClean] merged=%d | removedShort=%d | longFlag=%d\n", ...
    sum([subj.event.fix.mergedTo]'>0), ...
    sum([subj.event.fix.isRemovedShort]'), ...
    sum([subj.event.fix.isLongFix]'));

%% ===== [P3] BUILD resultsSubj (align results to subj.trial indexing) =====
resultsSubj = results;
resultsSubj.wordRects = cell(numel(subj.trial),1);
resultsSubj.words     = cell(numel(subj.trial),1);

resultsSubj.wordRects(mainTrialIdx) = results.wordRects;
if isfield(results,'words') && ~isempty(results.words)
    resultsSubj.words(mainTrialIdx) = results.words;
end

%% ===== [P4] MAKE wordTbl (all dwell-based measures depend on padding=30 labels) =====
wordTbl = makeWordTrialTable(subj, resultsSubj);

% (defense) if makeWordTrialTable returns trial indices 1..nMain, remap to subj indices
if ~isempty(wordTbl) && min(wordTbl.trial)==1 && max(wordTbl.trial) <= nMain
    wordTbl.trial = res2subj(wordTbl.trial);
end

% keep MAIN only
wordTbl = wordTbl(ismember(wordTbl.trial, mainTrialIdx), :);
fprintf("[P4 wordTbl] rows=%d | trials=%d\n", height(wordTbl), numel(unique(wordTbl.trial)));

%% ===== [P5] firstFixOnsetAbs (landing anchor helper; keep as official field) =====
st = nan(height(wordTbl),1);
ok = isfinite(wordTbl.trial) & wordTbl.trial>=1 & wordTbl.trial<=numel(subj.trial);
st(ok) = arrayfun(@(t) subj.trial(t).startTime, wordTbl.trial(ok));
wordTbl.firstFixOnsetAbs = wordTbl.firstFixOnset + st;
fprintf("[P5 onsetAbs] finite=%d/%d\n", sum(isfinite(wordTbl.firstFixOnsetAbs)), height(wordTbl));

%% ===== [P6] LANDING (OFFICIAL): BETTER STRICT (in-box first among wi-labeled) =====
fix = subj.event.fix;

landingRaw  = nan(height(wordTbl),1);
landingNorm = nan(height(wordTbl),1);
filled      = false(height(wordTbl),1);

for i=1:height(wordTbl)
    tr = wordTbl.trial(i);
    wi = wordTbl.wordIdx(i);

    if wordTbl.nFix(i)==0, continue; end
    if tr<1 || tr>numel(subj.trial) || ~isMain(tr), continue; end

    tw = subj2res(tr);
    if ~isfinite(tw), continue; end

    rects = results.wordRects{tw};
    if isempty(rects) || wi<1 || wi>size(rects,1), continue; end

    xL = rects(wi,1); xR = rects(wi,3);
    width = xR-xL;
    if ~isfinite(width) || width<=0, continue; end

    % dur-clean first
    if isfield(subj.trial,'fixIdxDurClean') && ~isempty(subj.trial(tr).fixIdxDurClean)
        fIdx = subj.trial(tr).fixIdxDurClean(:);
    else
        fIdx = subj.trial(tr).fixIdx(:);
    end
    fIdx = fIdx(fIdx>0 & fIdx<=numel(fix));
    if isempty(fIdx), continue; end

    wSeq = [fix(fIdx).word]';
    posAll = find(wSeq==wi);
    if isempty(posAll), continue; end

    got = false;
    for k = posAll(:)'  % chronological within-trial order
        f0 = fix(fIdx(k));
        xUse = f0.xCorr; if ~isfinite(xUse), xUse = f0.x; end

        if xUse>=xL && xUse<=xR
            landingRaw(i)  = xUse - xL;
            landingNorm(i) = (xUse - xL)/width;
            got = true;
            break;
        end
    end

    if got, filled(i)=true; end
end

wordTbl.landingRaw  = landingRaw;
wordTbl.landingNorm = landingNorm;

fprintf("[P6 landing] eligible(nFix>0)=%d | filled=%d (%.3f)\n", ...
    sum(wordTbl.nFix>0), sum(filled), sum(filled)/max(sum(wordTbl.nFix>0),1));
fprintf("[P6 landing] out-of-[0,1]=%d\n", ...
    sum(isfinite(wordTbl.landingNorm) & (wordTbl.landingNorm<0 | wordTbl.landingNorm>1)));

%% ===== [P7] TARGET / SPILLOVER FLAGS (subj-index aligned) =====
wordTbl.isTarget    = false(height(wordTbl),1);
wordTbl.isSpillover = false(height(wordTbl),1);

targetIdxMain = targetIdxPerTrial(:);          % expected length nMain
targetIdxSubj = zeros(numel(subj.trial),1);

for ii=1:numel(mainTrialIdx)
    tSubj = mainTrialIdx(ii);
    if ii<=numel(targetIdxMain) && isfinite(targetIdxMain(ii))
        targetIdxSubj(tSubj) = targetIdxMain(ii);
    else
        targetIdxSubj(tSubj) = 0;
    end
end

trU = unique(wordTbl.trial);
for t = trU(:)'
    tIdx = targetIdxSubj(t);
    if tIdx<=0, continue; end
    rows = (wordTbl.trial==t);
    pos  = wordTbl.wordIdx(rows);
    wordTbl.isTarget(rows)    = (pos==tIdx);
    wordTbl.isSpillover(rows) = (pos==(tIdx+1));
end

fprintf("[P7 target] isTarget=%d | isSpillover=%d\n", sum(wordTbl.isTarget), sum(wordTbl.isSpillover));

tRows = wordTbl.isTarget | wordTbl.isSpillover;
fprintf("[CHECK target+spill] rows=%d | nFix>0=%.3f | skipped=%.3f\n", ...
    sum(tRows), mean(wordTbl.nFix(tRows)>0), mean(wordTbl.skipped(tRows)));

%% ===== [P8] goPast (regression path duration) =====
fx = subj.event.fix;
fxWord   = [fx.word]';
fxOnset  = [fx.onset]';
fxOffset = [fx.offset]';

wordTbl.goPast = nan(height(wordTbl),1);

for i=1:height(wordTbl)
    tr = wordTbl.trial(i);
    wi = wordTbl.wordIdx(i);
    if ~isfinite(tr) || tr<1 || tr>numel(subj.trial), continue; end
    if ~isfinite(wi) || wi<=0, continue; end

    fixIdx = subj.trial(tr).fixIdx(:);
    fixIdx = fixIdx(fixIdx>0 & fixIdx<=numel(fxWord));
    if isempty(fixIdx), continue; end

    onThis = fixIdx(fxWord(fixIdx)==wi);
    if isempty(onThis), continue; end

    [~, loc] = ismember(onThis, fixIdx);
    loc(loc==0) = [];
    if isempty(loc), continue; end

    firstPos = min(loc);
    tStart   = fxOnset(fixIdx(firstPos));

    exitPos = NaN;
    for k = firstPos+1:numel(fixIdx)
        if fxWord(fixIdx(k)) > wi
            exitPos = k;
            break;
        end
    end

    if isnan(exitPos)
        tEnd = fxOffset(fixIdx(end));
    else
        tEnd = fxOffset(fixIdx(exitPos));
    end

    gp = tEnd - tStart;
    if isfinite(gp) && gp>8000, gp = NaN; end
    wordTbl.goPast(i) = gp;
end

fprintf("[P8 goPast] finite=%d/%d | median=%.0f\n", ...
    sum(isfinite(wordTbl.goPast)), height(wordTbl), median(wordTbl.goPast,'omitnan'));

%% ===== [P9] LANDING-ONLY FILTER SETS (create AFTER target & goPast) =====
wL  = wordTbl(wordTbl.nFix>0 & isfinite(wordTbl.landingNorm), :);
wL2 = wordTbl(wordTbl.nFix>0 & isfinite(wordTbl.landingNorm) & wordTbl.GD>=80, :);

fprintf("[P9 wL]  rows=%d (landing finite)\n", height(wL));
fprintf("[P9 wL2] rows=%d (landing finite & GD>=80)\n", height(wL2));

%% ===== [P10] EXPORT (optional) =====
assignin('base','subj',subj);
assignin('base','wordTbl',wordTbl);
assignin('base','wL',wL);
assignin('base','wL2',wL2);
assignin('base','paddingPx',paddingPx);

%% 8. summary 기반 참가자/시험 퀄리티 필터링 및 최종 저장

% ---- 8.0 설정값(원하면 나중에 조정 가능) ----
qcCfg = struct();
qcCfg.maxSkipRate        = 0.70;   % 단어 skip 비율이 70% 넘으면 이상 trial
qcCfg.minReadingTimeMs   = 300;    % 너무 짧은 읽기(문장 전체 TVT < 300ms)는 이상
qcCfg.maxReadingTimeMs   = 8000;   % 너무 긴 읽기(> 8초)는 이상
qcCfg.minComprehension   = 0.75;   % 참가자 정답률 75% 미만이면 제외 후보
qcCfg.minUsableTrialProp = 0.50;   % usable trial 비율이 50% 미만이면 제외 후보

% subj2res is already defined in [ANCHOR PATCH-CORE v1] as 1:nMain mapping
% ---- 8.1 trial-level summary 만들기 (wordTbl + results 기반) ----
trials = unique(wordTbl.trial);
mainSet = mainTrialIdx;
trials = trials(ismember(trials, mainSet));
nT = numel(trials);

trialSummary = table( ...
    'Size',[nT 9], ...
    'VariableTypes', {'string','double','double','double','double','double','double','double','logical'}, ...
    'VariableNames', {'subj','trial','nWords','nSkipped','skipRate','readingTimeMs','meanFFD','meanGD','acc'});

for ii = 1:nT
    t = trials(ii);

    % 이 trial의 word-level rows
    rows_t = wordTbl(wordTbl.trial == t, :);

    % 단어 수 / 스킵 수
    nWords   = height(rows_t);
    nSkipped = sum(rows_t.skipped);

    % 문장 전체 TVT = 해당 trial의 모든 단어 TVT 합
    tvtAll   = nansum(rows_t.TVT);

    % FFD / GD 평균 (단어들 중 NaN 아닌 것만)
    meanFFD = nanmean(rows_t.FFD);
    meanGD  = nanmean(rows_t.GD);

    % behavioral accuracy (mainscript results에서 가져오기)
    %   - Practice가 아닌 본 시행 기준: results.acc(i) 가 trial i 정오답 (논리)
    % subj trial -> results(main) index

    tw = subj2res(t);  % results index (1..224)
    assert(isfinite(tw), "Non-main trial entered: t=%d", t);

    if isfield(results,'acc') && isfinite(tw) && tw>=1 && tw<=numel(results.acc)
        acc = logical(results.acc(tw));
    else
        acc = true; % fallback
    end

    % subj ID (parseAscToStruct에서 subj.id에 넣어둔 파일명 사용)
    if isfield(subj, 'id')
        sid = string(subj.id);
    else
        sid = "PILOT3";   % fallback
    end

    trialSummary.subj(ii)          = sid;
    trialSummary.trial(ii)         = t;
    trialSummary.nWords(ii)        = nWords;
    trialSummary.nSkipped(ii)      = nSkipped;
    trialSummary.skipRate(ii)      = nSkipped / max(nWords,1);
    trialSummary.readingTimeMs(ii) = tvtAll;       % TVT 기반 읽기시간
    trialSummary.meanFFD(ii)       = meanFFD;
    trialSummary.meanGD(ii)        = meanGD;
    trialSummary.acc(ii)           = acc;
end

% ---- 8.2 trial-level 이상치 플래그 ----
badTrial = false(height(trialSummary),1);

% (1) skip 비율이 너무 높으면 이상 trial
badTrial = badTrial | trialSummary.skipRate > qcCfg.maxSkipRate;

% (2) 읽기 시간이 너무 짧거나 긴 trial
badTrial = badTrial | trialSummary.readingTimeMs < qcCfg.minReadingTimeMs;
badTrial = badTrial | trialSummary.readingTimeMs > qcCfg.maxReadingTimeMs;

% (3) 오답 trial 제외
badTrial = badTrial | ~trialSummary.acc;

trialSummary.badTrial = badTrial;

% ---- 8.3 participant-level 퀄리티 체크 (지금은 한 명 기준) ----
% 정답률 (badTrial 여부와 상관없이 전체 기준)
overallAcc = mean(double(trialSummary.acc));

% usable trial 비율 (badTrial이 아닌 trial 비율)
usableProp = mean(~trialSummary.badTrial);

badSubj = (overallAcc < qcCfg.minComprehension) || ...
          (usableProp < qcCfg.minUsableTrialProp);

fprintf('Participant %s: acc=%.2f, usableTrialProp=%.2f, badSubj=%d\n', ...
    trialSummary.subj(1), overallAcc, usableProp, badSubj);

% ---- 8.4 최종 clean summary 만들기 ----
if badSubj
    cleanSummary = trialSummary([],:);   % 전부 제외 (빈 table)
    warning('Participant %s flagged as badSubj → cleanSummary is empty.', trialSummary.subj(1));
else
    cleanSummary = trialSummary(~trialSummary.badTrial, :);
end

% ---- 8.5 최종 저장 ----
outFile = fullfile(baseDir, subjDir, 'eyeReading_cleanSummary.mat');
save(outFile, 'cleanSummary', 'trialSummary', 'qcCfg');

fprintf('[SAVE] cleanSummary saved to %s (nTrials=%d, nClean=%d)\n', ...
    outFile, height(trialSummary), height(cleanSummary));

%% A. 한 trial에서 time → word index scanpath 보기
k = 1;
t = mainTrialIdx(k);

% 1) 이 trial의 word-level 데이터만 추출
rows_t = wordTbl(wordTbl.trial == t, :);

% 2) 이 trial에서 실제로 fixation이 있었던 단어만 골라서,
%    그 단어들의 FFD start time이 있다고 가정 (makeWordTrialTable에서 필요하면 추가 가능)
%    일단은 TVT가 0이 아닌 단어들만 “읽힌 단어”로 간주하고, 
%    가짜 time index를 만든 버전:

readRows = rows_t(~isnan(rows_t.TVT) & rows_t.TVT > 0, :);

if isempty(readRows)
    warning('Trial %d: 읽힌 단어가 없습니다.', t);
else
    figure;
    % x축: 단어 순서대로 가짜 time index (1,2,3,…)
    % y축: wordIdx
    plot(1:height(readRows), readRows.wordIdx, '-o');
    xlabel('Fixation order (approx)');
    ylabel('Word index');
    title(sprintf('Trial %d: approx scanpath (word index)', t));
    grid on;
end

outWordFile = fullfile(baseDir, subjDir, 'wordTbl_PILOT3.mat');
save(outWordFile, 'wordTbl');

% === subj / design 정보도 따로 저장 (Stats용) ===
subjFile = fullfile(baseDir, subjDir, 'subj.mat');
save(subjFile, 'subj', 'designRow', 'Tmain');
fprintf('[SAVE] subj saved to %s\n', subjFile);