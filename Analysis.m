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
subjDir = '20251028_134927_Pilot3_교수님';
ascFileName = 'PILOT3.asc';
ascFile = fullfile(baseDir, subjDir, ascFileName);

disp(ascFile)   % 실제 경로 확인용
if ~isfile(ascFile)
    error('ASC 파일을 찾을 수 없습니다: %s', ascFile);
end

% --- mainscript 결과(ROI 정보) 불러오기 ---
elFile = fullfile(baseDir, subjDir, 'EL_demo.mat');
if ~isfile(elFile)
    error('EL_demo.mat not found: %s', elFile);
end

load(elFile, 'results', 'dp');
% results.wordRects{t} : 각 trial의 [nWords×4] (L T R B) ROI
% results.words{t}     : 각 trial의 단어 string 리스트(있다면)

% 2. ASC → MATLAB struct (sample + event)
%   → parseAscToStruct 내부에서 sample, fix/sacc/blink, MSG, START/END까지 다 파싱해서
%     subj.sample / subj.event.*에 정리함.
subj = parseAscToStruct(ascFile);

% 2.1 샘플 validity 플래그
subj = addSampleValidity(subj, [0 0 1920 1080], 50);

% 2.2 MSG 기반 trial 경계 정의 및 trial struct 만들기
subj = addTrialsFromMsg(subj);

%% 3. drift 기준점(xTrue, yTrue) 자동 계산 (main trial 첫 단어 평균)

% main trial 인덱스: id가 'TRIALID'로 시작하는 trial만 사용
isMain  = arrayfun(@(tr) startsWith(tr.id, 'TRIALID'), subj.trial);
mainIdx = find(isMain);

nMain = numel(mainIdx);
xList = NaN(nMain,1);
yList = NaN(nMain,1);

for k = 1:nMain
    t = mainIdx(k);               % subj.trial에서의 trial index
    if t > numel(results.wordRects)
        continue;
    end

    rects = results.wordRects{t}; % [nWords×4] (L T R B)
    if isempty(rects)
        continue;
    end

    rect = rects(1,:);            % 첫 단어 ROI [x1 y1 x2 y2]
    xList(k) = mean(rect([1 3])); % (x1+x3)/2
    yList(k) = mean(rect([2 4])); % (y1+y3)/2
end

xTrue = mean(xList, 'omitnan');
yTrue = mean(yList, 'omitnan');

fprintf('Drift 기준점(xTrue, yTrue) = (%.2f, %.2f)\n', xTrue, yTrue);

% 3. 좌표 sanity check 및 drift correction (문장 기준 좌표 보정)
%   % ---- drift correction 적용 ----
% 문장 시작점(첫 단어 or 왼쪽 fixation marker)의 이상적인 좌표를 지정.
% 실험 세팅에 맞게 xTrue/yTrue 바꿔줘야 함.
driftCfg = struct();
driftCfg.xTrue         = xTrue;   
driftCfg.yTrue         = yTrue;   
driftCfg.refEvent      = 'SENTENCE_ONSET';
driftCfg.minFixDur     = 60;      % ms
driftCfg.maxAbsOffsetX = 60;      % 허용 x 오차 (px)
driftCfg.maxAbsOffsetY = 40;      % 허용 y 오차 (px)

subj = applyDriftCorrection(subj, driftCfg);

% 이제 subj.sample.gxCorr / gyCorr, subj.event.fix.xCorr / yCorr,
% subj.trial(t).driftOffset / isBadDrift 등이 채워져 있음.
%
% ---- sanity check (시각 확인) ----
% 몇 개 trial을 골라 raw vs corrected time–gx 플롯을 비교해본다.
trialsToCheck = [1 5 10];   % 보고 싶은 trial 번호로 수정

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

% 5. fixation → word ROI 매핑
paddingPx = 2;
subj = mapFixationsToWordROIs(subj, results.wordRects, paddingPx);

% 6. Fixation duration 기반 클리닝
shortThresh = 80;
longThresh  = 800;
subj = cleanFixationDurations(subj, shortThresh, longThresh);

% 7. word × trial 지표 계산 (FFD, GD, TVT, skip, regressions 등)
wordTbl = makeWordTrialTable(subj, results);

% 결과 확인 예시
head(wordTbl)          % 앞 몇 줄 눈으로 확인
tvt = wordTbl.TVT;
fprintf('TVT: n=%d, mean=%.1f, sd=%.1f, min=%.1f, max=%.1f\n', ...
    sum(~isnan(tvt)), mean(tvt,'omitnan'), std(tvt,'omitnan'), ...
    min(tvt,[],'omitnan'), max(tvt,[],'omitnan'));

figure;
histogram(tvt, 20);
xlabel('TVT (ms)'); ylabel('Count');
title('Distribution of TVT');

%   - 전처리된 fixation·word 매핑을 바탕으로, 단어 w에 대해 대표적인 읽기 지표를 계산:
%       - first fixation duration (FFD):
%           → 해당 단어에서 첫 fixation 하나의 duration
%       - gaze duration (GD):
%           → 단어에 처음 진입한 후, 오른쪽으로 벗어날 때까지
%              그 단어 위 fixation duration의 합
%       - total viewing time (TVT):
%           → 모든 pass를 포함한 그 단어 위 fixation duration의 합
%       - skipping probability:
%           → 그 단어에 fixation이 한 번도 없는 trial 비율
%       - refixation probability:
%           → 한 단어에서 fixation이 2번 이상인 trial 비율
%       - regression path duration(= go-past time):
%           → 단어 w에 처음 들어간 순간부터 다시 그 단어 오른쪽으로 넘어갈 때까지
%             (왼쪽으로 회귀(regr.) 포함) 걸린 시간의 합
%
%   - 예시 summary table 스케치:
%       summary = table(subjID, trialID, wordIndex, ...
%                       FFD, GD, TVT, skip, refixProb, ...
%                       'VariableNames', {...});

%% 8. summary 기반 참가자/시험 퀄리티 필터링 및 최종 저장

% ---- 8.0 설정값(원하면 나중에 조정 가능) ----
qcCfg = struct();
qcCfg.maxSkipRate        = 0.70;   % 단어 skip 비율이 70% 넘으면 이상 trial
qcCfg.minReadingTimeMs   = 300;    % 너무 짧은 읽기(문장 전체 TVT < 300ms)는 이상
qcCfg.maxReadingTimeMs   = 8000;   % 너무 긴 읽기(> 8초)는 이상
qcCfg.minComprehension   = 0.75;   % 참가자 정답률 75% 미만이면 제외 후보
qcCfg.minUsableTrialProp = 0.50;   % usable trial 비율이 50% 미만이면 제외 후보

% ---- 8.1 trial-level summary 만들기 (wordTbl + results 기반) ----
trials = unique(wordTbl.trial);
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
    if isfield(results, 'acc') && numel(results.acc) >= t
        acc = logical(results.acc(t));
    else
        acc = true;  % 없으면 true로 두고 나중에 조정
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

% (3) (옵션) 오답 trial 제외하고 싶으면 아래 줄 주석 해제
% badTrial = badTrial | ~trialSummary.acc;

trialSummary.badTrial = badTrial;

% ---- 8.3 participant-level 퀄리티 체크 (지금은 한 명 기준) ----
% 정답률 (badTrial 여부와 상관없이 전체 기준)
overallAcc = mean(trialSummary.acc);

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
t = 1;   % 보고 싶은 trial 번호

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


