function subj = cleanFixationDurations(subj, shortThresh, longThresh)
% cleanFixationDurations  fixation duration 기반 클리닝
%
%   subj = cleanFixationDurations(subj, shortThresh, longThresh)
%
%   - shortThresh (ms): dur < shortThresh → 너무 짧은 fixation
%   - longThresh  (ms): dur > longThresh  → 너무 긴 fixation (flag만, 제거 X)
%
%   전제:
%       subj.event.fix(k).dur (ms 단위)
%       subj.trial(t).fixIdx = 이 trial의 fixation 인덱스들
%       (있다면) subj.event.fix(k).word = 단어 인덱스 (0 = unassigned)
%
%   결과:
%       subj.event.fix(k).isShortFix
%       subj.event.fix(k).isLongFix
%       subj.event.fix(k).isRemovedShort
%       subj.event.fix(k).mergedTo      (0이면 merge 안 됨)
%
%       subj.trial(t).fixIdxDurClean    (duration 기준 클린 인덱스)
%       subj.trial(t).longFixIdx        (긴 fixation 인덱스)

    if nargin < 2 || isempty(shortThresh)
        shortThresh = 80;   % 예: 60–80ms
    end
    if nargin < 3 || isempty(longThresh)
        longThresh = 800;   % 예: 800–1200ms
    end

    % (D) 작은 안전장치: fixation 이벤트가 없으면 바로 종료
    if ~isfield(subj, 'event') || ~isfield(subj.event, 'fix') || isempty(subj.event.fix)
        warning('cleanFixationDurations: no fixation events found. Skipping.');
        return;
    end

    % ----- 0) event 레벨 flag 초기화 -----
    nFix = numel(subj.event.fix);
    allDur = [subj.event.fix.dur];

    isShort = allDur < shortThresh;
    isLong  = allDur > longThresh;

    hasWord = isfield(subj.event.fix, 'word');

    for k = 1:nFix
        subj.event.fix(k).isShortFix    = isShort(k);
        subj.event.fix(k).isLongFix     = isLong(k);
        subj.event.fix(k).isRemovedShort = false;
        subj.event.fix(k).mergedTo      = 0;
    end

    % ----- 1) trial별로 짧은 fixation merge / 제거 -----
    nTrials = numel(subj.trial);

    for t = 1:nTrials
        idx = subj.trial(t).fixIdx;
        if isempty(idx)
            subj.trial(t).fixIdxDurClean = [];
            subj.trial(t).longFixIdx     = [];
            continue;
        end

        % 이 trial에 속한 fixation들
        fixIdx = idx(:)';                 % row vector
        fix    = subj.event.fix(fixIdx);  % struct array
        dur    = [fix.dur];

        if hasWord
            words = [fix.word];           % 같은 단어 여부 확인용
        else
            words = nan(size(dur));       % word 정보 없으면 merge 안 함
        end

        keep = true(size(fixIdx));        % 기본은 keep

        % 긴 fixation은 제거하지 않고 trial 정보에만 따로 저장
        longMask = dur > longThresh;
        subj.trial(t).longFixIdx = fixIdx(longMask);

        % --- 짧은 fixation 처리 ---
        for k = 1:numel(fixIdx)
            if dur(k) >= shortThresh
                continue;  % 짧지 않으면 패스
            end

            w = words(k);

            merged = false;

            % word 정보가 있고, 0이 아니면 merge 후보 찾기
            if ~isnan(w) && w > 0
                % 1) 바로 앞 fixation과 같은 단어이면 앞 fixation으로 merge
                if k > 1 && words(k-1) == w && dur(k-1) >= shortThresh
                    baseIdx = fixIdx(k-1);
                    curIdx  = fixIdx(k);

                    base = subj.event.fix(baseIdx);
                    cur  = subj.event.fix(curIdx);

                    base.onset  = min(base.onset,  cur.onset);
                    base.offset = max(base.offset, cur.offset);
                    base.dur    = base.offset - base.onset;

                    subj.event.fix(baseIdx) = base;

                    keep(k) = false;
                    subj.event.fix(curIdx).isRemovedShort = true;
                    subj.event.fix(curIdx).mergedTo       = baseIdx;
                    merged = true;
                % 2) 아니면, 바로 뒤 fixation과 같은 단어이면 뒤 fixation으로 merge
                elseif k < numel(fixIdx) && words(k+1) == w && dur(k+1) >= shortThresh
                    baseIdx = fixIdx(k+1);
                    curIdx  = fixIdx(k);

                    base = subj.event.fix(baseIdx);
                    cur  = subj.event.fix(curIdx);

                    base.onset  = min(base.onset,  cur.onset);
                    base.offset = max(base.offset, cur.offset);
                    base.dur    = base.offset - base.onset;

                    subj.event.fix(baseIdx) = base;

                    keep(k) = false;
                    subj.event.fix(curIdx).isRemovedShort = true;
                    subj.event.fix(curIdx).mergedTo       = baseIdx;
                    merged = true;
                end
            end

            % 3) merge 못 했으면 그냥 제거
            if ~merged
                keep(k) = false;
                subj.event.fix(fixIdx(k)).isRemovedShort = true;
            end
        end

        % --- 이 trial에서 duration 클리닝 후 사용할 fixation 인덱스 ---
        subj.trial(t).fixIdxDurClean = fixIdx(keep);
    end
end

%% 분석 결과
% 모든 trial의 원래 fixation index와 클린 후 index를 연결
% allOrig  = [];
% allClean = [];

% for t = 1:numel(subj.trial)
%     if isfield(subj.trial(t), 'fixIdx') && ~isempty(subj.trial(t).fixIdx)
%         allOrig  = [allOrig  subj.trial(t).fixIdx(:)'];          %#ok<AGROW>
%     end
%     if isfield(subj.trial(t), 'fixIdxDurClean') && ~isempty(subj.trial(t).fixIdxDurClean)
%         allClean = [allClean subj.trial(t).fixIdxDurClean(:)'];  %#ok<AGROW>
%     end
% end

% nOrig  = numel(allOrig);
% nClean = numel(allClean);

% propKept    = nClean / nOrig;       % 남은 비율 (유효)
% propRemoved = 1 - propKept;         % 잘려나간 비율

% [ nOrig, nClean, propKept, propRemoved ]

% ans =

%    1.0e+03 *

%     4.1420    3.8510    0.0009    0.0001

% nOrig ≈ 4142: 원래 trial별 fixIdx를 다 모았을 때 fixation 인덱스 총 개수
% nClean ≈ 3851:duration 클리닝 후 남은 fixation 개수 (fixIdxDurClean 기준)
% propKept ≈ 0.93: 유지된 비율 ≈ 92.97%
% propRemoved ≈ 0.07: 잘려나간 비율 ≈ 7.03%