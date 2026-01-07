function wordTbl = makeWordTrialTable(subj, results, trialIdx)
% trialIdx 없으면 전체 trial
% makeWordTrialTable  trial × word 수준의 읽기 지표(FFD, GD, TVT 등) 계산
%
%   wordTbl = makeWordTrialTable(subj, results)
%
%   전제:
%       - subj.trial(t).fixIdxDurClean : duration 클리닝 후 fixation 인덱스
%       - subj.event.fix(k).dur        : fixation duration (ms)
%       - subj.event.fix(k).onset      : fixation onset (절대 시간, ms)
%       - subj.event.fix(k).word       : word ROI 매핑 결과 (0 = unassigned)
%       - subj.trial(t).startTime      : trial 시작 시간 (절대 시간, ms)
%       - results.wordRects{t}         : trial t의 단어 ROI [nWords×4] (L T R B)
%       - (선택) results.words{t}      : trial t의 단어 문자열 리스트
%
%   출력(wordTbl 예시 컬럼):
%       trial           : trial 번호
%       wordIdx         : 문장 내 단어 인덱스(1..nWords)
%       wordStr         : 단어 문자열(있으면), 없으면 ""
%       nFix            : 해당 단어에 떨어진 fixation 개수
%       skipped         : 한 번도 fix 안 했으면 true
%       FFD             : First Fixation Duration (ms)
%       GD              : Gaze Duration (first-pass; ms)
%       TVT             : Total Viewing Time (ms)
%       firstFixOnset   : trial 시작 기준 첫 fixation onset (ms)

% ---- trialIdx 처리 (없으면 전체 trial) ----
if nargin < 3 || isempty(trialIdx)
    trialIdx = 1:numel(subj.trial);
end
trialIdx = unique(trialIdx(:)');                          % 중복 제거 + row
trialIdx = trialIdx(trialIdx>=1 & trialIdx<=numel(subj.trial));  % 범위 방어

rows = [];  % struct 배열로 모은 뒤 마지막에 table로 변환

% ---- trial loop: trialIdx만 사용 ----
for t = trialIdx
    % 이 trial에 ROI 정보가 없다면 스킵
    if ~isfield(results,'wordRects') || t > numel(results.wordRects) || isempty(results.wordRects{t})
        continue;
    end

    rects  = results.wordRects{t};
    nWords = size(rects,1);
    tr     = subj.trial(t);   % trial 정보 (startTime 포함)

    % duration 클리닝 반영된 fixation 인덱스 사용
    if ~isfield(subj.trial(t), 'fixIdxDurClean') || isempty(subj.trial(t).fixIdxDurClean)
        if isfield(subj.trial(t), 'fixIdx') && ~isempty(subj.trial(t).fixIdx)
            fixIdx = subj.trial(t).fixIdx(:)';
        else
            fixIdx = [];
        end
    else
        fixIdx = subj.trial(t).fixIdxDurClean(:)';
    end

        % (1) 이 trial에서 사용할 fixation index 선택
        if isempty(fixIdx)
            % fixation이 하나도 없는 trial이면 모든 단어를 skipped로 기록
            for w = 1:nWords
                row = makeEmptyWordRow(t, w, results);
                rows = [rows; row]; %#ok<AGROW>
            end
            continue;
        end
        
        % (2) 이 trial의 fixation struct 뽑기
        fix = subj.event.fix(fixIdx);

        % isValid 필드가 있다면, 여기서 한 번 더 필터
        if isfield(fix, 'isValid')
            validMask = [fix.isValid];
            fix    = fix(validMask);
            fixIdx = fixIdx(validMask);
        end

        % isValid까지 적용하고 나니까 이 trial에서 쓸 fixation이 0개일 수도 있음
    if isempty(fix)
        % 이 경우도 결국 "fixation이 하나도 없는 trial" 이라고 보고
        % 모든 단어를 skipped 로만 기록
        for w = 1:nWords
            row = makeEmptyWordRow(t, w, results);
            rows = [rows; row]; %#ok<AGROW>
        end
        continue;
    end

        % (3) onsetRel 계산 → trial 시작 기준 상대 시간(ms)으로 변환
         if isfield(fix, 'onset')
            onsetAbs = [fix.onset];
            [onsetAbs, sortIdx] = sort(onsetAbs);  % 시간 순서대로 정렬

            % fix, fixIdx 둘 다 같은 순서로 재정렬
            fix    = fix(sortIdx);
            fixIdx = fixIdx(sortIdx);
        else
            onsetAbs = nan(size(fixIdx));
        end

        dur = [fix.dur];

        % onset (절대 시간) → trial 시작 기준 상대 시간(ms)
        if isfield(tr, 'startTime') && ~isempty(tr.startTime)
            onsetRel = onsetAbs - tr.startTime;
        else
            onsetRel = onsetAbs;   % 혹시 startTime이 없으면 절대 시간 그대로
        end
        
        % (4) word index 가져오기
        if isfield(fix, 'word')
            wIdx = [fix.word];
        else
            wIdx = zeros(size(dur));  % word 정보 없으면 0으로
        end

        % 각 단어에 대해 지표 계산
        for w = 1:nWords
            wordFixMask      = (wIdx == w);
            thisFixIdxLocal  = find(wordFixMask);   % trial 내 fix 배열 기준 인덱스

            row = makeEmptyWordRow(t, w, results);

            if isempty(thisFixIdxLocal)
                % 한 번도 fix 안 함 → skipped = true (기본값 유지)
                % firstFixOnset은 NaN 그대로 둠
            else
                row.skipped = false;
                row.nFix    = numel(thisFixIdxLocal);

                % --- TVT (Total Viewing Time) ---
                row.TVT = sum(dur(thisFixIdxLocal));

                % --- FFD (First Fixation Duration) ---
                firstIdx = thisFixIdxLocal(1);
                row.FFD  = dur(firstIdx);

                % --- GD (Gaze Duration: first-pass) ---
                %   첫 진입 후 연속으로 머무는 구간 합
                runEnd = firstIdx;
                while runEnd < numel(fix) && wIdx(runEnd+1) == w
                    runEnd = runEnd + 1;
                end
                row.GD = sum(dur(firstIdx:runEnd));

                % --- firstFixOnset (trial 시작 기준 단어 첫 fixation onset) ---
                if all(isnan(onsetRel))
                    row.firstFixOnset = NaN;
                else
                    % 이 단어에 대해 첫 fixation의 onset (상대시간)
                    row.firstFixOnset = onsetRel(firstIdx);
                end
            end

            rows = [rows; row]; %#ok<AGROW>
        end
    end

    % struct 배열 → table
    if isempty(rows)
        wordTbl = table();
    else
        wordTbl = struct2table(rows);
    end
end

% ------------------------------------------------------------
function row = makeEmptyWordRow(t, w, results)
% trial/word에 대한 기본 값 채우기

    row.trial   = t;
    row.wordIdx = w;

    % word 문자열이 있으면 함께 저장
    if isfield(results, 'words') && numel(results.words) >= t ...
            && ~isempty(results.words{t}) && numel(results.words{t}) >= w
        ws = results.words{t}{w};
        if isstring(ws), ws = char(ws); end
        row.wordStr = string(ws);
    else
        row.wordStr = "";
    end

    row.nFix          = 0;
    row.skipped       = true;
    row.FFD           = NaN;
    row.GD            = NaN;
    row.TVT           = NaN;
    row.firstFixOnset = NaN;   % 기본값은 NaN
end


