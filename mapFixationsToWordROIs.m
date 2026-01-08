function subj = mapFixationsToWordROIs(subj, wordRectsCell, paddingPx)
% mapFixationsToWordROIs
%   각 trial의 단어 ROI(사각형)에 fixation을 할당해서
%   subj.event.fix(k).word 를 채워주는 함수.
%
%   subj = mapFixationsToWordROIs(subj, wordRectsCell, paddingPx)
%
%   입력:
%     - subj: ASC 파싱 + trial 분류 + (선택) drift / 클리닝까지 끝난 struct
%     - wordRectsCell{t}: t번째 trial의 [nWords × 4] (L T R B) 좌표
%         * mainscript에서 만든 results.wordRects랑 동일 형식 가정
%     - paddingPx: 경계 여유 픽셀 (예: 2~3 픽셀)
%
%   전제:
%     - subj.trial(t).fixIdx 또는 subj.trial(t).fixIdxDurClean 이 있음
%     - subj.event.fix(k).xCorr / .yCorr 가 있으면 그걸 우선 사용, 없으면 x / y 사용
%
%   결과:
%     - subj.event.fix(k).word : 0 = 아무 단어에도 안 들어감, 1..nWords
%     - subj.trial(t).wordRects : padding 적용된 [nWords × 4]
%     - subj.trial(t).wordFixIdx{w} : 단어 w에 속한 fixation 인덱스 리스트

    if nargin < 3 || isempty(paddingPx)
        paddingPx = 0;   % 기본은 padding 없음
    end

    % --- 0) 모든 fixation의 word 필드를 0으로 초기화 ---
    nFixAll = numel(subj.event.fix);
    for k = 1:nFixAll
        subj.event.fix(k).word = 0;
    end

    nTrials = numel(subj.trial);

    % --- 1) trial별로 ROI + fixation 매핑 ---
    for t = 1:nTrials

        % wordRectsCell 보다 trial 수가 많으면 남은 trial은 스킵
        if t > numel(wordRectsCell) || isempty(wordRectsCell{t})
            subj.trial(t).wordRects   = [];
            subj.trial(t).wordFixIdx  = {};
            subj.trial(t).nWords      = 0;
            continue;
        end

        rects = wordRectsCell{t};   % [nWords × 4], [L T R B]
        if size(rects,2) ~= 4
            error('wordRectsCell{%d}는 [nWords x 4] (L T R B) 여야 합니다.', t);
        end

        % --- padding 적용 ---
        L = rects(:,1) - paddingPx;
        T = rects(:,2) - paddingPx;
        R = rects(:,3) + paddingPx;
        B = rects(:,4) + paddingPx;
        rectsPad = [L T R B];

        subj.trial(t).wordRects  = rectsPad;
        nWords                   = size(rectsPad,1);
        subj.trial(t).nWords     = nWords;
        subj.trial(t).wordFixIdx = cell(nWords,1);  % 나중에 채움

        % --- 이 trial에서 사용할 fixation 인덱스 선택 ---
        if isfield(subj.trial(t), 'fixIdxDurClean') && ~isempty(subj.trial(t).fixIdxDurClean)
            fixIdx = subj.trial(t).fixIdxDurClean(:)';  % duration 클리닝 반영
        else
            fixIdx = subj.trial(t).fixIdx(:)';          % 기본 fixIdx
        end

        if isempty(fixIdx)
            continue;
        end

        % (B) 안전장치: fixation index가 subj.event.fix 범위를 벗어나지 않는지 확인
        maxFixIdx = numel(subj.event.fix);
        if any(fixIdx < 1 | fixIdx > maxFixIdx)
            error('mapFixationsToWordROIs: trial %d has invalid fixation index (max=%d).', ...
                  t, maxFixIdx);
        end

        % --- fixation 하나씩 word ROI에 매핑 ---
        for k = 1:numel(fixIdx)
            fi = fixIdx(k);

            % 1) drift 보정 좌표가 있으면 그걸 우선 사용, 없으면 raw 사용
            if isfield(subj.event.fix(fi), 'xCorr')
                fx = subj.event.fix(fi).xCorr;
                fy = subj.event.fix(fi).yCorr;
            else
                fx = subj.event.fix(fi).x;
                fy = subj.event.fix(fi).y;
            end

            % 2) fixation validity (cleanFixations 결과)를 존중해서
            %    invalid fixation이면 어떤 단어에도 매핑하지 않음
            if isfield(subj.event.fix(fi), 'isValid') && ~subj.event.fix(fi).isValid
                subj.event.fix(fi).word = 0;
                continue;
            end

            % 3) 좌표 NaN이면 매핑하지 않음
            if isnan(fx) || isnan(fy)
                subj.event.fix(fi).word = 0;
                continue;
            end

            % 4) 어떤 단어 ROI 안에 들어가는지 검사 (L <= x <= R && T <= y <= B)
            inX = (fx >= rectsPad(:,1)) & (fx <= rectsPad(:,3));
            inY = (fy >= rectsPad(:,2)) & (fy <= rectsPad(:,4));
            inWord = find(inX & inY);

            if isempty(inWord)
                % 어떤 단어에도 안 들어감
                subj.event.fix(fi).word = 0;
            else
                if numel(inWord)==1
                    w = inWord;
                else
                    cx = (rectsPad(inWord,1) + rectsPad(inWord,3)) / 2;
                    cy = (rectsPad(inWord,2) + rectsPad(inWord,4)) / 2;
                    d2 = (fx - cx).^2 + (fy - cy).^2;
                    [~,ix] = min(d2);   % 동률이면 첫 후보(작은 w)
                    w = inWord(ix);
                end

                subj.event.fix(fi).word = w;
                subj.trial(t).wordFixIdx{w}(end+1,1) = fi;
            end
        end
    end
end
