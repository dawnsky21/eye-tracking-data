function subj = cleanFixations(subj, screenRect, lineYRange, blinkMargin, shortDurThresh)
% cleanFixations  blink / out-of-range / duration 기반 fixation 전처리
%
%   subj = cleanFixations(subj, screenRect, lineYRange, blinkMargin, shortDurThresh)
%
%   screenRect     = [xMin yMin xMax yMax], 예: [0 0 1920 1080]
%   lineYRange     = [yLineMin yLineMax], 문장 라인 y 범위 (없으면 [] 또는 생략)
%   blinkMargin    = blink 앞뒤로 몇 ms까지 blink 관련으로 볼지 (기본 50)
%   shortDurThresh = blink 직후 매우 짧은 fixation 기준 (기본 70 ms)
%
%   결과:
%       subj.event.fix(k).isValid
%       subj.event.fix(k).isBlinkRelated
%       subj.event.fix(k).isOutScreen
%       subj.event.fix(k).isOffLine
%       subj.event.fix(k).isShortAfterBlink
%       subj.event.fixValidIdx / fixInvalidIdx

    if nargin < 2 || isempty(screenRect)
        screenRect = [0 0 1920 1080];
    end
    if nargin < 3
        lineYRange = [];
    end
    if nargin < 4 || isempty(blinkMargin)
        blinkMargin = 50;    % ms
    end
    if nargin < 5 || isempty(shortDurThresh)
        shortDurThresh = 70; % ms
    end

        % ---- screenRect / lineYRange sanity check ----
    if numel(screenRect) ~= 4 || any(screenRect(3:4) <= screenRect(1:2))
        error('cleanFixations: screenRect seems invalid. Expected [xMin yMin xMax yMax].');
    end

    if ~isempty(lineYRange)
        if numel(lineYRange) ~= 2 || lineYRange(2) <= lineYRange(1)
            error('cleanFixations: lineYRange seems invalid. Expected [yMin yMax].');
        end
    end

    if ~isfield(subj, 'event') || ~isfield(subj.event, 'fix') || isempty(subj.event.fix)
        warning('cleanFixations: no fixation event found.');
        return;
    end

    fix = subj.event.fix;
    nFix = numel(fix);

    % ---- 좌표: drift 보정된 값이 있으면 그걸 우선 사용 ----
    x = nan(nFix,1);
    y = nan(nFix,1);
    for k = 1:nFix
        if isfield(fix(k), 'xCorr')
            x(k) = fix(k).xCorr;
            y(k) = fix(k).yCorr;
        else
            x(k) = fix(k).x;
            y(k) = fix(k).y;
        end
    end

    onset = [fix.onset]';
    offset = [fix.offset]';
    dur = [fix.dur]';

    % ---- blink 구간 정보 (EBLINK) ----
    blinkStart = [];
    blinkEnd   = [];
    if isfield(subj.event, 'blink') && ~isempty(subj.event.blink)
        blinkStart = [subj.event.blink.onset]';
        blinkEnd   = [subj.event.blink.offset]';
    end

    % blinkMargin 포함한 윈도우
    if ~isempty(blinkStart)
        blinkStartM = blinkStart - blinkMargin;
        blinkEndM   = blinkEnd   + blinkMargin;
    else
        blinkStartM = [];
        blinkEndM   = [];
    end

    % ---- 플래그 초기화 ----
    isValid           = true(nFix,1);
    isBlinkRelated    = false(nFix,1);
    isOutScreen       = false(nFix,1);
    isOffLine         = false(nFix,1);
    isShortAfterBlink = false(nFix,1);

    % =========================================================
    % 4.1 blink 구간 처리
    %   - EBLINK ± blinkMargin과 겹치는 fixation → isBlinkRelated = true, invalid
    %   - blink 직후 매우 짧은 fixation(dur < shortDurThresh) → artifact으로 invalid
    % =========================================================
    if ~isempty(blinkStartM)
        for k = 1:nFix
            % (1) blink 구간과 겹치는지 확인
            overlap = (onset(k) <= blinkEndM) & (offset(k) >= blinkStartM);
            if any(overlap)
                isBlinkRelated(k) = true;
            end
        end

        % (2) blink 직후 very short fixation
        % blink가 끝난 시점과 fixation onset 사이의 시간 차이를 보고 판단
        shortAfterWindow = 100; % ms 이내면 "직후"로 간주 (원하면 조정 가능)
        for k = 1:nFix
            dt = onset(k) - blinkEnd;        % 각 blink 끝 이후 경과 시간
            dt(dt < 0) = inf;                % blink 전에 시작한 건 무시
            minDt = min(dt);
            if minDt >= 0 && minDt <= shortAfterWindow && dur(k) < shortDurThresh
                isShortAfterBlink(k) = true;
            end
        end
    end

    % =========================================================
    % 4.2 화면 밖·비정상 fixation 처리
    %   - screenRect 밖이면 제거
    %   - lineYRange 크게 벗어나면 "off-line"으로 태깅 (기본은 invalid로 두지 않음)
    % =========================================================
    xMin = screenRect(1);  yMin = screenRect(2);
    xMax = screenRect(3);  yMax = screenRect(4);

    outScreen = x < xMin | x > xMax | y < yMin | y > yMax | ...
                isnan(x) | isnan(y);

    isOutScreen(outScreen) = true;

    % 문장 라인의 y 범위를 크게 벗어나는 fixation
    if ~isempty(lineYRange)
        yLineMin = lineYRange(1);
        yLineMax = lineYRange(2);
        offLine = (y < yLineMin) | (y > yLineMax);
        isOffLine(offLine) = true;
        % 전략:
        %   - 여기서는 off-line을 "unassigned 후보"로만 태깅하고,
        %     완전히 제거하지는 않는다 (isValid는 유지).
        %   - 필요하면 아래에서 isValid(offLine) = false; 로 바꿔도 됨.
    end

    %   → 현재 구현: off-line fixation도 나중에 regressions 등 분석에 쓸 수 있게 valid로 남겨둠.

    % =========================================================
    % 최종 validity 결정
    %   - blink 관련 + blink 직후 초단기 + 화면 밖 → invalid
    %   - off-line은 태그만 하고 기본은 valid로 둠
    % =========================================================
    isValid(isBlinkRelated)    = false;
    isValid(isShortAfterBlink) = false;
    isValid(isOutScreen)       = false;

    % 필요하면 한 줄로 바꿀 수 있음:
    % isValid = ~(isBlinkRelated | isShortAfterBlink | isOutScreen);

    % ---- 결과를 subj.event.fix에 저장 ----
    for k = 1:nFix
        subj.event.fix(k).isValid           = isValid(k);
        subj.event.fix(k).isBlinkRelated    = isBlinkRelated(k);
        subj.event.fix(k).isOutScreen       = isOutScreen(k);
        subj.event.fix(k).isOffLine         = isOffLine(k);
        subj.event.fix(k).isShortAfterBlink = isShortAfterBlink(k);
    end

    subj.event.fixValidIdx   = find(isValid);
    subj.event.fixInvalidIdx = find(~isValid);
end

%% 분석 결과
% >> subj = cleanFixations(subj, [0 0 1920 1080], [500 620], 50, 70);

% sum([subj.event.fix.isValid])
% sum([subj.event.fix.isBlinkRelated])
% sum([subj.event.fix.isOutScreen])
% sum([subj.event.fix.isShortAfterBlink])

% ans = 3653 최종적으로 분석에 쓸 수 있는 fixation 개수
% ans = 485  전체 fixation 중 깜빡임(EBLINK ± 50ms) 구간과 겹친 fixation 개수
% ans = 2    화면 밖 또는 NaN 좌표였던 fix 개수
% ans = 24   blink 끝난 직후(100ms 이내) + duration < 70ms인 초단기 fix 개수