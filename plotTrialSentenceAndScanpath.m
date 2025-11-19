function plotTrialSentenceAndScanpath(subj, results, t, useCorrected)
% plotTrialSentenceAndScanpath
%   한 trial의 문장(ROI) 위에 fixation 경로를 그려주는 함수
%
%   plotTrialSentenceAndScanpath(subj, results, t)
%   plotTrialSentenceAndScanpath(subj, results, t, useCorrected)
%
%   useCorrected = true  → drift 보정 좌표(xCorr,yCorr) 사용(있을 때)
%                  false → raw 좌표(x,y) 사용

    if nargin < 4 || isempty(useCorrected)
        useCorrected = false;
    end

    % ---------- 1. ROI / 단어 정보 ----------
    if t > numel(results.wordRects) || isempty(results.wordRects{t})
        error('results.wordRects{%d}가 비어 있습니다.', t);
    end

    rects = results.wordRects{t};   % [nWords x 4] (L T R B)
    nWords = size(rects,1);

    if isfield(results, 'words') && numel(results.words) >= t ...
            && ~isempty(results.words{t})
        words = results.words{t};
    else
        words = arrayfun(@(k) sprintf('W%d',k), 1:nWords, 'UniformOutput', false);
    end

    % ---------- 2. 이 trial에서 사용할 fixation ----------
    if isfield(subj.trial, 'fixIdxDurClean') && ~isempty(subj.trial(t).fixIdxDurClean)
        fixIdx = subj.trial(t).fixIdxDurClean(:)';
    elseif isfield(subj.trial, 'fixIdxValid') && ~isempty(subj.trial(t).fixIdxValid)
        fixIdx = subj.trial(t).fixIdxValid(:)';
    else
        fixIdx = subj.trial(t).fixIdx(:)';
    end

    if isempty(fixIdx)
        warning('Trial %d: fixation이 없습니다.', t);
        return;
    end

    F = subj.event.fix(fixIdx);

    hasCorr = isfield(F, 'xCorr') && isfield(F, 'yCorr');
    if useCorrected && hasCorr
        gx = [F.xCorr];
        gy = [F.yCorr];
        gazeType = 'Corrected';
    else
        gx = [F.x];
        gy = [F.y];
        gazeType = 'Raw';
    end

    % ---------- 3. 그림 ----------
    figure; hold on;
    title(sprintf('Trial %d scanpath on sentence (%s gaze)', t, gazeType));
    xlabel('X (pixels)');
    ylabel('Y (pixels)');
    set(gca,'YDir','reverse');   % 화면 좌표계와 맞추기

    xlim([0 1920]);
    ylim([0 1080]);
    grid on;

    % (1) 단어 ROI와 텍스트를 실제 위치에 그리기
for w = 1:nWords
    L = rects(w,1); T = rects(w,2);
    R = rects(w,3); B = rects(w,4);
    wWidth  = R - L;
    wHeight = B - T;

    % ROI 박스
    rectangle('Position',[L T wWidth wHeight], ...
              'EdgeColor',[0.8 0.8 0.8]);

    % 단어 텍스트 (ROI 중앙)
    cx = (L + R)/2;
    cy = (T + B)/2;
    text(cx, cy, words{w}, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontSize',18, 'Color','k');
end

    % (2) 단어 중심 위치 계산 (ROI 중앙)
    cx = (rects(:,1) + rects(:,3)) / 2;
    cy = (rects(:,2) + rects(:,4)) / 2;

    % (3) scanpath: fixation 위치와 순서 표시
    plot(gx, gy, '-o', 'LineWidth',1.5);

    for k = 1:numel(F)
        text(gx(k), gy(k)-10, sprintf('%d',k), ...
            'HorizontalAlignment','center', ...
            'Color','b');
    end

    hold off;
end
