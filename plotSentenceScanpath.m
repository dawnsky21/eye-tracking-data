function plotSentenceScanpath(subj, results, t, useCorr)
% plotSentenceScanpath
%   한 trial(t)에 대해
%   - 문장(단어들)을 화면 좌표에 맞게 한 줄로 찍고
%   - 그 위에 fixation scanpath를 겹쳐 보여주는 함수
%
% 사용 예:
%   t = 50;
%   plotSentenceScanpath(subj, results, t, true);   % 보정 좌표(xCorr,yCorr) 사용
%
% 입력
%   subj    : Analysis에서 만든 subj struct
%   results : EL_demo.mat 에서 읽어온 results (wordRects, words 포함)
%   t       : 보고 싶은 trial 번호
%   useCorr : true이면 xCorr/yCorr, false이면 raw x/y 사용

    if nargin < 4
        useCorr = true;
    end

    % ----- 1) trial 존재 & ROI / 단어 정보 확인 -----
    if t > numel(subj.trial)
        error('trial %d 는 subj.trial 범위를 벗어납니다.', t);
    end

    if t > numel(results.wordRects) || isempty(results.wordRects{t})
        error('trial %d 에 wordRects가 없습니다.', t);
    end

    wr = results.wordRects{t};   % [nWords x 4] (L T R B)
    nWords = size(wr,1);

    if isfield(results, 'words') && numel(results.words) >= t ...
            && ~isempty(results.words{t})
        ws = results.words{t};
    else
        % 단어 문자열이 없으면 인덱스 번호만 표시
        ws = arrayfun(@(w) sprintf('%d',w), 1:nWords, 'uni', 0);
    end

    % ----- 2) 이 trial에서 사용할 fixation 인덱스 선택 -----
    if isfield(subj.trial, 'fixIdxDurClean') && ~isempty(subj.trial(t).fixIdxDurClean)
        fixIdx = subj.trial(t).fixIdxDurClean(:)';   % duration 클리닝 반영
    elseif isfield(subj.trial, 'fixIdxValid') && ~isempty(subj.trial(t).fixIdxValid)
        fixIdx = subj.trial(t).fixIdxValid(:)';
    else
        fixIdx = subj.trial(t).fixIdx(:)';
    end

    if isempty(fixIdx)
        warning('trial %d 에 fixation이 없습니다.', t);
        return;
    end

    F = subj.event.fix(fixIdx);

    % 좌표 선택 (보정 좌표가 있으면 그걸 쓸지 여부)
    hasCorr = isfield(F, 'xCorr') && isfield(F, 'yCorr');
    if useCorr && hasCorr
        fx = [F.xCorr];
        fy = [F.yCorr];
        coordLabel = 'Corrected gaze';
    else
        fx = [F.x];
        fy = [F.y];
        coordLabel = 'Raw gaze';
    end

    % duration도 나중에 마커 크기 조절용으로 사용 (선택)
    if isfield(F, 'dur')
        dur = [F.dur];
    else
        dur = ones(size(fx)) * 100;
    end

    % ----- 3) Figure 세팅 -----
    figure; hold on;

    % 화면 해상도(실험에 맞게 수정 가능)
    scrW = 1920;
    scrH = 1080;

    % 축을 화면 좌표처럼: (0,0) = 좌상단이 되도록 Y 방향 뒤집기
    axis([0 scrW 0 scrH]);
    set(gca, 'YDir', 'reverse');

    % ----- 4) 문장(단어) 그리기 -----
    for w = 1:nWords
        rect = wr(w,:);          % [L T R B]
        xCenter = (rect(1) + rect(3)) / 2;
        yCenter = (rect(2) + rect(4)) / 2;

        % 단어 텍스트 표시
        text(xCenter, yCenter, ws{w}, ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontSize', 20);
    end

    % 단어 ROI 박스를 보고 싶으면(옵션) 주석 해제
    % for w = 1:nWords
    %     rectangle('Position', [wr(w,1), wr(w,2), wr(w,3)-wr(w,1), wr(w,4)-wr(w,2)], ...
    %               'EdgeColor',[0.8 0.8 0.8]);
    % end

    % ----- 5) fixation scanpath 그리기 -----
    % 라인(이동 경로)
    plot(fx, fy, '-o', 'LineWidth', 1.5, 'MarkerSize', 6);

    % duration 비례 marker size (선택)
    %   80~800ms 사이 정도 가정하고 scale
    s = 10 + 0.05 * dur;    % 너무 크면 계수 줄이면 됨
    scatter(fx, fy, s, 'filled', 'MarkerFaceAlpha', 0.5);

    % fixation 순서 번호도 표시하고 싶다면:
    for i = 1:numel(fx)
        text(fx(i), fy(i) - 20, num2str(i), ...
            'Color','b', 'FontSize', 8, ...
            'HorizontalAlignment','center');
    end

    % 제목 / 라벨
    title(sprintf('Trial %d sentence + scanpath (%s)', t, coordLabel), ...
        'Interpreter','none');
    xlabel('X (px)');
    ylabel('Y (px)');
    grid on;

    hold off;
end
