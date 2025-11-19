function plotTrialScanpathOnImage(subj, results, t, imgFile, useCorr)
% plotTrialScanpathOnImage
%   PTB 화면 캡처 이미지를 깔고, 그 위에 해당 trial의 fixation scanpath를 그림.
%
%   plotTrialScanpathOnImage(subj, results, t, imgFile, useCorr)
%
%   입력:
%     subj    : Analysis에서 만든 subj struct
%     results : EL_demo.mat에서 로드한 results (단어 정보 등, 여기서는 거의 안 씀)
%     t       : 보고 싶은 trial 번호 (정수)
%     imgFile : 이 trial 화면 캡처 이미지 경로 (png/jpg 등)
%     useCorr : true면 drift correction 된 xCorr/yCorr 사용,
%               false면 raw x/y 사용
%
%   전제:
%     - subj.trial(t).fixIdxDurClean (또는 fixIdx) 존재
%     - subj.event.fix(k).x / .y (그리고 선택적으로 xCorr / yCorr)
%     - subj.trial(t).startTime 존재 (addTrialsFromMsg에서 생성)

    if nargin < 5 || isempty(useCorr)
        useCorr = true;   % 기본은 보정 좌표 사용
    end

    if ~isfile(imgFile)
        error('이미지 파일을 찾을 수 없습니다: %s', imgFile);
    end

    % ----- 1) 이미지 불러오기 -----
    img = imread(imgFile);
    [imgH, imgW, ~] = size(img);

    % ----- 2) 이 trial에서 사용할 fixation 인덱스 선택 -----
    if t > numel(subj.trial)
        error('trial %d는 subj.trial 범위를 벗어납니다.', t);
    end

    tr = subj.trial(t);

    if isfield(tr, 'fixIdxDurClean') && ~isempty(tr.fixIdxDurClean)
        fixIdx = tr.fixIdxDurClean(:)';   % duration 클리닝 반영
    elseif isfield(tr, 'fixIdxValid') && ~isempty(tr.fixIdxValid)
        fixIdx = tr.fixIdxValid(:)';
    else
        fixIdx = tr.fixIdx(:)';
    end

    if isempty(fixIdx)
        warning('Trial %d: fixation이 없습니다.', t);
        return;
    end

    F = subj.event.fix(fixIdx);

    % ----- 3) 좌표와 시간 가져오기 -----
    if useCorr && isfield(F, 'xCorr') && isfield(F, 'yCorr')
        fx = [F.xCorr];
        fy = [F.yCorr];
    else
        fx = [F.x];
        fy = [F.y];
    end

    ft = [F.onset];                         % 절대 시간(ms)
    t0 = tr.startTime;                      % trial 시작 시간(ms)
    tRel = ft - t0;                         % trial 기준 상대 시간(ms)

    % NaN / 화면 밖 fixation 제거(선택)
    valid = ~isnan(fx) & ~isnan(fy) & ...
            fx >= 0 & fx <= imgW & ...
            fy >= 0 & fy <= imgH;
    fx = fx(valid);
    fy = fy(valid);
    tRel = tRel(valid);

    if isempty(fx)
        warning('Trial %d: 화면 안에 있는 fixation이 없습니다.', t);
        return;
    end

    % ----- 4) 그림 그리기 -----
    figure;
    imagesc([0 imgW], [0 imgH], img);
    axis image;
    set(gca, 'YDir','reverse');   % 화면 좌표계와 동일하게 (위가 작은 y)
    hold on;

    % fixation 순서대로 선 + 점
    plot(fx, fy, '-o', 'LineWidth', 1.5, 'MarkerSize', 6);

    % 각 fixation에 번호 및 시간(ms) 레이블 달기
    for k = 1:numel(fx)
        txt = sprintf('%d\n%.0fms', k, tRel(k));
        text(fx(k)+5, fy(k)-5, txt, ...
            'Color','b', 'FontSize',8, 'FontWeight','bold');
    end

    % 제목
    if isfield(results, 'words') && numel(results.words) >= t ...
            && ~isempty(results.words{t})
        sentStr = strjoin(results.words{t}, ' ');
    else
        sentStr = '';
    end

    title(sprintf('Trial %d sentence + scanpath (%s gaze)\n%s', ...
        t, ternary(useCorr,'Corrected','Raw'), sentStr), ...
        'Interpreter','none');

    xlabel('X (px)');
    ylabel('Y (px)');
    grid on;
end

% --- 작은 헬퍼 함수 (MATLAB에는 ternary 연산자가 없으니 직접 정의) ---
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
