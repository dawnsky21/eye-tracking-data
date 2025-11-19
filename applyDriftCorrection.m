function subj = applyDriftCorrection(subj, driftCfg)
% applyDriftCorrection  trial 단위 drift correction 및 bad drift trial 플래그
%
%   subj = applyDriftCorrection(subj, driftCfg)
%
%   필수/옵션 필드 (driftCfg):
%       .xTrue        : 기준 단어/marker의 이상적인 x 좌표 (스칼라 또는 [nTrial×1])
%       .yTrue        : 기준 단어/marker의 이상적인 y 좌표 (스칼라 또는 [nTrial×1])
%       .refEvent     : 기준 이벤트 이름 (지금 구현은 'SENTENCE_ONSET'만 사용)
%       .minFixDur    : 기준 fixation 최소 duration (ms), 기본 60
%       .maxAbsOffsetX: |x_true - x_obs| 허용 최대 픽셀, 기본 60 (글자 1~2칸 정도 가정)
%       .maxAbsOffsetY: |y_true - y_obs| 허용 최대 픽셀, 기본 40 (line 높이 일부)
%
%   출력:
%       subj.sample.gxCorr / gyCorr      : drift 보정된 sample 좌표
%       subj.event.fix.xCorr / yCorr     : drift 보정된 fixation 좌표
%       subj.event.sacc.sxCorr / syCorr / exCorr / eyCorr : 보정된 saccade 좌표
%       subj.trial(t).driftOffset        : [dx dy] (true - obs)
%       subj.trial(t).driftRefPos        : [x_obs y_obs]
%       subj.trial(t).isBadDrift         : 너무 많이 틀어진 trial 여부

 % 이미 drift correction이 적용된 경우, 한 번 더 적용하지 않기 위한 가드
    if isfield(subj.sample, 'hasDriftCorr') && subj.sample.hasDriftCorr
        warning('applyDriftCorrection: drift already applied to this subj. Skipping.');
        return;
    end
    
% --------- 기본 파라미터 채우기 ---------
    if ~isfield(driftCfg, 'refEvent') || isempty(driftCfg.refEvent)
        driftCfg.refEvent = 'SENTENCE_ONSET';
    end
    if ~isfield(driftCfg, 'minFixDur') || isempty(driftCfg.minFixDur)
        driftCfg.minFixDur = 60;    % ms
    end
    if ~isfield(driftCfg, 'maxAbsOffsetX') || isempty(driftCfg.maxAbsOffsetX)
        driftCfg.maxAbsOffsetX = 60;    % px
    end
    if ~isfield(driftCfg, 'maxAbsOffsetY') || isempty(driftCfg.maxAbsOffsetY)
        driftCfg.maxAbsOffsetY = 40;    % px
    end

    nTrials = numel(subj.trial);

    % xTrue / yTrue: 스칼라면 전체 trial에 복제
    if isscalar(driftCfg.xTrue)
        xTrueAll = repmat(driftCfg.xTrue, nTrials, 1);
    else
        xTrueAll = driftCfg.xTrue(:);
        if numel(xTrueAll) < nTrials
            xTrueAll(end+1:nTrials,1) = xTrueAll(end);  % 마지막 값으로 채우기
        end
    end
    if isscalar(driftCfg.yTrue)
        yTrueAll = repmat(driftCfg.yTrue, nTrials, 1);
    else
        yTrueAll = driftCfg.yTrue(:);
        if numel(yTrueAll) < nTrials
            yTrueAll(end+1:nTrials,1) = yTrueAll(end);
        end
    end

    % --------- 보정 좌표 필드 초기화 (처음 한 번만) ---------
    if ~isfield(subj.sample, 'gxCorr')
        subj.sample.gxCorr = subj.sample.gx;
        subj.sample.gyCorr = subj.sample.gy;
    end

    if ~isempty(subj.event.fix)
        if ~isfield(subj.event.fix, 'xCorr')
            for k = 1:numel(subj.event.fix)
                subj.event.fix(k).xCorr = subj.event.fix(k).x;
                subj.event.fix(k).yCorr = subj.event.fix(k).y;
            end
        end
    end

    if ~isempty(subj.event.sacc)
        if ~isfield(subj.event.sacc, 'sxCorr')
            for k = 1:numel(subj.event.sacc)
                subj.event.sacc(k).sxCorr = subj.event.sacc(k).sx;
                subj.event.sacc(k).syCorr = subj.event.sacc(k).sy;
                subj.event.sacc(k).exCorr = subj.event.sacc(k).ex;
                subj.event.sacc(k).eyCorr = subj.event.sacc(k).ey;
            end
        end
    end

    % --------- trial 루프 ---------
    for t = 1:nTrials
        tr = subj.trial(t);

        % 기준 이벤트 시간 선택 (지금은 SENTENCE_ONSET만)
        switch driftCfg.refEvent
            case 'SENTENCE_ONSET'
                tRef = tr.sentenceOnset;
            otherwise
                % 필요하면 나중에 QUESTION_ONSET 등 추가
                tRef = tr.sentenceOnset;
        end

        % 기준 이벤트가 없으면 skip
        if isnan(tRef)
            subj.trial(t).driftOffset   = [NaN NaN];
            subj.trial(t).driftRefPos   = [NaN NaN];
            subj.trial(t).isBadDrift    = true;  % 기준 없음 → 신뢰도 낮음으로 표시
            continue;
        end

        % 이 trial에서 사용 가능한 fixation 리스트
        fixIdx = tr.fixIdx;
        if isempty(fixIdx)
            subj.trial(t).driftOffset   = [NaN NaN];
            subj.trial(t).driftRefPos   = [NaN NaN];
            subj.trial(t).isBadDrift    = true;
            continue;
        end

        fixAll = subj.event.fix(fixIdx);

        % sentence onset 이후의 fixation 중에서 minFixDur 이상인 첫 fixation 찾기
        onsetVec = [fixAll.onset];
        durVec   = [fixAll.dur];

        candMask = (onsetVec >= tRef) & (durVec >= driftCfg.minFixDur);
        candIdx  = find(candMask, 1, 'first');

        % 없으면 그냥 onset 기준으로 가장 가까운 fixation 사용 (fallback)
        if isempty(candIdx)
            candIdx = find(onsetVec >= tRef, 1, 'first');
        end
        if isempty(candIdx)
            % 정말 없으면 이 trial은 drift 구할 수 없음
            subj.trial(t).driftOffset   = [NaN NaN];
            subj.trial(t).driftRefPos   = [NaN NaN];
            subj.trial(t).isBadDrift    = true;
            continue;
        end

        refFix = fixAll(candIdx);
        xObs = refFix.x;
        yObs = refFix.y;

        xTrue = xTrueAll(t);
        yTrue = yTrueAll(t);

        dx = xTrue - xObs;
        dy = yTrue - yObs;

        % trial 정보에 기록
        subj.trial(t).driftOffset = [dx dy];
        subj.trial(t).driftRefPos = [xObs yObs];

        % 3.3: 너무 많이 틀어진 trial 플래그
        errX = abs(dx);
        errY = abs(dy);
        isBad = (errX > driftCfg.maxAbsOffsetX) | ...
                (errY > driftCfg.maxAbsOffsetY);
        subj.trial(t).isBadDrift = isBad;

        % 3.2: 실제 좌표 보정 (sample / fix / sacc)
        % sample: 이 trial 시간 구간에 포함되는 index
        sampleIdx = tr.sampleIdx;
        subj.sample.gxCorr(sampleIdx) = subj.sample.gxCorr(sampleIdx) + dx;
        subj.sample.gyCorr(sampleIdx) = subj.sample.gyCorr(sampleIdx) + dy;

        % fixation: 이 trial에 속한 fixation index
        for k = 1:numel(fixIdx)
            fi = fixIdx(k);
            subj.event.fix(fi).xCorr = subj.event.fix(fi).xCorr + dx;
            subj.event.fix(fi).yCorr = subj.event.fix(fi).yCorr + dy;
        end

        % saccade: 이 trial에 속한 saccade index
        saccIdx = tr.saccIdx;
        for k = 1:numel(saccIdx)
            si = saccIdx(k);
            subj.event.sacc(si).sxCorr = subj.event.sacc(si).sxCorr + dx;
            subj.event.sacc(si).syCorr = subj.event.sacc(si).syCorr + dy;
            subj.event.sacc(si).exCorr = subj.event.sacc(si).exCorr + dx;
            subj.event.sacc(si).eyCorr = subj.event.sacc(si).eyCorr + dy;
        end
    end
% 이 시점까지 왔다면 drift correction 1회 완료
    subj.sample.hasDriftCorr = true;
end
