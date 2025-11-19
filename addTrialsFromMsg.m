%   - MSG 라인 중 TRIALID, PRACTICE_TRIALID, SENTENCE_ONSET, QUESTION_ONSET,
%     RESPONSE, TRIAL_END 등 반복적으로 사용하는 코드 목록을 정리한다.
%   - 각 trial에 대해:
%       - trialStartTime  = 해당 trial의 TRIALID (또는 PRACTICE_TRIALID) 시간
%       - sentenceOnset   = SENTENCE_ONSET 시간
%       - questionOnset   = QUESTION_ONSET 시간 (있다면)
%       - responseTime    = RESPONSE 시간 (있다면)
%       - trialEndTime    = TRIAL_END 또는 다음 trial 시작 직전 시간
%   - 이 시간 범위를 기준으로 sample과 event를 trial별로 분류:
%       subj(s).trial(t).samples = 해당 trial 시간 구간 내 sample
%       subj(s).trial(t).event   = 해당 trial 시간 구간 내 fix/sacc/blink/msg

function subj = addTrialsFromMsg(subj)
% addTrialsFromMsg  MSG 이벤트를 이용해 trial struct 생성
%
%   subj = addTrialsFromMsg(subj)
%
%   subj.event.msg(k).time / .text 를 읽어서
%   - TRIALID / PRACTICE_TRIALID  → trial 시작
%   - TRIAL_END                   → trial 끝
%   - SENTENCE_ONSET              → 문장 시작
%   - QUESTION_ONSET              → 질문 시작 (있다면)
%   - RESPONSE                    → 응답 시간 (있다면)
%   을 찾아서 subj.trial(t) 를 만든다.

    % ----- MSG 타임라인 정리 -----
    msg   = subj.event.msg;
    mTime = [msg.time];
    mText = {msg.text};

    % helper: 이 MSG가 특정 토큰으로 "시작"하는지 확인
    %  예: 'TRIALID T001_main' 은 isMsg(...,'TRIALID') == true
    isMsg = @(txt, token) strncmp(txt, token, length(token));

    % 1) trial 시작 메시지 인덱스 (TRIALID / PRACTICE_TRIALID)
    isTrialStart = cellfun(@(s) isMsg(s, 'TRIALID') | isMsg(s, 'PRACTICE_TRIALID'), mText);
    trialStartIdx = find(isTrialStart);
    nTrials = numel(trialStartIdx);

    % trial struct 틀
    trial = struct('id', [], 'startTime', [], 'endTime', [], ...
                   'sentenceOnset', [], 'questionOnset', [], 'responseTime', [], ...
                   'sampleIdx', [], 'fixIdx', [], 'saccIdx', [], 'blinkIdx', [], ...
                   'msgIdx', []);

    % ----- 이벤트 시간 미리 벡터로 정리 (성능/가독성용) -----
    sampTime = subj.sample.time;

    if isfield(subj.event, 'fix') && ~isempty(subj.event.fix)
        fixOn  = [subj.event.fix.onset];
        fixOff = [subj.event.fix.offset];
    else
        fixOn  = [];
        fixOff = [];
    end

    if isfield(subj.event, 'sacc') && ~isempty(subj.event.sacc)
        saccOn  = [subj.event.sacc.onset];
        saccOff = [subj.event.sacc.offset];
    else
        saccOn  = [];
        saccOff = [];
    end

    if isfield(subj.event, 'blink') && ~isempty(subj.event.blink)
        blinkOn  = [subj.event.blink.onset];
        blinkOff = [subj.event.blink.offset];
    else
        blinkOn  = [];
        blinkOff = [];
    end

    % ----- trial 루프 -----
    for t = 1:nTrials
        iStart = trialStartIdx(t);
        tStart = mTime(iStart);
        trialIdText = mText{iStart};

        % 1-2) trial end 찾기: 같은 trial 시작–다음 trial 시작 사이 TRIAL_END 중 마지막
        if t < nTrials
            nextStartTime = mTime(trialStartIdx(t+1));
            inThisWindow = mTime > tStart & mTime < nextStartTime;
        else
            % 마지막 trial은 파일 끝까지
            inThisWindow = mTime > tStart;
        end

        isTrialEnd = inThisWindow & cellfun(@(s) isMsg(s, 'TRIAL_END'), mText);
        endIdx = find(isTrialEnd);

        if ~isempty(endIdx)
            tEnd = mTime(endIdx(end));  % 이 trial의 마지막 TRIAL_END 시간
        else
            % TRIAL_END가 없으면 fallback: 다음 trial 시작 직전 or 샘플 끝
            if t < nTrials
                tEnd = nextStartTime - 1;
            else
                tEnd = sampTime(end);
            end
            % 경고: TRIAL_END 누락
            warning('addTrialsFromMsg: trial %d (%s) has no TRIAL_END. Using fallback tEnd=%d.', ...
                    t, trialIdText, tEnd);
        end

        % 2) SENTENCE_ONSET, QUESTION_ONSET, RESPONSE 찾기
        inTrialMsg = mTime >= tStart & mTime <= tEnd;

        sentMask = inTrialMsg & cellfun(@(s) isMsg(s, 'SENTENCE_ONSET'), mText);
        qMask    = inTrialMsg & cellfun(@(s) isMsg(s, 'QUESTION_ONSET'), mText);
        respMask = inTrialMsg & cellfun(@(s) isMsg(s, 'RESPONSE'),       mText);

        sentTimes = mTime(sentMask);
        qTimes    = mTime(qMask);
        respTimes = mTime(respMask);

        if ~isempty(sentTimes)
            sentenceOnset = sentTimes(1);
        else
            sentenceOnset = NaN;
            warning('addTrialsFromMsg: trial %d (%s) has no SENTENCE_ONSET within [%d, %d].', ...
                    t, trialIdText, tStart, tEnd);
        end

        if ~isempty(qTimes)
            questionOnset = qTimes(1);
        else
            questionOnset = NaN;
        end

        if ~isempty(respTimes)
            responseTime = respTimes(1);
        else
            responseTime = NaN;
        end

        % 3) 이 trial 구간에 속한 sample / fix / sacc / blink 인덱스
        %    - sample: time이 [tStart, tEnd] 안에 있는 샘플
        sampleIdx = find(sampTime >= tStart & sampTime <= tEnd);

        %    - fix/sacc/blink: "이벤트 구간"이 trial 구간과 조금이라도 겹치는 것
        %      (onset만이 아니라 onset–offset 전체를 고려)
        if ~isempty(fixOn)
            fixIdx = find(fixOn <= tEnd & fixOff >= tStart);
        else
            fixIdx = [];
        end

        if ~isempty(saccOn)
            saccIdx = find(saccOn <= tEnd & saccOff >= tStart);
        else
            saccIdx = [];
        end

        if ~isempty(blinkOn)
            blinkIdx = find(blinkOn <= tEnd & blinkOff >= tStart);
        else
            blinkIdx = [];
        end

        msgIdx = find(inTrialMsg);

        % 4) struct 채우기
        trial(t).id            = trialIdText;
        trial(t).startTime     = tStart;
        trial(t).endTime       = tEnd;
        trial(t).sentenceOnset = sentenceOnset;
        trial(t).questionOnset = questionOnset;
        trial(t).responseTime  = responseTime;

        trial(t).sampleIdx = sampleIdx;
        trial(t).fixIdx    = fixIdx;
        trial(t).saccIdx   = saccIdx;
        trial(t).blinkIdx  = blinkIdx;
        trial(t).msgIdx    = msgIdx;
    end

    subj.trial = trial;
end
