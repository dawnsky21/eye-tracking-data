function subj = addReadWindowFromMsg(subj)
% addReadWindowFromMsg: 본시행(TRIALID)에서 읽기 구간(SENTENCE_ONSET~PROMPT_ONSET) 인덱스 생성

startToken = "SENTENCE_ONSET_RIGHTABC";  % 현재 데이터의 실제 SENTENCE_ONSET
endToken   = "PROMPT_ONSET";

msgText = string({subj.event.msg.text}); msgText = msgText(:);
mTime   = [subj.event.msg.time]';

fxOn  = [subj.event.fix.onset]';
fxOff = [subj.event.fix.offset]';
scOn  = [subj.event.sacc.onset]';
scOff = [subj.event.sacc.offset]';
bkOn  = [subj.event.blink.onset]';
bkOff = [subj.event.blink.offset]';

for t = 1:numel(subj.trial)
    tr = subj.trial(t);

    isMain = startsWith(string(tr.id), "TRIALID", "IgnoreCase", true);
    if ~isMain, continue; end

    inTrial = (mTime >= tr.startTime) & (mTime <= tr.endTime);

    ixS = find(inTrial & startsWith(msgText, startToken), 1, "first");
    ixP = find(inTrial & startsWith(msgText, endToken),   1, "first");

    % 기본값
    subj.trial(t).readStart = NaN;
    subj.trial(t).readEnd   = NaN;
    subj.trial(t).readSampleIdx = [];
    subj.trial(t).readFixIdx    = [];
    subj.trial(t).readSaccIdx   = [];
    subj.trial(t).readBlinkIdx  = [];
    subj.trial(t).readMsgIdx    = [];

    if isempty(ixS) || isempty(ixP), continue; end

    readStart = mTime(ixS);
    readEnd   = mTime(ixP);
    if readEnd <= readStart, continue; end

    subj.trial(t).readStart = readStart;
    subj.trial(t).readEnd   = readEnd;

    % sample: global sample time에서 구간 필터 (OK)
    subj.trial(t).readSampleIdx = find(subj.sample.time >= readStart & subj.sample.time <= readEnd);

    % msg: trial msg 중 read window에 해당
    subj.trial(t).readMsgIdx = find(inTrial & (mTime >= readStart) & (mTime <= readEnd));

    % fix/sacc/blink: "해당 trial의 idx"에서만 필터 (정석)
    fi = tr.fixIdx(:);    fi = fi(fi>0 & fi<=numel(fxOn));
    si = tr.saccIdx(:);   si = si(si>0 & si<=numel(scOn));
    bi = tr.blinkIdx(:);  bi = bi(bi>0 & bi<=numel(bkOn));

    subj.trial(t).readFixIdx   = fi(fxOn(fi) <= readEnd & fxOff(fi) >= readStart);
    subj.trial(t).readSaccIdx  = si(scOn(si) <= readEnd & scOff(si) >= readStart);
    subj.trial(t).readBlinkIdx = bi(bkOn(bi) <= readEnd & bkOff(bi) >= readStart);
end
end
