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

    isMain = startsWith(string(tr.id), "TRIALID");
    if ~isMain, continue; end

    inTrial = (mTime >= tr.startTime) & (mTime <= tr.endTime);

   ixS = find(inTrial & startsWith(msgText, startToken), 1, "first");
   ixP = find(inTrial & startsWith(msgText, endToken),   1, "first");

    if isempty(ixS) || isempty(ixP)
        subj.trial(t).readStart = NaN;
        subj.trial(t).readEnd   = NaN;
        subj.trial(t).readSampleIdx = [];
        subj.trial(t).readFixIdx    = [];
        subj.trial(t).readSaccIdx   = [];
        subj.trial(t).readBlinkIdx  = [];
        subj.trial(t).readMsgIdx    = [];
        continue;
    end

    readStart = mTime(ixS);
    readEnd   = mTime(ixP);

    if readEnd <= readStart
        subj.trial(t).readStart = readStart;
        subj.trial(t).readEnd   = readEnd;
        subj.trial(t).readSampleIdx = [];
        subj.trial(t).readFixIdx    = [];
        subj.trial(t).readSaccIdx   = [];
        subj.trial(t).readBlinkIdx  = [];
        subj.trial(t).readMsgIdx    = [];
        continue;
    end

    subj.trial(t).readStart = readStart;
    subj.trial(t).readEnd   = readEnd;

    subj.trial(t).readSampleIdx = find(subj.sample.time >= readStart & subj.sample.time <= readEnd);
    subj.trial(t).readFixIdx    = find(fxOn <= readEnd & fxOff >= readStart);
    subj.trial(t).readSaccIdx   = find(scOn <= readEnd & scOff >= readStart);
    subj.trial(t).readBlinkIdx  = find(bkOn <= readEnd & bkOff >= readStart);
    subj.trial(t).readMsgIdx    = find(inTrial & (mTime >= readStart) & (mTime <= readEnd));
end
end
