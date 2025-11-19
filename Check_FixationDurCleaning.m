load('eyeReading_cleanIntermediate.mat', 'subj');  % 예: 중간 결과 불러오기

shortThresh = 80;
longThresh  = 800;
subj = cleanFixationDurations(subj, shortThresh, longThresh);

t = 10;
useFixIdx = subj.trial(t).fixIdxDurClean;
fix = subj.event.fix(useFixIdx);

longIdx = subj.trial(t).longFixIdx;
longFix = subj.event.fix(longIdx);

% 여기서 fix / longFix 내용 inspect, plot 등등
