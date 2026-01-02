function subj = addSampleValidity(subj, screenRect, blinkMargin)
% addSampleValidity  샘플 레벨 validity 플래그 추가
%
%   subj = addSampleValidity(subj, screenRect, blinkMargin)
%
%   screenRect = [xMin yMin xMax yMax], 예: [0 0 1920 1080]
%   blinkMargin = blink 앞뒤로 몇 ms까지 invalid로 볼지 (예: 50)

    if nargin < 2 || isempty(screenRect)
        screenRect = [0 0 1920 1080];
    end
    if nargin < 3 || isempty(blinkMargin)
        blinkMargin = 50;   % ms
    end

    % screenRect 형태 검사
    if any(screenRect(3:4) <= screenRect(1:2))
        error('addSampleValidity: screenRect seems invalid. Expected [xMin yMin xMax yMax].');
    end

    t   = subj.sample.time;
    gx  = subj.sample.gx;
    gy  = subj.sample.gy;
    pup = subj.sample.pupil;

    % 샘플 필드 길이 일관성 체크
    n = numel(t);
    if numel(gx) ~= n || numel(gy) ~= n || numel(pup) ~= n
        error('addSampleValidity: sample fields have inconsistent lengths.');
    end

    isValid     = true(n,1);
    isBlink     = false(n,1);
    isSacc      = false(n,1);
    isOutRange  = false(n,1);
    isBadPupil  = false(n,1);

    % 1) pupil 기준
    badPupil = pup <= 0 | isnan(pup);
    isBadPupil(badPupil) = true;
    isValid(badPupil)    = false;

    % 2) 화면 범위 기준
    xMin = screenRect(1);  yMin = screenRect(2);
    xMax = screenRect(3);  yMax = screenRect(4);

    outRange = gx < xMin | gx > xMax | gy < yMin | gy > yMax | ...
               isnan(gx) | isnan(gy);
    isOutRange(outRange) = true;
    isValid(outRange)    = false;

    % 3) blink 구간 (EBLINK)
    if isfield(subj.event, 'blink') && ~isempty(subj.event.blink)
        for b = 1:numel(subj.event.blink)
            on  = subj.event.blink(b).onset  - blinkMargin;
            off = subj.event.blink(b).offset + blinkMargin;
            idx = t >= on & t <= off;
            isBlink(idx) = true;
            isValid(idx) = false;
        end
    end

    % 4) (옵션) saccade 구간 (ESACC) – 필요하면 invalid로 포함
    if isfield(subj.event, 'sacc') && ~isempty(subj.event.sacc)
        for s = 1:numel(subj.event.sacc)
            on  = subj.event.sacc(s).onset;
            off = subj.event.sacc(s).offset;
            idx = t >= on & t <= off;
            isSacc(idx) = true;
            % 샘플 기반 분석에서는 여기도 invalid로 둘 수 있음:
            % isValid(idx) = false;
        end
    end

    % 5) 결과를 subj.sample에 저장
    subj.sample.isValid    = isValid;
    subj.sample.isBlink    = isBlink;
    subj.sample.isSacc     = isSacc;
    subj.sample.isOutRange = isOutRange;
    subj.sample.isBadPupil = isBadPupil;
end

%% 분석 결과
% whos subj
%  Name      Size               Bytes  Class     Attributes
%
%  subj      1x1             30272928  struct              
%
% screenRect  = [0 0 1920 1080];
% blinkMargin = 50;
%
% subj = addSampleValidity(subj, screenRect, blinkMargin);
% subj.sample
%
% ans = 

%  다음 필드를 포함한 struct:

%          time: [704299×1 double]
%            gx: [704299×1 double]
%            gy: [704299×1 double]
%         pupil: [704299×1 double]
%       isValid: [704299×1 logical]
%       isBlink: [704299×1 logical]
%        isSacc: [704299×1 logical]
%    isOutRange: [704299×1 logical]
%    isBadPupil: [704299×1 logical]

% sum(~subj.sample.isValid)      % invalid 샘플 개수
% sum(subj.sample.isBlink)       % blink 샘플 개수
% sum(subj.sample.isOutRange)    % 화면 밖 샘플 개수
% sum(subj.sample.isBadPupil)    % pupil<=0/NaN 샘플 개수

% ans = 57534
% ans = 57411
% ans = 2084
% ans = 38087

% mean(subj.sample.isValid)      % 전체 샘플 중 valid 비율 (0~1)

% ans = 0.9183

% idx = find(~subj.sample.isValid, 1);   % 첫 번째 invalid 샘플 인덱스
% [subj.sample.time(idx) ...
% subj.sample.gx(idx) ...
% subj.sample.gy(idx) ...
% subj.sample.pupil(idx) ...
% subj.sample.isBlink(idx) ...
% subj.sample.isOutRange(idx) ...
% subj.sample.isBadPupil(idx)]

%ans = 1.0e+07 *

%  1 ~ 7번 열

%    1.1132   -0.0001   -0.0000    0.0001         0    0.0000         0
