% 이 함수는 ASC 텍스트를 라인 단위로 읽어서:
%   - 숫자로 시작하는 줄 → sample (time, gx, gy, pupil)
%   - SFIX / EFIX         → fixation event
%   - SSACC / ESACC       → saccade event
%   - SBLINK / EBLINK     → blink event
%   - MSG                 → 메타 이벤트 (trial 시작/끝, 문장/질문/응답 등)
%   - START / END         → recording block 시작/끝
% 로 나누어 subj.sample / subj.event.* struct를 만든다.

function subj = parseAscToStruct(ascFile)

% parseAscToStruct  EyeLink ASC 파일을 읽어서 sample + event struct로 변환
%
%   subj = parseAscToStruct(ascFile)
%
%   - subj.sample.time, gx, gy, pupil
%   - subj.event.msg(time, text)
%   - subj.event.fix(onset, offset, dur, x, y)
%   - subj.event.sacc(onset, offset, dur, sx, sy, ex, ey)
%   - subj.event.blink(onset, offset, dur)

    %% 0) 파일 통째로 읽어서 줄 단위 cell array로 분리
    txt   = fileread(ascFile);                 % 전체 텍스트
    lines = regexp(txt, '\r\n|\n', 'split');   % 줄 분리
    nLines = numel(lines);

    % 저장용 변수
    sampleTime = [];
    sampleGx   = [];
    sampleGy   = [];
    samplePup  = [];

    msgList  = struct('time', {}, 'text', {});
    fixList  = struct('onset', {}, 'offset', {}, 'dur', {}, 'x', {}, 'y', {});

    % saccade / blink / recording
    saccList  = struct('onset', {}, 'offset', {}, 'dur', {}, ...
                       'sx', {}, 'sy', {}, 'ex', {}, 'ey', {});
    blinkList = struct('onset', {}, 'offset', {}, 'dur', {});
    recStartTimes = [];   % START 이벤트 시간들
    recEndTimes   = [];   % END   이벤트 시간들

    % ---- NEW: 설정/메타 ----
    sampleCfgLine = string("");
    sampleCoordType = "UNKNOWN";   % "GAZE" / "HREF" / "UNKNOWN"
    displayRect     = [NaN NaN NaN NaN]; % [L T R B], DISPLAY_COORDS에서 얻음

    % ---- 파싱 진단용 카운터 ----
    nEFIX    = 0; nEFIX_ok    = 0;
    nEBLINK  = 0; nEBLINK_ok  = 0;
    nSamp    = 0; nSamp_ok    = 0;

    badSampleLines = {};  % 샘플 파싱 실패 라인 저장용

    %% 1) 줄 단위 파싱
    for i = 1:nLines
        tline = strtrim(lines{i});
        if isempty(tline), continue; end

        % ===== NEW (2): DISPLAY_COORDS 파싱 (스킵 전에 먼저) =====
        if startsWith(tline, 'DISPLAY_COORDS')
            tok = regexp(tline, '^DISPLAY_COORDS\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', 'tokens', 'once');
            if ~isempty(tok)
                displayRect = cellfun(@str2double, tok);  % [L T R B]
            end
            continue;
        end


        % ===== NEW (1): SAMPLES 라인 저장 + coordType 판정 (스킵 전에 먼저) =====
        if startsWith(tline, 'SAMPLES')
            sampleCfgLine = string(tline);

            % "GAZE"/"HREF" 키워드 기반으로 결정 (설정이 곧 정의)
            if contains(sampleCfgLine, "GAZE", "IgnoreCase", true)
                sampleCoordType = "GAZE";
            elseif contains(sampleCfgLine, "HREF", "IgnoreCase", true)
                sampleCoordType = "HREF";
            else
                sampleCoordType = "UNKNOWN";
            end
            continue;
        end

        % 1) 헤더/설정 라인 → 스킵 (START/END는 빼기)
        if startsWith(tline, '**') || startsWith(tline, 'PUPIL') || ...
           startsWith(tline, 'EVENTS') || ...
           startsWith(tline, 'INPUT')  || ...
           startsWith(tline, 'PRESCALER') || startsWith(tline, 'VPRESCALER')
            continue;
        end

        % 1.5) START / END → recording block 경계
        if startsWith(tline, 'START')
            tok = regexp(tline, '^START\s+(\d+)', 'tokens', 'once');
            if ~isempty(tok)
                recStartTimes(end+1,1) = str2double(tok{1});
            end
            continue;
        end

        if startsWith(tline, 'END')
            tok = regexp(tline, '^END\s+(\d+)', 'tokens', 'once');
            if ~isempty(tok)
                recEndTimes(end+1,1) = str2double(tok{1});
            end
            continue;
        end

        % 2) MSG 라인
        if startsWith(tline, 'MSG')
            tok = regexp(tline, '^MSG\s+(\d+)\s+(.*)$', 'tokens', 'once');
            if ~isempty(tok)
                t    = str2double(tok{1});
                text = strtrim(tok{2});
                msgList(end+1).time = t;
                msgList(end).text   = text;
            end
            continue;
        end

        % 3) Fixation event (EFIX만 사용)
        if startsWith(tline, 'SFIX') || startsWith(tline, 'EFIX')
            if startsWith(tline, 'EFIX')
                nEFIX = nEFIX + 1;  % EFIX 라인 카운트

                % 예: "EFIX L  10956118 10956220 102 957.6 537.6 1000.0"
                C = textscan(tline, 'EFIX %*s %f %f %f %f %f');
                if ~any(cellfun(@isempty, C))
                    nEFIX_ok = nEFIX_ok + 1;  % 파싱 성공

                    tStart = C{1};
                    tEnd   = C{2};
                    dur    = C{3};
                    x      = C{4};
                    y      = C{5};
                    fixList(end+1).onset  = tStart;
                    fixList(end).offset   = tEnd;
                    fixList(end).dur      = dur;
                    fixList(end).x        = x;
                    fixList(end).y        = y;
                end
            end
            continue;
        end

        % 4) Saccade event (ESACC)
        if startsWith(tline, 'SSACC') || startsWith(tline, 'ESACC')
            if startsWith(tline, 'ESACC')
                % 예(일반형):
                % ESACC L  10956090 10956114 24 960.0 540.0 957.6 537.6 0.75
                %            1        2      3   4      5      6      7    8
                C = textscan(tline, 'ESACC %*s %f %f %f %f %f %f %f %f');
                if ~any(cellfun(@isempty, C(1:8)))
                    tStart = C{1};   % onset
                    tEnd   = C{2};   % offset
                    dur    = C{3};   % duration
                    sx     = C{4};   % start x
                    sy     = C{5};   % start y
                    ex     = C{6};   % end x
                    ey     = C{7};   % end y
                    % C{8} 는 peak velocity 등 (원하면 나중에 사용)

                    saccList(end+1).onset = tStart;
                    saccList(end).offset  = tEnd;
                    saccList(end).dur     = dur;
                    saccList(end).sx      = sx;
                    saccList(end).sy      = sy;
                    saccList(end).ex      = ex;
                    saccList(end).ey      = ey;
                end
            end
            continue;
        end

        % 5) BLINK (EBLINK)
        if startsWith(tline, 'SBLINK') || startsWith(tline, 'EBLINK')
            if startsWith(tline, 'EBLINK')
                nEBLINK = nEBLINK + 1;  % EBLINK 라인 카운트

                % 예: "EBLINK L  10956050 10956090 40"
                %                   1        2      3
                C = textscan(tline, 'EBLINK %*s %f %f %f');
                if ~any(cellfun(@isempty, C))
                    nEBLINK_ok = nEBLINK_ok + 1;  % 파싱 성공

                    tStart = C{1};
                    tEnd   = C{2};
                    dur    = C{3};
                    blinkList(end+1).onset  = tStart;
                    blinkList(end).offset   = tEnd;
                    blinkList(end).dur      = dur;
                end
            end
            continue;
        end

        % 6) 숫자로 시작하는 줄 → sample 후보
        if ~isempty(regexp(tline, '^\d', 'once'))
            nSamp = nSamp + 1;

    [ok, t, gx, gy, pup] = parseSample4(tline);

    if ok 
       nSamp_ok = nSamp_ok + 1; 
       sampleTime(end+1,1) = t; 
       sampleGx(end+1,1) = gx; % '.'였으면 NaN 
       sampleGy(end+1,1) = gy; % '.'였으면 NaN 
       samplePup(end+1,1) = pup; % '.'였으면 NaN 
    else 
       badSampleLines{end+1,1} = tline; 
    end 
    continue; 
    end

        % 그 외 라인은 일단 무시
    end   % for i = 1:nLines

     % ---- 파싱 요약 ----
    fprintf('EFIX parsed:   %d / %d lines\n', nEFIX_ok,   nEFIX);
    fprintf('EBLINK parsed: %d / %d lines\n', nEBLINK_ok, nEBLINK);
    fprintf('SAMPLE parsed: %d / %d lines\n', nSamp_ok,   nSamp);

    if any(isnan(displayRect))
        warning('DISPLAY_COORDS not found in ASC. subj.displayRect remains NaN.');
    end

    % ---- NEW: 설정 요약 ----
    fprintf('[ASC CFG] coordType=%s | displayRect=[%g %g %g %g]\n', ...
        char(sampleCoordType), displayRect(1),displayRect(2),displayRect(3),displayRect(4));

    if strlength(sampleCfgLine) > 0
        fprintf('[ASC CFG] %s\n', char(sampleCfgLine));
    end

    if ~isempty(badSampleLines)
        fprintf('Non-sample numeric lines (not parsed as samples): %d\n', numel(badSampleLines));
        disp('--- Bad SAMPLE lines (first few) ---');
        disp(badSampleLines(1:min(5,end)));
    end

    %% 2) subj struct로 정리
    [~, name, ~] = fileparts(ascFile);
    subj.id = name;

    subj.sample.time   = sampleTime;
    subj.sample.gx     = sampleGx;
    subj.sample.gy     = sampleGy;
    subj.sample.pupil  = samplePup;

    % ===== NEW =====
    subj.sample.cfgLine   = sampleCfgLine;
    subj.sample.coordType = sampleCoordType;
    subj.displayRect      = displayRect;
    % ===============

    subj.event.msg   = msgList;
    subj.event.fix   = fixList;
    subj.event.sacc  = saccList;
    subj.event.blink = blinkList;

    % recording block 정보 (스타트/엔드 쌍으로 만들기)
    nRec = min(numel(recStartTimes), numel(recEndTimes));
    recList = struct('startTime', {}, 'endTime', {});
    for r = 1:nRec
        recList(r).startTime = recStartTimes(r);
        recList(r).endTime   = recEndTimes(r);
    end
    subj.event.rec = recList;

    function [ok, t, gx, gy, pup] = parseSample4(tline)
    % time gx gy pupil 4개만 안전 파싱
    % '.' 또는 비수치 토큰 -> NaN 허용

    tok = regexp(strtrim(tline), '\s+', 'split');
    if numel(tok) < 4
        ok = false; t=NaN; gx=NaN; gy=NaN; pup=NaN; return;
    end

    tok4 = tok(1:4);

    % (A) time 토큰은 정수(ms)만 허용
    % "12345" OK, "1.2", "1e-5", "0.0001," 전부 탈락
    if isempty(regexp(tok4{1}, '^\d+$', 'once'))
        ok = false; t=NaN; gx=NaN; gy=NaN; pup=NaN; return;
    end

    % (B) gx/gy/pupil만 콤마 제거 (time은 건드리지 않음)
    for j = 2:4
        tok4{j} = regexprep(tok4{j}, ',', '');
    end

    v = nan(1,4);
    for j = 1:4
        if strcmp(tok4{j}, '.') || strcmpi(tok4{j}, 'NA') || strcmpi(tok4{j}, 'NaN')
            v(j) = NaN;
        else
            v(j) = str2double(tok4{j}); % 숫자면 숫자, 아니면 NaN
        end
    end

    % time은 숫자여야 “샘플”로 인정
    ok = ~isnan(v(1));
    t   = v(1);
    gx  = v(2);
    gy  = v(3);
    pup = v(4);
    end
end

%% 분석 결과
% >> subj = parseAscToStruct(ascFile);

% EFIX parsed:   4142 / 4142 lines
% EBLINK parsed: 390 / 390 lines
% SAMPLE parsed: 704299 / 704308 lines
% 경고: Some SAMPLE lines could not be parsed. Check ASC
% format or sample format (time gx gy pupil). 
% > parseAscToStruct (205번 라인) 
% --- Bad SAMPLE lines (first few) ---
%     {'6.1723e-06, -1.5413e-05'}
%     {'1.0666e-05, -3.4343e-05'}
%     {'5.2101e-06, -1.3849e-05'}
%     {'3.064e-05,  4.7819e-06' }
%     {'4.8805e-05, -1.3954e-05'}

% EyeLink sample 라인(time gx gy pupil)이 아니고,
% 6.1723e-06, -1.5413e-05 처럼
% 숫자 두 개 + 콤마 구조의 “계수/파라미터”처럼 보이는 라인

% 실제 시선 샘플이 아니라, 파일 안에 섞여 있는 ‘이상한 수치 줄’ 9개를 우리가 잘 잡아낸 것이고,
% EyeLink가 때때로 내뱉는 계산 계수/설정 잔여물일 가능성이 크다.
% → 분석에는 안 쓰는 정보라 무시해도 된다.