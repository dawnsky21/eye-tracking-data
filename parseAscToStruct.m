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

    % ---- 파싱 진단용 카운터 ----
    nEFIX    = 0; nEFIX_ok    = 0;
    nEBLINK  = 0; nEBLINK_ok  = 0;
    nSamp    = 0; nSamp_ok    = 0;

    badSampleLines = {};  % 샘플 파싱 실패 라인 저장용

    %% 1) 줄 단위 파싱
    for i = 1:nLines
        tline = strtrim(lines{i});
        if isempty(tline)
            continue;
        end

        % 1) 헤더/설정 라인 → 스킵 (START/END는 빼기)
        if startsWith(tline, '**') || startsWith(tline, 'PUPIL') || ...
           startsWith(tline, 'EVENTS') || startsWith(tline, 'SAMPLES') || ...
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
        if ~isempty(regexp(tline(1), '\d', 'once'))

            % 지금 ASC 포맷: time gx gy pupil 네 개
            C = textscan(tline, '%f %f %f %f');

            if ~any(cellfun(@isempty, C))
                % → 진짜 샘플 포맷(4개 숫자)일 때만 샘플로 인정
                nSamp    = nSamp + 1;      % 실제 샘플 라인 수
                nSamp_ok = nSamp_ok + 1;   % 파싱 성공

                sampleTime(end+1,1) = C{1};
                sampleGx(end+1,1)   = C{2};
                sampleGy(end+1,1)   = C{3};
                samplePup(end+1,1)  = C{4};
            else
                % 샘플 포맷이 아닌 숫자줄 (예: "6.17e-06, -1.54e-05")
                badSampleLines{end+1,1} = tline;
            end

            continue;
        end

        % 그 외 라인은 일단 무시
    end   % for i = 1:nLines

    % ---- 파싱 요약/경고 ----
    fprintf('EFIX parsed:   %d / %d lines\n', nEFIX_ok,   nEFIX);
    fprintf('EBLINK parsed: %d / %d lines\n', nEBLINK_ok, nEBLINK);
    fprintf('SAMPLE parsed: %d / %d lines\n', nSamp_ok,   nSamp);

    if nEFIX_ok < nEFIX
        warning('Some EFIX lines could not be parsed. Check ASC format.');
    end
    if nEBLINK_ok < nEBLINK
        warning('Some EBLINK lines could not be parsed. Check ASC format.');
    end
    if nSamp_ok < nSamp
        warning('Some SAMPLE lines could not be parsed. Check ASC format or sample format (time gx gy pupil).');
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