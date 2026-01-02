% revised on Oct.20.2025
% this script was tested only in the dummy mode
% see if EL functions (drift correction, recalibration, saving files, etc) work well before using any parts of this script

function visualsearch_EL_demo

    % ===== SAFETY RESET (코드 맨 위, PTB 초기화 전에) =====
    ListenChar(0);
    try
        KbQueueRelease(-1);     % 모든 키보드 장치의 queue 강제 해제
    catch
    end

    try
        Eyelink('Shutdown');    % 이전 run에서 Eyelink가 비정상 종료됐을 때 대비
    catch
    end

    Screen('CloseAll');         % 혹시 남아 있는 윈도우 정리
    ShowCursor;
    Priority(0);
    % =====================================================


 % --- 프리앰블: 시작 전 동기화 셀프체크 ---
    sync = ptb_sync_safety_preamble();

    % (선택) 통과 여부 정책
    if ~sync.pass
        % 여기서 중단할지, 경고만 하고 진행할지 결정
        warning('[SYNC] Preflight failed → SkipSyncTests=1로 진행(저정확). 환경 점검 권장.');
        % 원하는 경우: return;  % ← 아예 중단
    end

% === 키보드 베이스라인 초기화 (루프 밖, 파일 상단) ===
% === 키 이름/장치 준비 ===
KbName('UnifyKeyNames');

% 1) KeyNames 불러와서 '문자형 행벡터'만 안전 추출
knRaw = KbName('KeyNames');
if iscell(knRaw)
    keynames = knRaw( cellfun(@(x) ischar(x) && (isrow(x) || isempty(x)), knRaw) );
    keynames = keynames(~cellfun(@isempty, keynames));
elseif isstring(knRaw)
    keynames = cellstr(knRaw(:));
    keynames = keynames(~cellfun(@isempty, keynames));
elseif ischar(knRaw)   % 드물게 char matrix로 오는 환경 대비
    keynames = cellstr(knRaw);
    keynames = keynames(~cellfun(@isempty, keynames));
else
    keynames = {};
end

% 2) 모두 소문자로 통일 (후속 비교에서 혼선 방지)
keynames = lower(strtrim(keynames));

% 3) 소문자 후보군과의 안전 교집합 헬퍼
safeIntersect = @(cands) unique(intersect(lower(cands(:).'), keynames, 'stable'));

% 4) 후보 키 집합 (NumLock ON/OFF, 화살표 맵핑까지 모두 포함)
cand_RESP1 = {'KP_1','NUMPAD1','1','1!','END','KP_END'};                % "1"
cand_RESP2 = {'KP_2','NUMPAD2','2','2@','DOWNARROW','KP_DOWN'};         % "2"
cand_ENTER = {'KP_Enter','NumpadEnter','return','enter'};               % 키패드/일반 엔터
cand_SPACE = {'space'};
cand_ADMIN = {'escape','r'};                                            % ESC / 재캘

% 5) 실제 존재하는 키만 선택
RESP1_KEYS = safeIntersect(cand_RESP1);
RESP2_KEYS = safeIntersect(cand_RESP2);
PAD_ENTER  = safeIntersect(cand_ENTER);
KBD_SPACE  = safeIntersect(cand_SPACE);
KBD_ADMIN  = safeIntersect(cand_ADMIN);

% 6) 폴백 (환경에 따라 일부 이름이 없을 수 있음)
if isempty(RESP1_KEYS), RESP1_KEYS = {'1','1!'}; end
if isempty(RESP2_KEYS), RESP2_KEYS = {'2','2@'}; end
if isempty(PAD_ENTER),  PAD_ENTER  = {'return'}; end

% === (선택) 대표 키코드 1개씩 뽑아 KbCheck 분기에 사용하고 싶다면 === 
keyA   = KbName(RESP1_KEYS{1});   % "1" 계열(자연스러움)
keyL   = KbName(RESP2_KEYS{1});   % "2" 계열(어색함)
keyEsc = KbName('escape');
recal  = KbName('r');

% === (권장) 여러 후보 전부 허용하려면 집합 코드도 만들어 두기 ===
codes_RESP1 = KbName(cellstr(RESP1_KEYS));   % 예: any(kc(codes_RESP1))
codes_RESP2 = KbName(cellstr(RESP2_KEYS));
codes_ADMIN = KbName(cellstr(KBD_ADMIN));

% === 장치 선택: 키패드/키보드 분리 ===
padIdx = pick_keypad_device(8);   % 8초 동안 외장 키패드 아무 키를 눌러 선택 
% ← 로컬 함수로 추가 필요(아래 설명)
kbdIdx = -1;                     % 기본 키보드(전체)  → ESC/R/SPACE/1/2 받기

% 두 개의 입력 핸들(Queue 사용)
kbPad = init_keyboard(struct('useKbQueueCheck',true,'restrict',true,'deviceIndex',padIdx));
kbKb  = init_keyboard(struct('useKbQueueCheck',true,'restrict',true,'deviceIndex',kbdIdx));

kbPad.clear();  kbKb.clear();
kbPad.debounce(); kbKb.debounce();

kb = kbKb;   % 이후 기존의 kb.only / kb.wait 호출은 '키보드' 큐를 사용

try
s = RandStream('mt19937ar','Seed','Shuffle'); %난수 발생기(PRNG) 스트림을 만드는 코드, 매번 다른 난수열이 나옴
RandStream.setGlobalStream(s); %방금 만든 스트림 s를 전역 난수원으로 지정, 기본 난수 함수(rand, randn)들이 이 스트림을 사용하도록 함
KbName('UnifyKeyNames'); %운영체제별 키 이름 차이 통일

% === Save to .../VisualSearch/results/<yyyymmdd_HHMMSS> ===
thisFile    = mfilename('fullpath');           % .../VisualSearch/scripts/visualsearch_EL_demo.m
projectRoot = fileparts(fileparts(thisFile));  % .../VisualSearch
stamp       = datestr(now,'yyyymmdd_HHMMSS');  % 예: 20250904_143210
dirname = fullfile(projectRoot, 'results', stamp);
if ~isfolder(dirname), mkdir(dirname); end
fprintf('[SAVE DIR] %s\n', dirname); 
 
bgColor   = [128 128 128]/255; %배경색=중간 회색
 
% timing
fixdur = 1; % required fixation duration to display/remove the sentence stimulus (in sec) %시선 고정 시간 1초로 설정

%% Eyelink initial settings
dummymode = 0; % 0=real, 1=dummy
isReal    = false;

% 1) 모드 결정
if dummymode == 0
    ok = EyelinkInit(0);  % 실장비 연결 시도

    if ~ok
        % 아직 EDF 안 만들었으므로 filetransfer는 의미 없음
        % 그냥 안내 문구 + 에러로 종료
        Screen('Preference','TextEncodingLocale','UTF-8');
        try
            % 만약 아직 w 안 열었다면 OpenWindow 전에 에러가 나니,
            % 여기서는 간단히 command window만 써도 됨
            fprintf('[MODE] EyeLinkInit 실패. 트래커/케이블/네트워크 상태를 확인하세요.\n');
        end
        error('[MODE] EyeLinkInit failed. EyeLink not connected.');
    end

    if Eyelink('IsConnected') ~= 1
        fprintf('[MODE] EyeLink not connected after init.\n');
        error('[MODE] EyeLink not connected. 실험을 중단합니다.');
    end

    isReal = true;
else
    EyelinkInit(1);
end

% 2) EDF 파일명 결정
if isReal
    prompt  = {'Enter tracker EDF file name (1–8 A-Z/0-9)'};
    answer  = inputdlg(prompt, 'Create EDF file', 1, {'DEMO'});
    if isempty(answer) || isempty(answer{1})
        edfFile = 'DEMO';
    else
        edfFile = upper(regexprep(answer{1}, '[^A-Z0-9]', ''));
        if isempty(edfFile), edfFile = 'DEMO'; end
        if numel(edfFile) > 8, edfFile = edfFile(1:8); end
    end
else
    edfFile = 'DUMMY';
end
if isReal
    fprintf('[MODE] REAL | EDF=%s\n', edfFile);
else
    fprintf('[MODE] DUMMY | EDF=%s\n', edfFile);
end

% 3) EDF 오픈은 "실모드에서만"
if isReal
    openerror = Eyelink('Openfile', edfFile);
    if openerror
        fprintf('Cannot create EDF file: %s\n', edfFile);
        cleanup; return;
    end

    % 이 줄을 여기(실모드 블록 안)로 이동
    preambleText = sprintf('Visual search demo: This session was initiated at approximately: %s', ...
                           datetime('now','TimeZone','local'));
    Eyelink('Command','add_file_preamble_text "%s"', preambleText);
end   

if dummymode == 0 % 실모드(장비 연결)
    % Get EyeLink tracker and software version
    [v, vs]=Eyelink('GetTrackerVersion'); % [v, vs]: 트래커 버전 정보 읽기, v : 숫자 버전 (예: 5), vs: 문자열 설명 (예: 'EYELINK 1000 Plus 5.15') 
    fprintf('Running experiment on a %s version %d\n', vs, v); % 어떤 모델/버전에서 돌고 있는지 콘솔에 출력
    vsn = regexp(vs,'\d','match'); % vs 문자열에서 숫자 한 자리씩 뽑아 셀 배열로 저장
else 
    % dummy mode won't return tracker version so placeholders are needed.
    v = 0; % 더미 모드: 실제 값을 못 읽으니 플레이스홀더로 v=0, vsn{1}='0' 넣어 둠
    vsn{1} = '0';
end

% make sure that we get gaze data from the Eyelink(Eyelink에서 시선 데이터가 제대로 수집되도록 확인)
% Select which events are saved in the EDF file. Here: include everything(어떤 이벤트를 EDF 파일에 저장할지 선택, 오프라인 분석용, 여기서는 전부 포함)
Eyelink('Command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT'); % LEFT, RIGHT: 왼쪽/오른쪽 눈의 이벤트 모두 포함; FIXATION: 고정 이벤트(시작/끝); SACCADE: 도약 이벤트; BLINK: 눈깜빡임 이벤트; MESSAGE: Eyelink('Message',...)로 찍는 타임스탬프·마커; BUTTON: 버튼박스 입력; INPUT: 아날로그/디지털 입력(있을 때)
% Select which events are available in real time for gaze-contingent experiments. Here: include everything(실시간으로 받을 이벤트 선택, 응시 연동 자극 제어용, 여기서도 전부 포함)
Eyelink('Command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');

if v == 3 && str2double(vsn{1}) == 5 % As of September 2023, v = 3, vs = EYELINK CL 5.50
    % Check tracker version and include 'HTARGET' to save head target sticker data for supported eye trackers
    Eyelink('Command', 'file_sample_data  = LEFT,RIGHT,GAZE,GAZERES,PUPIL,HREF,AREA,HTARGET,BUTTON,STATUS,INPUT'); % file_sample_data : EDF(오프라인) 에 저장할 연속 샘플 항목; 주요 샘플 항목: LEFT, RIGHT : 좌/우안 데이터, GAZE : 시선 좌표, GAZERES : 시선 좌표 해상도 메타, PUPIL : 동공 크기, HREF : 머리 기준 좌표계(HREF) 데이터, AREA : 동공 면적, HTARGET : 헤드 타깃(스티커) 관련 데이터(해당 트래커에서 지원 시), BUTTON : 버튼 입력, STATUS : 유효성/상태 비트, INPUT : 외부 입력(디지털/아날로그)
    Eyelink('Command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,HTARGET,STATUS,INPUT'); % link_sample_data : 실시간 링크(온라인 제어)로 받을 연속 샘플 항목(대역폭 절약 위해 보통 더 적게)
else
    Eyelink('Command', 'file_sample_data  = LEFT,RIGHT,GAZE,GAZERES,PUPIL,HREF,AREA,BUTTON,STATUS,INPUT');
    Eyelink('Command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT');
end

% === 실험 참가자 정보 수집 ===
SONAID = input('소나 아이디를 입력해주세요: ', 's');

validGender = false;
while ~validGender
    gender = input('성별을 입력해주세요 (M/F): ', 's');
    if strcmpi(gender, 'M') || strcmpi(gender, 'F')
        validGender = true;
    else
        disp('Invalid input. Please enter M for Male or F for Female.');
    end
end

age = input('만나이를 입력해주세요: ');

validHand = false;
while ~validHand
    hand = input('주로 사용하는 손은 무엇인가요? (R/L): ', 's');
    if strcmpi(hand, 'R') || strcmpi(hand, 'L')
        validHand = true;
    else
        disp('Invalid input. Please enter R for Right hand or L for Left hand.');
    end
end

result.info = struct('studentID', SONAID, ...
                     'gender',     gender, ...
                     'age',        age, ...
                     'hand',       hand);

% 결과 저장 파일명 설정
saveFileName = ['results_Exp_' SONAID '.mat'];

%% Open Screen
whichscreen = max(Screen('Screens'));

% ==== (선택) 개발용 반투명 디버그 창 ====
% 실제 실험 땐 주석
%PsychDebugWindowConfiguration;

% ==== OpenWindow에 넘길 디스플레이 설정 (cm 단위 보장) ====
dp.screenNum  = max(Screen('Screens'));
dp.bkColor    = bgColor;          % [128 128 128] 유지 가능

% 필수: cm 단위 보장 (직접 지정 or 자동 추정) 
if ~isfield(dp, 'width') || isempty(dp.width) || dp.width < 10
    % 하드웨어에서 물리 폭(mm) 읽어와 cm로 변환 (가능하면 이 경로 사용)
    try
        [wmm, ~] = Screen('DisplaySize', dp.screenNum);  % mm
        if ~isempty(wmm) && wmm > 0
            dp.width = wmm / 10;    % cm
        else
            dp.width = 60;          % fallback: 60 cm (모니터 물리폭 직접 기입 권장)
        end
    catch
        dp.width = 60;              % fallback
    end
end
if ~isfield(dp, 'dist') || isempty(dp.dist) || dp.dist < 10
    dp.dist = 98;                   % cm (실험 장비 실제 거리로 교체 권장)
end

% (선택) 스테레오 쓰면 지정
if ~isfield(dp,'stereoMode'), dp.stereoMode = 0; end

% (1) 텍스트/렌더링 선설정 — 반드시 OpenWindow 이전
Screen('Preference','TextRenderer', 1);
Screen('Preference','TextEncodingLocale','UTF-8');  
Screen('Preference','DefaultFontSize', 40);

% --- 폰트: D2Coding 고정 (+안전 폴백) ---
desiredFont = 'D2Coding';
try
    fList = listfonts;
    if any(strcmpi(fList, desiredFont))
        pickedFont = desiredFont;
    else
        % 폴백 후보 (원하는 순서대로)
        fallbacks = {'NanumGothicCoding','Malgun Gothic','Arial'};
        hit = find(ismember(lower(fallbacks), lower(fList)), 1, 'first');
        if ~isempty(hit)
            pickedFont = fallbacks{hit};
        else
            pickedFont = 'Arial';  % 최종 폴백
        end
        warning('"%s" not found. Falling back to "%s".', desiredFont, pickedFont);
    end
catch
    pickedFont = desiredFont;  % listfonts 실패 시에도 D2Coding 시도
end
pickedRenderer = 1;  % FreeType

% ==== 창 열기 (여기서 ifi/프레임/ppd/센터 다 채워짐) ====
dp.width = 60;   % cm
dp.dist  = 98;   % cm
dp = OpenWindow(dp);

% --- 여기 추가: 현재 창 크기로 좌표계 동기화 ---
rect = dp.wRect;                 % [0 x0 y0 x1 y1] 형식일 경우는 아래처럼 width/height만 사용
scrW = rect(3); scrH = rect(4);
dp.resolution = [scrW scrH];     % 이후 grid/좌표 계산도 같은 기준 쓰게 동기화

if dummymode == 0                        % 실모드에서만 호스트에 통보
    Eyelink('Command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, scrW-1, scrH-1);
    Eyelink('Message','DISPLAY_COORDS %ld %ld %ld %ld',        0, 0, scrW-1, scrH-1);
end

% 핸들 꺼내 쓰기
w = dp.wPtr;
fprintf('Window pointer w = %d\n', w);
if w <= 0 || Screen('WindowKind', w) <= 0
    error('Invalid window handle detected. Psychtoolbox window was not opened correctly.');
end
rect = dp.wRect;
cx   = dp.cx;
cy   = dp.cy;

% ==== 텍스트 상태를 D2Coding으로 보장 ====
Screen('Preference','TextRenderer', pickedRenderer);
Screen('TextFont',  w, pickedFont);   % ← D2Coding(설치됨)
Screen('TextSize',  w, 40);
Screen('TextStyle', w, 0);

% ==== picked* 변수 미리 지정 (나중 코드 호환용) ====
pickedRenderer = 1;

% (2) 개발용 타이밍 스킵/디버그 창
%PsychDebugWindowConfiguration;      % 반투명 디버그 창

DEBUG_FONT = false;  % 개발 중에만 true

% 디버그가 필요할 때만 주석 해제해서 사용
% disp(which('Screen','-all'));
% disp(which('DrawFormattedText','-all'));

% % to linearize the monitorf
% load('BenQ_20220808.mat'); 정확한 광도 제어가 필요한 실험에서만 주석 해제해서 사용. 파일(BenQ_20220808.mat)은 실측한 모니터별로 달라.
% gammaTable = repmat(luminanceGamma',1,3);모니터가 입력값 0~255에 대해 실제 밝기가 선형이 아니기 때문에, 감마 테이블로 선형에 가깝게 교정, 측정해 둔 감마 곡선(luminanceGamma)을 불러와 RGB 3채널로 복제 → LoadNormalizedGammaTable로 적용.
% oldGamma = Screen('LoadNormalizedGammaTable', w, gammaTable);

% === ABC 위치(가로) 픽셀 환산 ===
ppcmX          = dp.resolution(1) / dp.width;  % px per cm (가로)
ABC_OFFSET_CM  = 12;                           
ABC_OFFSET_PX = round(ABC_OFFSET_CM * ppcmX); % 픽셀 오프셋

% === 텍스트 왼쪽 기준 X (화면 왼쪽에서 8cm 지점) ===
TEXT_LEFT_MARGIN_CM = 12;              % 원하면 5~12cm 사이로 조정
TEXT_LEFT_X         = round(TEXT_LEFT_MARGIN_CM * ppcmX);

%% === (삽입) 엑셀 읽기 + 4블록(56개) 순서 만들기 ===
Screen('Preference','TextEncodingLocale','UTF-8');
fprintf('WindowKind at TextSize: %d\n', Screen('WindowKind', w));

% 윈도우가 살아있는지, 폰트/렌더러/사이즈가 설정돼있는지 한번에 보장
[w, rect, cx, cy] = ensureWindow(w, whichscreen, bgColor, 40, pickedFont, pickedRenderer);

% (선택) 개발 중 점검
assert(Screen('WindowKind', w) > 0, 'Window handle lost before Excel/Text section.');

xlsxPath = fullfile(pwd, 'Experimental stimulus_실험용_수정본7.xlsx');

if ~isfile(xlsxPath)
    error('Stimulus excel not found: %s', xlsxPath);
end
reqSheets = {'Main','Practice'};
[~, sheets] = xlsfinfo(xlsxPath);
if any(~ismember(reqSheets, sheets))
    error('Excel missing required sheets. Found: %s', strjoin(sheets, ', '));
end
% --- 검사 끝 ---

% 실제 읽기
T = readtable(xlsxPath, 'Sheet','Main', 'TextType','string');

% --- target_idx / target_word 존재 여부 검사 ---
needVars = ["sentence","freq","valence","is_catch","target_word","target_idx"];
if any(~ismember(needVars, string(T.Properties.VariableNames)))
    error('엑셀 Main 시트에 필요한 변수들이 없습니다. 필요 변수: %s', strjoin(needVars, ', '));
end

% (선택) 나중에 쓰기 편하도록 이름 맞춰두기
T.targetWord   = T.target_word;
if isnumeric(T.target_idx)
    T.targetWordIdx = double(T.target_idx);
else
    T.targetWordIdx = str2double(string(T.target_idx));
end

% --- (중요) non-catch인데 target_idx 비어 있지 않은지 검사 ---
badTarget = find(~T.is_catch & isnan(T.targetWordIdx));
if ~isempty(badTarget)
    error('Non-catch trial 중 target_idx가 NaN인 행이 있습니다. row index: %s', mat2str(badTarget));
end

% --- 데이터 최소량·셀 커버리지 검사 ---
if height(T)==0
    error('Main table empty: 실험용 메인 시트(Main)가 비어 있습니다.');
end

cells = ["HF","LF"];   % 빈도 조건
vals  = ["P","N","U"]; % 정서 조건 (Positive/Negative/Neutral)

for f = cells
    for v = vals
        n = nnz(~T.is_catch & T.freq==f & T.valence==v);
        if n == 0
            error('No items for cell %s_%s (빈 셀: %s-%s).', f, v, f, v);
        end
    end
end

% === 본실험 테이블 T 정규화 ===
T.sentence = string(T.sentence);
T.freq    = upper(strtrim(string(T.freq)));
T.valence = upper(strtrim(string(T.valence)));
if islogical(T.is_catch)
    % ok
elseif isnumeric(T.is_catch)
    T.is_catch = T.is_catch ~= 0;
else
    xs = lower(strtrim(string(T.is_catch)));
    T.is_catch = xs=="1" | xs=="y" | xs=="yes" | xs=="true" | xs=="t";
end
T.is_catch = logical(T.is_catch);

% === [진단 출력: 여기에 추가] ===
fprintf('[CHECK] rows=%d | noncatch=%d | catch=%d\n', height(T), sum(~T.is_catch), sum(T.is_catch)); 

uf = unique(strtrim(string(T.freq)));
uv = unique(strtrim(string(T.valence)));
disp(table(uf, 'VariableNames',{'unique_freq'}));
disp(table(uv, 'VariableNames',{'unique_valence'}));

if any(~T.is_catch)
    G = groupsummary(T(~T.is_catch,:), {'freq','valence'});
    disp(G(:, {'freq','valence','GroupCount'}));
end
% === [진단 끝] ===

% === Practice sheet 읽기 ===
T_prac = readtable(xlsxPath, 'Sheet','Practice', 'TextType','string');

% 컬럼 정리(본 실험과 동일 형식 보장)
T_prac.sentence = string(T_prac.sentence);
T_prac.freq     = string(T_prac.freq);
T_prac.valence  = string(T_prac.valence);
% is_catch robust 처리
if islogical(T_prac.is_catch)
    % ok
elseif isnumeric(T_prac.is_catch)
    T_prac.is_catch = T_prac.is_catch ~= 0;
else
    xs = lower(strtrim(string(T_prac.is_catch)));
    T_prac.is_catch = xs=="1" | xs=="y" | xs=="yes" | xs=="true" | xs=="t";
end
T_prac.is_catch = logical(T_prac.is_catch);

% 문장 정제
T_prac.sentence = arrayfun(@cleanSentence, T_prac.sentence);

if height(T_prac)==0
    warning('Practice sheet empty: 연습 시트가 비어 있어 연습 블록을 건너뜁니다.');
    nPractice = 0;
end

% === Practice 역균등화(10개 전부 사용) ===
% 아이디어: catch가 인접하지 않게, 동일 freq/valence 3연속 방지 (본 로직 축약판)
idxAll = (1:height(T_prac))';
ok=false; tries=0;
freqOfP = @(idx) string(T_prac.freq(idx));
valOfP  = @(idx) string(T_prac.valence(idx));
isCatchP= @(idx) T_prac.is_catch(idx);

while ~ok && tries<2000
    tries=tries+1;
    seq = idxAll(randperm(numel(idxAll)));  % 전부 섞기
    ok = true;
    % [연습] 맨 앞 catch 금지만 유지
if isCatchP(seq(1)), ok=false; end

    % 동일 freq/valence 3연속 금지
    if ok
        for i_=1:numel(seq)-2
            f3=freqOfP(seq(i_:i_+2)); v3=valOfP(seq(i_:i_+2));
            if all(f3=="HF") || all(f3=="LF") || all(v3=="P") || all(v3=="N") || all(v3=="U")
                ok=false; break;
            end
        end
    end
end
PracOrder = seq(:);
Practice  = T_prac(PracOrder,:);
nPractice = height(Practice);


T.sentence = string(T.sentence);
T.freq     = string(T.freq);
T.valence  = string(T.valence);

% ← 여기! 총량 점검
fprintf('non-catch total=%d | catch total=%d\n', sum(~T.is_catch), sum(T.is_catch));

% 디버그 프린트(처음 1개 미리보기)
ix = find(~T.is_catch & strlength(T.sentence)>1, 1, 'first');
if ~isempty(ix)
    fprintf('[XL] example sentence: "%s"\n', T.sentence(ix));
else
    fprintf('[XL] non-catch sentence preview not found.\n');
end

% 6 셀(HF/LF × P/N/U)과 캐치 인덱스
cells = ["HF_P","HF_N","HF_U","LF_P","LF_N","LF_U"];
cellMask = containers.Map();
for k = 1:numel(cells)
    key = char(cells(k));           % ← char로 통일
    f = extractBefore(cells(k),"_");
    v = extractAfter(cells(k),"_");
    cellMask(key) = find(~T.is_catch & T.freq==f & T.valence==v);
end
idxCatch = find(T.is_catch);

% 디버그 출력
for k = 1:numel(cells)
    key = char(cells(k));           % ← 조회도 char
    fprintf('%s: %d\n', key, numel(cellMask(key)));
end
fprintf('CATCH: %d\n', numel(idxCatch));
fprintf('NON-CATCH TOTAL: %d\n', sum(cellfun(@numel, values(cellMask))));

% --- 4블록 설계 ---
nBlocks   = 4;
blockSize = 56;

% 기본 패턴(31개 → 8,8,8,7)을 보유 개수에 맞춰 일반화
cellCounts = zeros(1,numel(cells));
for k = 1:numel(cells)
    key = char(cells(k));
    cellCounts(k) = numel(cellMask(key));
end

% --- 이후 섞기 제약(동일 freq/valence 3연속 금지)은 기존 코드 그대로 ---

% ==== 역균등화 기반 블록 구성 (56 trials × 4 blocks) ====

rng('shuffle');

% 1) 각 셀(HF/LF × P/N/U)별 인덱스 풀을 미리 무작위화
pools   = cell(1, numel(cells));
offs    = zeros(1, numel(cells));   % 셀별 오프셋
for c = 1:numel(cells)
    key = char(cells(c));
    v = cellMask(key);
    pools{c} = v(randperm(numel(v)));
end

% 2) 캐치도 전체를 미리 무작위화
idxCatch = idxCatch(randperm(numel(idxCatch)));

% 3) 블록별로 '역균등화' 순서 만들기
orderIdx = cell(1, nBlocks);

% 캐치를 블록에 최대한 균등 분배
totalCatch = numel(idxCatch);
baseC = floor(totalCatch / nBlocks);
remC  = totalCatch - baseC*nBlocks;
catchPerBlock = baseC*ones(1,nBlocks);
if remC > 0
    addHere = randperm(nBlocks, remC);
    catchPerBlock(addHere) = catchPerBlock(addHere) + 1;
end

% 각 블록의 비-캐치 필요 개수
needNonPerBlock = blockSize - catchPerBlock;

% === 여기서 총량 점검을 실행 ===
needNonTotal = sum(needNonPerBlock);
totalNon = sum(cellfun(@numel, pools));

if totalNon < needNonTotal
    error('비-캐치 총량 부족: 필요=%d, 보유=%d', needNonTotal, totalNon);
end

% --- is_catch robust: 숫자/문자 모두 처리 (단 한 번만) ---
for b = 1:nBlocks
    needNon = needNonPerBlock(b);
    seq = zeros(1,0);

    % 남은 비-캐치 수 확인
    remainNonGlobal = sum(cellfun(@(x,o) numel(x)-o, pools, num2cell(offs)));
    if remainNonGlobal < needNon
        error('비-캐치 문장 수가 부족합니다: 필요=%d, 남은=%d', needNon, remainNonGlobal);
    end

    % 라운드로빈으로 비-캐치 채우기
    while numel(seq) < needNon
        progressed = false;
        for c = 1:numel(cells)
            if numel(seq) >= needNon, break; end
            if offs(c) < numel(pools{c})
                offs(c) = offs(c) + 1;
                seq(end+1) = pools{c}(offs(c)); %#ok<AGROW>
                progressed = true;
            end
        end
        if ~progressed, break; end
    end

    % 캐치 삽입(첫/끝 상관없음)
    kcatch = catchPerBlock(b);
    if kcatch > 0
        if numel(idxCatch) < kcatch
            error('캐치 문장 수가 부족합니다: 필요=%d, 남은=%d', kcatch, numel(idxCatch));
        end
        insertPos = round(linspace(1, numel(seq)+1, kcatch));
        for kk = 1:kcatch
            pos = insertPos(kk);
            seq = [seq(1:pos-1), idxCatch(1), seq(pos:end)]; %#ok<AGROW>
            idxCatch(1) = [];
            insertPos = insertPos + 1;
        end
    end

    %% === [PATCH] 동일 freq/valence 3연속 금지(사후 보정; 첫/끝 허용) ===
    maxTries = 200; tries = 0;
    while tries < maxTries
        tries = tries + 1;
        fixed = false;

        i_ = 1;
        while i_ <= numel(seq) - 2
            f3 = T.freq(seq(i_:i_+2));
            v3 = T.valence(seq(i_:i_+2));
            if ~(all(f3==f3(1)) || all(v3==v3(1))), i_ = i_ + 1; continue; end

            targetF = f3(1); targetV = v3(1);
            swapped = false;

            % 뒤쪽 탐색
            for j = (i_+3):numel(seq)
                if ~(T.freq(seq(j))==targetF || T.valence(seq(j))==targetV)
                    cand = seq; [cand(i_+1), cand(j)] = deal(cand(j), cand(i_+1));
                    bad = false; L = numel(cand); left = max(1,i_-1); right = min(L-2,i_+3);
                    for t = left:right
                        ff = T.freq(cand(t:t+2)); vv = T.valence(cand(t:t+2));
                        if all(ff==ff(1)) || all(vv==vv(1)), bad = true; break; end
                    end
                    if ~bad, seq = cand; swapped = true; fixed = true; break; end

                    cand = seq; [cand(i_+2), cand(j)] = deal(cand(j), cand(i_+2));
                    bad = false; L = numel(cand); left = max(1,i_-1); right = min(L-2,i_+3);
                    for t = left:right
                        ff = T.freq(cand(t:t+2)); vv = T.valence(cand(t:t+2));
                        if all(ff==ff(1)) || all(vv==vv(1)), bad = true; break; end
                    end
                    if ~bad, seq = cand; swapped = true; fixed = true; break; end
                end
            end


        % 앞쪽 탐색
            if ~swapped && i_ > 1
                for j = 1:(i_-1)
                    if ~(T.freq(seq(j))==targetF || T.valence(seq(j))==targetV)
                        cand = seq; [cand(i_+1), cand(j)] = deal(cand(j), cand(i_+1));
                        bad = false; L = numel(cand); left = max(1,i_-1); right = min(L-2,i_+3);
                        for t = left:right
                            ff = T.freq(cand(t:t+2)); vv = T.valence(cand(t:t+2));
                            if all(ff==ff(1)) || all(vv==vv(1)), bad = true; break; end
                        end
                        if ~bad, seq = cand; swapped = true; fixed = true; break; end

                        cand = seq; [cand(i_+2), cand(j)] = deal(cand(j), cand(i_+2));
                        bad = false; L = numel(cand); left = max(1,i_-1); right = min(L-2,i_+3);
                        for t = left:right
                            ff = T.freq(cand(t:t+2)); vv = T.valence(cand(t:t+2));
                            if all(ff==ff(1)) || all(vv==vv(1)), bad = true; break; end
                        end
                        if ~bad, seq = cand; swapped = true; fixed = true; break; end
                    end
                end
            end

                    % 내부 교환
            if ~swapped
                cand = seq; [cand(i_+1), cand(i_+2)] = deal(cand(i_+2), cand(i_+1));
                bad = false; L = numel(cand); left = max(1,i_-1); right = min(L-2,i_+3);
                for t = left:right
                    ff = T.freq(cand(t:t+2)); vv = T.valence(cand(t:t+2));
                    if all(ff==ff(1)) || all(vv==vv(1)), bad = true; break; end
                end
                if ~bad, seq = cand; fixed = true; end
            end

            i_ = i_ + 1;
        end

        if ~fixed, break; end
    end

    % ---- 여기서 블록 길이 확인 + 저장(필수) ----
    if numel(seq) ~= blockSize
        error('블록 %d의 trial 수(%d)가 blockSize(%d)와 다릅니다.', b, numel(seq), blockSize);
    end
    orderIdx{b} = seq(:);
end  % ← for b 루프 닫기

% 4) 최종 결합
finalOrder = vertcat(orderIdx{:});
assert(numel(finalOrder) == nBlocks*blockSize, ...
    'finalOrder 길이(%d)가 기대값(%d)과 다릅니다.', numel(finalOrder), nBlocks*blockSize);

Final      = T(finalOrder, :);
Final.designRow = finalOrder(:);   % 이 trial이 참조한 엑셀 Main row 번호

% 각 블록의 실제 길이에 맞춰 block 라벨 생성
blockLen   = cellfun(@numel, orderIdx).';          % [nBlocks×1]
Final.block = repelem((1:nBlocks)', blockLen);      % 길이 = height(Final)
numTrial    = height(Final);

% (선택) 무결성 체크
assert(sum(blockLen)==height(Final), 'blockLen 합과 Final 높이가 불일치합니다.');

for b=1:nBlocks
    fprintf('Block %d: total=%d (non-catch ~%d, catch ~%d)\n', ...
        b, blockLen(b), sum(~T.is_catch(orderIdx{b})), sum(T.is_catch(orderIdx{b})));
end

Final.sentence = arrayfun(@cleanSentence, Final.sentence);

% (assert는 DrawFormattedText 호출 직전에 두는 걸 권장)
% assert(Screen('WindowKind', w) > 0, 'Window handle lost before DFT (instruction).');

% key setting
keyA   = KbName('1');     % 비-캐치(true/자연스러움) 정답
keyL   = KbName('2');     % 캐치(어색함) 정답
keyEsc = KbName('ESCAPE');
recal  = KbName('r');     % (기존 유지)

% colours
bk = [0, 0, 0]; % 검정색

% output files
results = struct(); % 각 트라이얼의 반응시간, 정오답, 시선좌표 요약 등을 results.trial(t).RT 같은 식으로 저장할 그릇

% stimulus properties
ss = (dp.ppd*1) /2; % stimsize 
setsize = 10; % 검색 배열에 놓을 아이템 개수

%% EYELINK SETTINGS: Calibration, 실모드(dummymode=0)일 때 캘리브레이션 준비 
% Provide Eyelink with details about the graphics environment and perform some initializations. The information is returned
% in a structure that also contains useful defaults and control codes (e.g. tracker state bit and Eyelink key values).
if dummymode == 0
    el = EyelinkInitDefaults(w);
    el.backgroundcolour        = bgColor;
    el.calibrationtargetcolour = [0 0 0];

    % EyeLink 비프/타깃 삑소리 완전 비활성화
    el.feedbackbeep = 0;
    el.targetbeep   = 0;

    % 오디오 안 쓸 거면 이 줄은 삭제 또는 주석
    % InitializePsychSound(1);

    EyelinkUpdateDefaults(el);

    Eyelink('Command','screen_pixel_coords = %ld %ld %ld %ld', ...
        0, 0, dp.resolution(1)-1, dp.resolution(2)-1);
    Eyelink('Message','DISPLAY_COORDS %ld %ld %ld %ld', ...
        0, 0, dp.resolution(1)-1, dp.resolution(2)-1);
    Eyelink('Command','calibration_type = HV9');

    HideCursor;

    Eyelink('SetOfflineMode'); 
    WaitSecs(0.1);
    EyelinkDoTrackerSetup(el);
else
    ShowCursor('Arrow', w);
end

u8 = @(s) uint8(unicode2native(s,'UTF-8'));   % 한글 안전 출력 헬퍼

%% Instruction:
% === 텍스트 상태 가드 ===
assert(Screen('WindowKind', w) > 0, 'Window handle lost before DFT (instruction).');
Screen('Preference','TextRenderer', pickedRenderer);
Screen('TextFont', w, pickedFont);
Screen('TextSize', w, 40);
Screen('TextStyle', w, 0);

curSz = Screen('TextSize', w);
fprintf('[TEXT] font="%s" renderer=%d size=%d\n', pickedFont, pickedRenderer, curSz);
if pickedRenderer == 1
    tb = Screen('TextBounds', w, uint8(unicode2native('가나다 ABC 0123','UTF-8')));
else
    tb = Screen('TextBounds', w, '가나다 ABC 0123');
end
fprintf('[TEXT] bbox height=%d px, width=%d px\n', tb(4)-tb(2), tb(3)-tb(1));

Screen('TextSize', w, 50);  % 안내문 크기

instr = [ ...
    '<문장 읽기 과제>\n' ...
    '속으로 조용히 제시된 문장을 읽어주세요.\n' ...
    '시선은 평소에 글을 읽듯이 자연스럽게 움직이면 됩니다.\n' ...
    '문장을 읽는 동안 눈 깜박임은 최대한 자제해 주세요.\n' ...
    '문장에 전체 문맥과 맞지 않는 부분이 있는지 판단해주세요.\n' ...
    '문장이 의미적으로 자연스러우면 "1번" 버튼, 어색하면 "2번" 버튼을 눌러주세요.\n' ...
    '준비되면 Enter 버튼을 눌러 연습 시행을 시작하세요.' ...
    ];

DrawFormattedText(w, u8(instr), 'center','center', bk, 0, [], [], 2.5, 0);
Screen('Flip', w);                 % 안내 화면 온
WaitSecs(0.05);
if dummymode==0, Eyelink('Message','Instruction'); end

% === INSTRUCTION: ENTER/ESC 대기 (KbQueue 사용, 최소 설정) ===
kb.clear();      % 잔여 정리
kb.debounce();

allow = unique([KBD_SPACE, PAD_ENTER, {'return','escape'}]);  % 'return'은 일반 엔터 폴백
kb.only(allow);

pressed = kb.wait(allow, inf);   % 문자열 반환 가정
pressed = char(pressed);

kb.clear();
kb.debounce();

if strcmpi(pressed,'escape')
    if dummymode==0
        Eyelink('SetOfflineMode'); WaitSecs(0.5);
        Eyelink('CloseFile'); filetransfer(dirname, edfFile, dummymode);
    end
    cleanup; return;
end

Screen('Flip', w);  % 진행1 a

% [ANCHOR INSTR-WAIT END]

% === [ANCHOR AFTER-INSTR] Instruction 종료 직후: PTB 윈도우 유효성 가드 ===
if ~Screen('WindowKind', w)
    fprintf('[SAFE] PTB window handle lost after instruction. Reopening...\n');
    [w, rect, cx, cy] = ensureWindow(w, whichscreen, bgColor, 40, pickedFont, pickedRenderer);
end

% --- 다음 블록: 키 인덱스 캐시 ---
keyA   = KbName('1');
keyL   = KbName('2');
keyEsc = KbName('ESCAPE');
recal  = KbName('r');

%% === PRACTICE BLOCK (연습 시행; 본 시행 전에 10개 모두) ===
if nPractice > 0
    Eyelink('Message','PRACTICE_START N=%d', nPractice);

    % 연습 결과 저장용(분리 보관 권장)
    results.practice = struct();
    results.practice.RT      = nan(nPractice,1);
    results.practice.resp    = strings(nPractice,1);
    results.practice.acc     = false(nPractice,1);
    results.practice.isCatch = false(nPractice,1);

    for p = 1:nPractice
        % 창/텍스트 상태 보장
        [w, rect, cx, cy] = ensureWindow(w, whichscreen, bgColor, 40, pickedFont, pickedRenderer);

        % === ABC 좌표(픽셀) 미리 계산 ===
% (ABC_OFFSET_PX는 이미 위에서 계산되어 있음)
leftX  = ABC_OFFSET_PX;
yMid   = cy;

        Eyelink('Message','PRACTICE_TRIALID %d', p);
        Eyelink('Command','record_status_message "PRACTICE %d/%d"', p, nPractice);

        % 표시 좌표 계산(ABC는 본 시행과 동일 규칙)
        grid   = make_grid(dp.resolution, ss);
        nSlots = size(grid,2);
        if nSlots < setsize, setsize = nSlots; end
        grididx = randperm(nSlots, setsize);
        locs    = grid(:, grididx);
        tloc    = randperm(setsize, 1);
        tframe  = round([locs(1, tloc) - ss, locs(2, tloc) - ss, ...
                         locs(1, tloc) + ss, locs(2, tloc) + ss]);

        % 실모드면 드리프트 보정(선택: 연습 첫 trial에만 해도 됨)
        if dummymode == 0
    if p == 1
        % 혹시 이전에 Recording이 켜져 있었다면 안전하게 끄고 오프라인으로
        try Eyelink('StopRecording'); end
        Eyelink('SetOfflineMode');  WaitSecs(0.1);

        Eyelink('Message','PRACTICE_DRIFT_BEGIN %d', p);
        EyelinkDoDriftCorrection(el, cx, cy);
        Eyelink('Message','PRACTICE_DRIFT_END %d', p);

        % 드리프트 UI에서 복귀 직후 참가자 화면 버퍼를 깨끗이 복원
        Screen('FillRect', w, bgColor);
        Screen('Flip', w);
    end

    Eyelink('Command','clear_screen %d',7);
        
    % (선택) 중앙 십자 & 타깃 박스(호스트 시각화용)
    Eyelink('Command','draw_cross %d %d 15', cx, cy);
    Eyelink('Command','draw_box %d %d %d %d %d', ...
        tframe(1), tframe(2), tframe(3), tframe(4), 15);

    % === 호스트 오버레이: 왼쪽 ABC 위치 시각화 ===
hostDrawABC(leftX, yMid, dp.ppd, 'outer',1.5, 'cross',0.18, 'inner',0.18);
end

Eyelink('StartRecording'); WaitSecs(0.1);
% if dummymode == 0
%     HideCursor;                  % 실모드: 혹시 보이는 커서 다시 숨김
% else
%     ShowCursor('Arrow', w);      % 더미: 커서 보이도록
% end

outerDeg = 0.84;  
crossDeg = 0.10;  
innerDeg = 0.10;

% === [추가] 참가자 화면에 "왼쪽 ABC"를 먼저 실제로 표시 (초기 1회 Flip) ===
Screen('FillRect', w, bgColor);
draw_ABC_fixation(w, leftX, yMid, dp.ppd, ...
    'outer', outerDeg, 'cross', crossDeg, 'inner', innerDeg, ...
    'circlecolor', [0 0 0], 'crosscolor', bgColor);
[~, abcLeftOnsetP] = Screen('Flip', w);
Eyelink('Message','PRACTICE_ABC_LEFT_ONSET');

% === 좌측 ABC dwell check (실모드/더미 통합) ===
phase_log('BEGIN','LEFT_DWELL', 'trial', p);

holdSec  = 1.0;
winPx    = 75;
timeoutS = 10;  % 권장: 10초마다 재무장(=t0 초기화 후 10초를 다시 셈, 실험 흐름 유지). 누적 카운트가 필요하면 별도 로직 추가가 필요함.
t0       = GetSecs;
tEnter   = NaN;

if dummymode == 0
    Eyelink('Message','PRACTICE_LEFT_FIX_CHECK');
    eye_used = Eyelink('EyeAvailable'); 
    if eye_used == 2, eye_used = 0; end
    if eye_used < 0, eye_used = 0; end

    kb.only({'escape','r'});   % dwell 동안 ESC/R만 허용
    kb.debounce();             % 잔여 키 해제
    t0 = GetSecs;              % (디바운스 후) 타이머 리셋
    
    try KbQueueStop(-1);  catch, end
    try KbQueueFlush(-1); catch, end
    RestrictKeysForKbCheck([]); 
    DisableKeysForKbCheck([]);

    while true
        if ~Screen('WindowKind', w)
        [w, rect, cx, cy] = ensureWindow(w, whichscreen, bgColor, 40, pickedFont, pickedRenderer);
        Screen('FillRect', w, bgColor); Screen('Flip', w);
    else
        rect = Screen('Rect', w);
    end
        % --- ESC / 재캘 ---
        [down,~,kc] = KbCheck(-3);
        if down
            if kc(keyEsc)
                phase_log('END','LEFT_DWELL', 'trial', p, 'reason', 'ESC');
                kb.clear();
                Eyelink('SetOfflineMode'); WaitSecs(0.5);
                Eyelink('CloseFile'); filetransfer(dirname, edfFile, dummymode);
                cleanup; return;

            elseif kc(recal)
                Eyelink('Message','PRACTICE_LEFT_FIX_MANUAL_RECAL');
                Eyelink('SetOfflineMode'); WaitSecs(0.1);
                EyelinkDoTrackerSetup(el);
                Eyelink('StartRecording'); WaitSecs(0.1);
                eye_used = Eyelink('EyeAvailable'); 
                if eye_used == 2, eye_used = 0; end
                if eye_used < 0, eye_used = 0; end

                % 화면 복구
                Screen('FillRect', w, bgColor);
                draw_ABC_fixation(w, leftX, yMid, dp.ppd, ...
                    'outer', outerDeg, 'cross', crossDeg, 'inner', innerDeg, ...
                    'circlecolor', [0 0 0], 'crosscolor', bgColor);
                Screen('Flip', w);

                % 입력/타이머 재무장
                kb.only({'escape','r'});
                kb.debounce();
                tEnter = NaN; t0 = GetSecs;
                continue;
            end
        end

        % --- 타임아웃 → 단순 리셋(재캘 없음) ---
        if (GetSecs - t0) > timeoutS
            Eyelink('Message','PRACTICE_LEFT_FIX_TIMEOUT');
            fprintf('[LEFT DWELL] 10s timeout → 화면 리셋 + 타이머 재무장(재캘은 R 키 수동).\n');


            % 화면만 다시 그려주고
            Screen('FillRect', w, bgColor);
            draw_ABC_fixation(w, leftX, yMid, dp.ppd, ...
                'outer', outerDeg, 'cross', crossDeg, 'inner', innerDeg, ...
                'circlecolor', [0 0 0], 'crosscolor', bgColor);
            Screen('Flip', w);

            % 키 / 타이머 리셋 후 다시 dwell 기다림
            kb.only({'escape','r'});
            kb.debounce();
            tEnter = NaN; 
            t0     = GetSecs;
            continue;
        end


        % --- 시선 샘플 → dwell 판정 ---
        if Eyelink('NewFloatSampleAvailable') > 0
            evt = Eyelink('NewestFloatSample');
            x = evt.gx(eye_used+1); y = evt.gy(eye_used+1);
            if ~isnan(x) && ~isnan(y) && x~=el.MISSING_DATA && y~=el.MISSING_DATA
                inWin = abs(x-leftX)<=winPx && abs(y-yMid)<=winPx;
                if inWin
                    if isnan(tEnter), tEnter = GetSecs; end
                    if (GetSecs - tEnter) >= holdSec
                        Eyelink('Message','PRACTICE_LEFT_FIX_OK');
                        phase_log('END','LEFT_DWELL', 'trial', p, 'reason', 'OK');
                        break;
                    end
                else
                    tEnter = NaN;
                end
            end
        end
        WaitSecs(0.005);  % 루프 안에서만 짧게 쉬기
    end

    kb.clear();  % 정상 종료 후 한 번만 정리

else
    % === 더미모드: 마우스 “응시” 판정 ===
    ShowCursor('Arrow', w);
    RestrictKeysForKbCheck([]);                 
    while KbCheck(-3); WaitSecs(0.01); end

    holdSec  = 1.0;  winPx = 75;  timeoutS = 10;
    t0 = GetSecs;  tEnter = NaN;

    while true
        [down, ~, kc] = KbCheck(-3);
        if down
            if kc(keyEsc)
                phase_log('END','LEFT_DWELL', 'trial', p, 'reason', 'ESC');
                RestrictKeysForKbCheck([]); cleanup; return;
            elseif kc(recal)
                % 화면 복구 + 타이머 리셋
                Screen('FillRect', w, bgColor);
                draw_ABC_fixation(w, leftX, yMid, dp.ppd, ...
                    'outer', outerDeg, 'cross', crossDeg, 'inner', innerDeg, ...
                    'circlecolor', [0 0 0], 'crosscolor', bgColor);
                Screen('Flip', w);
                t0 = GetSecs; tEnter = NaN; safeKbReleaseWait;
                continue;
            end
        end

        % [ANCHOR P-LEFT-DUMMY] GetMouse 직전
% === (가드) PTB 윈도우 핸들 복구 ===
if ~Screen('WindowKind', w)
            [w, rect, cx, cy] = ensureWindow(w, whichscreen, bgColor, 40, pickedFont, pickedRenderer);
        end

        % 마우스 기반 “응시”
        [mx, my] = GetMouse;
        inWin = ~isnan(mx) && ~isnan(my) && ...
                abs(mx-leftX) <= winPx && abs(my-yMid) <= winPx;

        if inWin
            if isnan(tEnter), tEnter = GetSecs; end
            if (GetSecs - tEnter) >= holdSec
                phase_log('END','LEFT_DWELL', 'trial', p, 'reason', 'OK');
                break;
            end
        else
            tEnter = NaN;
        end

% 타임아웃: 화면 복구
        if (GetSecs - t0) > timeoutS
            % 화면 복구 및 리셋
            Screen('FillRect', w, bgColor);
            draw_ABC_fixation(w, leftX, yMid, dp.ppd, ...
                'outer', outerDeg, 'cross', crossDeg, 'inner', innerDeg, ...
                'circlecolor', [0 0 0], 'crosscolor', bgColor);
            Screen('Flip', w);
            fprintf('[LEFT DWELL DUMMY] 10s timeout → 화면 리셋, 타이머 재무장(R 키 수동 재캘 가능).\n');
            tEnter = NaN; t0 = GetSecs;
            safeKbReleaseWait;
        end

        WaitSecs(0.005);
    end
end

kb.clear();   % [MOVE] 루프 정상 종료 후 한 번만 해제/플러시

        % 문장 + 오른쪽 ABC
Screen('FillRect', w, bgColor);
thisSentP  = Practice.sentence(p);
thisSentU8 = uint8(unicode2native(char(thisSentP), 'UTF-8'));

% ← 문장 바이트 준비 '직후'에 오른쪽 ABC x좌표 계산
%    3 는 온점에서 오른쪽으로 3cm, outerDeg는 ABC 외곽 지름(도)과 동일 값 사용
rightX = computeRightABCx_fromDot(w, char(thisSentP), TEXT_LEFT_X, dp, 3, outerDeg);

DrawFormattedText(w, thisSentU8, TEXT_LEFT_X, 'center', bk, 0, [], [], 1);
draw_ABC_fixation(w, rightX, yMid, dp.ppd, ...
    'outer', outerDeg, 'cross', crossDeg, 'inner', innerDeg, ...
    'circlecolor', [0 0 0], 'crosscolor', bgColor);

[~, SearchOnsetP] = Screen('Flip', w);

% Eyelink 로그(Practice 전용 라벨/변수 유지)
msg = sprintf('PRACTICE_SENTENCE_ONSET %s_%s CATCH=%d', ...
    char(Practice.freq(p)), char(Practice.valence(p)), int32(Practice.is_catch(p)));
Eyelink('Message', msg);

% === PRACTICE: 오른쪽 ABC dwell 게이트 (필수) ===
phase_log('BEGIN','RIGHT_DWELL', 'trial', p);  % ← 연습이면 p, 본시행이면 i

holdSec  = 1.0;
winPx    = 75;
timeoutS = 10;
t0       = GetSecs;
tEnter   = NaN;

if dummymode == 0
    Eyelink('Message','PRACTICE_RIGHT_FIX_CHECK');
    eye_used = Eyelink('EyeAvailable'); 
    if eye_used == 2, eye_used = 0; end
    if eye_used < 0, eye_used = 0; end
    dwellPauseOffsetSec = 0;  % 재캘 시간 누적 (RightDwell 보정용)

    % dwell 동안 ESC/R만 허용 + 잔여 키 해제 + 타이머 리셋
    kb.only({'escape','r'});
    kb.debounce();
    t0 = GetSecs; tEnter = NaN;

    try KbQueueStop(-1);  catch, end
    try KbQueueFlush(-1); catch, end
    RestrictKeysForKbCheck([]); 
    DisableKeysForKbCheck([]);

     while true
        % --- ESC / 재캘 ---
        [down, ~, kc] = KbCheck(-3);
        if down
            if kc(keyEsc)
                phase_log('END','RIGHT_DWELL', 'trial', p, 'reason','ESC');
                kb.clear();
                Eyelink('SetOfflineMode'); WaitSecs(0.5);
                Eyelink('CloseFile'); filetransfer(dirname, edfFile, dummymode);
                cleanup; return;

            elseif kc(recal)
                Eyelink('Message','PRACTICE_RIGHT_FIX_MANUAL_RECAL');
                tPause = GetSecs;
                Eyelink('SetOfflineMode'); WaitSecs(0.1);
                EyelinkDoTrackerSetup(el);
                Eyelink('StartRecording'); WaitSecs(0.1);
                eye_used = Eyelink('EyeAvailable'); 
                if eye_used == 2, eye_used = 0; end 
                if eye_used < 0, eye_used = 0; end

                % 화면 복원(문장 + 오른쪽 ABC)
                Screen('FillRect', w, bgColor);
                DrawFormattedText(w, thisSentU8, TEXT_LEFT_X, 'center', bk, 0, [], [], 1);
                draw_ABC_fixation(w, rightX, yMid, dp.ppd, ...
                    'outer', outerDeg, 'cross', crossDeg, 'inner', innerDeg, ...
                    'circlecolor', [0 0 0], 'crosscolor', bgColor);
                Screen('Flip', w);

                % 시간 보정 + 입력 재무장 + 타이머 리셋
                dwellPauseOffsetSec = dwellPauseOffsetSec + (GetSecs - tPause);
                kb.only({'escape','r'});
                kb.debounce();
                tEnter = NaN; t0 = GetSecs;
                continue;
            end
        end

% === PRACTICE: 오른쪽 ABC dwell 게이트 (필수) ===
%%% PRACTICE: RIGHT_DWELL 응답 루프(1/2/ESC/R)

        % 타임아웃 → 단순 리셋(재캘 없음)
        if (GetSecs - t0) > timeoutS
            Eyelink('Message','PRACTICE_RIGHT_FIX_TIMEOUT');

            % 화면 복원만 하고
            Screen('FillRect', w, bgColor);
            DrawFormattedText(w, thisSentU8, TEXT_LEFT_X, 'center', bk, 0, [], [], 1);
            draw_ABC_fixation(w, rightX, yMid, dp.ppd, ...
                'outer', outerDeg, 'cross', crossDeg, 'inner', innerDeg, ...
                'circlecolor', [0 0 0], 'crosscolor', bgColor);
            Screen('Flip', w);

            % 타이머 리셋 후 다시 dwell 기다림
            t0     = GetSecs; 
            tEnter = NaN;
            continue;
        end

        % 시선 샘플 → dwell 판정
        if Eyelink('NewFloatSampleAvailable') > 0
            evt = Eyelink('NewestFloatSample');
            x = evt.gx(eye_used+1); y = evt.gy(eye_used+1);
            if ~isnan(x) && ~isnan(y) && x~=el.MISSING_DATA && y~=el.MISSING_DATA
                inWin = abs(x-rightX)<=winPx && abs(y-yMid)<=winPx;
                if inWin
                    if isnan(tEnter), tEnter = GetSecs; end
                    if (GetSecs - tEnter) >= holdSec
                        rightDwellSec = max(0, (GetSecs - tEnter) - dwellPauseOffsetSec);
                        Eyelink('Message','PRACTICE_RIGHT_FIX_OK');
                        Eyelink('Message','!V TRIAL_VAR RightDwellMs %d', round(rightDwellSec*1000));
                        phase_log('END','RIGHT_DWELL', 'trial', p, 'reason','OK', 'dwellMs', round(rightDwellSec*1000));
                        break;
                    end
                else
                    tEnter = NaN;
                end
            end
        end
        WaitSecs(0.005);
    end

    kb.clear();  % 정상 종료 후 한 번만 정리

else
    % === 더미모드: 마우스 “응시” 판정 ===
    ShowCursor('Arrow', w);

    % (더미에서도 동일한 위생: ESC/R 체크, 루프 내부 WaitSecs)
    while true
        [down,~,kc] = KbCheck(-3);
        if down
            if kc(keyEsc)
                phase_log('END','RIGHT_DWELL', 'trial', p, 'reason','ESC');
                RestrictKeysForKbCheck([]); cleanup; return;
            elseif kc(recal)
                % 화면 복구 + 타이머 리셋 (문장 + 오른쪽 ABC)
                Screen('FillRect', w, bgColor);
                DrawFormattedText(w, thisSentU8, TEXT_LEFT_X, 'center', bk, 0, [], [], 1);
                draw_ABC_fixation(w, rightX, yMid, dp.ppd, ...
                    'outer', outerDeg, 'cross', crossDeg, 'inner', innerDeg, ...
                    'circlecolor', [0 0 0], 'crosscolor', bgColor);
                Screen('Flip', w);
                t0 = GetSecs; tEnter = NaN; safeKbReleaseWait;
                continue;
            end
        end

        % [ANCHOR P-RIGHT-DUMMY] GetMouse 직전
    % === (가드) PTB 윈도우 핸들 복구 ===
    if ~Screen('WindowKind', w)
        fprintf('[SAFE] PTB window handle lost before GetMouse — reopening...\n');
        [w, rect, cx, cy] = ensureWindow(w, whichscreen, bgColor, 40, pickedFont, pickedRenderer);
    end

        % 기존 마우스 dwell 판정 
        [mx, my, w, rect, cx, cy] = safeGetMouse(w, whichscreen, bgColor, pickedFont, pickedRenderer, 40);
        inWin = ~isnan(mx) && ~isnan(my) && ...
                abs(mx-rightX) <= winPx && abs(my-yMid) <= winPx;
        if inWin
            if isnan(tEnter), tEnter = GetSecs; end
            if (GetSecs - tEnter) >= holdSec
                rightDwellSec = (GetSecs - tEnter);
                Eyelink('Message','PRACTICE_RIGHT_FIX_OK DUMMY');
                Eyelink('Message','!V TRIAL_VAR RightDwellMs %d', round(rightDwellSec*1000));
                phase_log('END','RIGHT_DWELL', 'trial', p, 'reason','OK', 'dwellMs', round(rightDwellSec*1000));
                break;
            end
        else
            tEnter = NaN;
        end

         if (GetSecs - t0) > timeoutS
            % 화면 재그림
            Screen('FillRect', w, bgColor);
            DrawFormattedText(w, thisSentU8, TEXT_LEFT_X, 'center', bk, 0, [], [], 1);
            draw_ABC_fixation(w, rightX, yMid, dp.ppd, ...
                'outer', outerDeg, 'cross', crossDeg, 'inner', innerDeg, ...
                'circlecolor', [0 0 0], 'crosscolor', bgColor);
            Screen('Flip', w);
            t0 = GetSecs; tEnter = NaN;
        end
        WaitSecs(0.005);
    end
end

% === PROMPT 화면 띄우고, 여기부터 1/2/ESC/R 받기 (Practice) ===
Screen('FillRect', w, bgColor);
Screen('TextSize', w, 40);
promptTxt = ['방금 제시된 문장이 의미적으로 자연스러웠나요?\n' ...
             '자연스러웠다: 1, 어색했다: 2'];
DrawFormattedText(w, u8(promptTxt), 'center','center', bk, 70, [], [], 1.4);
[~, PromptOnsetP] = Screen('Flip', w);
Eyelink('Message','PRACTICE_PROMPT_ONSET');

% --- 기존 safeKbReleaseWait / RestrictKeysForKbCheck 제거 ---
kbPad.clear();  kbKb.clear();
kbPad.debounce();  kbKb.debounce();

kbPad.only([RESP1_KEYS, RESP2_KEYS]);
kbKb.only([RESP1_KEYS, RESP2_KEYS, KBD_ADMIN]);   % KBD_ADMIN = {'escape','r'}

pauseOffsetSec = 0;
name = ""; tResp = NaN;

while name==""
    % keypad: 1/2만
    [nP, tP] = kbPad.wait([RESP1_KEYS, RESP2_KEYS], 0.02);
    % keyboard: 1/2 + ESC/R
    [nK, tK] = kbKb.wait([RESP1_KEYS, RESP2_KEYS, KBD_ADMIN], 0.02);

    % 1) ESC: 즉시 종료
     if (nK~="") && strcmpi(nK,'escape')
        if dummymode==0
            Eyelink('SetOfflineMode'); WaitSecs(0.5);
            Eyelink('CloseFile'); filetransfer(dirname, edfFile, dummymode);
        end
        cleanup; return;
    end

    % 2) 재캘(R): 장비 재설정 → 화면 복구 → 시간 보정 → 다시 대기
    if (nK~="") && strcmpi(nK,'r') && dummymode==0
    tPauseStart = GetSecs;
    Eyelink('Message','!!! PRACTICE_RECAL !!!');
    try Eyelink('StopRecording'); end
    Eyelink('SetOfflineMode');  WaitSecs(0.1);
    EyelinkDoTrackerSetup(el);
    Eyelink('StartRecording');  WaitSecs(0.1);

    % 프롬프트 다시 그리기
    Screen('FillRect', w, bgColor);
    DrawFormattedText(w, u8(promptTxt), 'center','center', bk, 70, [], [], 1.4);
    Screen('Flip', w);

    % RT 보정 및 큐 초기화/재설정
    pauseOffsetSec = pauseOffsetSec + (GetSecs - tPauseStart);
    kbPad.clear();  kbKb.clear();
    kbPad.debounce(); kbKb.debounce();
    kbPad.only([RESP1_KEYS, RESP2_KEYS]);
    kbKb.only([RESP1_KEYS, RESP2_KEYS, KBD_ADMIN]);  % KBD_ADMIN={'escape','r'}

    continue;  % 다시 대기
    end
    
    % 3) 1 또는 2 응답
if name==""
    if nP~=""
        name  = nP;  tResp = tP;                   % keypad 1/2
    elseif (nK~="") && ~any(strcmpi(char(nK), KBD_ADMIN))
        name  = nK;  tResp = tK;                   % keyboard 1/2 (ESC/R 제외)
    end
end

n = char(string(name));   % '1','2' 등으로 강제 (string/char 혼선 방지)

% 키 분류(환경에 RESP* 세트가 있으면 사용, 없으면 '1'/'2' 폴백)
if exist('RESP1_KEYS','var') && ~isempty(RESP1_KEYS)
    is1 = any(strcmpi(n, RESP1_KEYS));
else
    is1 = strcmpi(n,'1');
end
if exist('RESP2_KEYS','var') && ~isempty(RESP2_KEYS)
    is2 = any(strcmpi(n, RESP2_KEYS));
else
    is2 = strcmpi(n,'2');
end

if is1 || is2
    % 시간 소스 안전 선택: tResp → t → GetSecs
    if exist('tResp','var') && ~isempty(tResp) && isfinite(tResp)
        tUse = tResp;
    elseif exist('t','var') && ~isempty(t) && isfinite(t)
        tUse = t;
    else
        tUse = GetSecs;
    end

    rtMs  = round( (tUse - PromptOnsetP - pauseOffsetSec) * 1000 );  % 본실험은 PromptOnset 사용
    isCat = Practice.is_catch(p) ~= 0;
    if is1, resp='1'; acc = ~isCat; else, resp='2'; acc = isCat; end

    results.practice.RT(p,1)      = rtMs/1000; 
    results.practice.resp(p,1)    = string(resp);
    results.practice.acc(p,1)     = acc;
    results.practice.isCatch(p,1) = isCat;

    Eyelink('Message', sprintf('PRACTICE_RESPONSE %s %d ACC=%d IS_CATCH=%d', ...
        char(resp), int32(rtMs), int32(acc), int32(isCat)));

    % (선택) 피드백
    Screen('FillRect', w, bgColor);
    if acc, fbText='맞았습니다'; fbColor=[0 200 0]; fbTag='CORRECT';
    else,  fbText='틀렸습니다'; fbColor=[230 0 0]; fbTag='INCORRECT'; end
    Eyelink('Message',['PRACTICE_FEEDBACK_ON ' fbTag]);
    Screen('TextSize', w, 50);
    DrawFormattedText(w, u8(fbText), 'center','center', fbColor, 0, [], [], 1.2);
    Screen('Flip', w);  WaitSecs(1.0);
    Screen('FillRect', w, bgColor); Screen('Flip', w);
    Eyelink('Message','PRACTICE_FEEDBACK_OFF');

    kb.clear();   % 이 블록 종료 시 한 번만 정리
    break;
end
        end

        WaitSecs(0.1);
        Eyelink('StopRecording');
        Eyelink('Message','PRACTICE_TRIAL_RESULT 0');

        % 연습 체크포인트 (예: 5 trial마다)
if mod(p,5)==0
    try
        save(fullfile(dirname,'practice_checkpoint.mat'),'results','dp','p');
    end
end
% === 루프 끝 (end of for p = 1:nPractice) ===
end  % <-- 연습 trial 루프 닫힘

% === PRACTICE BLOCK 종료 ===
Eyelink('Message','PRACTICE_END');

% [ANCHOR PRACTICE->MAIN REC SYNC]
% === PRACTICE BLOCK 끝난 직후, 본 시행 진입하기 직전 ===
if dummymode==0
    try Eyelink('StopRecording'); end
    Eyelink('SetOfflineMode'); WaitSecs(0.1);
end

% === [SAFETY] 본시행 진입 전: 입력/커서 완전 리셋 ===
kb.clear();          % 잔여 큐/제한 정리
kb.debounce();       % 눌려 있던 키 완전 해제
if dummymode~=0
    ShowCursor('Arrow', w);
else
    HideCursor;
end

    % 본 시행 안내 화면
    Screen('FillRect', w, bgColor);
    DrawFormattedText(w, u8('연습시행이 끝났습니다.\n이제 본 실험을 시작합니다.\nEnter 버튼을 누르면 시작합니다.'), ...
                      'center','center', bk, 60, [], [], 1.4);
    Screen('Flip', w);

% === INSTRUCTION: enter/ESC 대기 ===
kb.clear();            % 잔여 제한/큐 정리
kb.debounce();         % 눌린 키가 남아 있으면 전부 떼질 때까지 대기

allow = unique([KBD_SPACE, PAD_ENTER, {'return','escape'}]);  % 'return'은 일반 엔터 폴백
kb.only(allow);
[name] = kb.wait(allow, inf);  % 블로킹 대기

kb.clear();            % 제한 해제 / 큐 플러시
kb.debounce();         % (권장) 다음 단계로 넘어가기 전 키 떼기

if strcmpi(name,'escape')
    if dummymode==0
        Eyelink('SetOfflineMode'); WaitSecs(0.5);
        Eyelink('CloseFile'); filetransfer(dirname, edfFile, dummymode);
    end
    cleanup; return;
end

% 여기 도달했다면 space 또는 (키패드/메인) enter였다 → 화면 정리 후 다음 단계
Screen('FillRect', w, bgColor);   % ← 이 줄 추가
Screen('Flip', w);                % ← 이제 진짜로 클리어됨

% === [SAFETY] 본시행 진입 전 화면/좌표/링크 재점검 ===
if ~Screen('WindowKind', w)
    fprintf('[SAFE] PTB window handle lost — reopening...\n');
    [w, rect, cx, cy] = ensureWindow(w, whichscreen, bgColor, 40, pickedFont, pickedRenderer);
else
    rect = Screen('Rect', w);          % ← 현재 열린 창 기준으로 다시 얻기
    [cx, cy] = RectCenter(rect);
end

% grid/ABC/tframe 계산이 모두 동일한 기준을 보도록 동기화
dp.resolution = [rect(3) rect(4)];

% (실모드만) 호스트 좌표계도 현재 창 크기로 동기화
if exist('dummymode','var') && dummymode==0
    Eyelink('Command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, rect(3)-1, rect(4)-1);
    Eyelink('Message','DISPLAY_COORDS %ld %ld %ld %ld',        0, 0, rect(3)-1, rect(4)-1);

    % 연결 상태만 확인 (중간 재초기화/Initialize는 EDF 리스크가 있어서 생략)
    if Eyelink('IsConnected') ~= 1
        warning('[EL] Link down after practice. 재캘(R) 또는 세션 재시작 권장.');
    end
end
end

%% Main experiment: 각 트라이얼을 시작·표시·기록하는 전체 흐름
%numTrial = 10; % 총 10회 반복. i가 현재 트라이얼 번호

% === 키보드 베이스라인 초기화 (루프 밖, 파일 상단 추천) === 
KbName('UnifyKeyNames');
RestrictKeysForKbCheck([]);   % 혹시 남은 잠금 해제
DisableKeysForKbCheck([]);    % 혹시 비활성 목록 초기화

keyA   = KbName('1');
keyL   = KbName('2');
keyEsc = KbName('ESCAPE');
recal  = KbName('r');

% [A0] ESC 워치독: RestrictKeysForKbCheck와 무관하게 ESC만 따로 감시
try KbQueueRelease(-1); catch, end
escList = zeros(1,256); escList(keyEsc) = 1;

% Block 시작 trial 번호 사전 계산
blockFirstTrial = accumarray(Final.block, (1:height(Final))', [], @min);

for i = 1:numTrial
    if i == 1
        blockStart = true;                 % 첫 trial은 무조건 블록 시작으로 처리
    else
        blockStart = (Final.block(i) ~= Final.block(i-1));
    end

    % --- 트라이얼 시작 시 창 유효성 보장 ---
    [w, rect, cx, cy] = ensureWindow(w, whichscreen, bgColor, 40, pickedFont, pickedRenderer);

    % >>> 여기 추가 <<<
    checkLinkOrAbort('TRIAL_START', sprintf('trial %d', i), ...
        dirname, edfFile, dummymode, w, bgColor, pickedFont, pickedRenderer);

     % 여기서 rect를 동기화 (이 줄 바로 아래 삽입)
    rect = Screen('Rect', w);
    dp.resolution = [rect(3) rect(4)];   % grid, ABC 오프셋, tframe 모두 같은 기준

    % --- ABC 좌표(픽셀) 계산 (LEFT/RIGHT에서 같이 씀) ---
leftX  = ABC_OFFSET_PX;
yMid   = cy;

% (5) 블록 시작이면 중앙 '+' 1초 표시 — 키 입력 요구 없음
if blockStart
    Screen('FillRect', w, bgColor);
    colPlus = [0 0 0];
    Screen('DrawLine', w, colPlus, cx-25, cy, cx+25, cy, 5);
    Screen('DrawLine', w, colPlus, cx, cy-25, cx, cy+25, 5);
    [~, plusOnset] = Screen('Flip', w);  % plusOnset 안 쓰면 출력만 하고 무시해도 됨
    Eyelink('Message','BLOCK_FIXATION_ON %d', Final.block(i));

% === [PATCH] 블록 시작 '+' 화면에서 fresh 키 입력 대기 ===
% (이전 run/phase에서 눌려 있던 키 상태가 남아있어도 여기서 정리)

 % 키는 안 기다리고, +를 1초 동안만 보여줌
    WaitSecs(1.0);
end

    fprintf('Trial %d 시작, WindowKind=%d\n', i, Screen('WindowKind', w));
    % Sending a 'TRIALID' message to mark the start of a trial in Data Viewer. 
    % This is different than the start of recording message START that is logged when the trial recording begins. 
    % The viewer will not parse any messages, events, or samples, that exist in the data file prior to this message.       
    Eyelink('Message', 'TRIALID %d', i); % TRIALID: Data Viewer에서 여기부터 '트라이얼 i'로 인식(이전 메시지·샘플은 파싱 대상X)
    % This supplies the title at the bottom of the eyetracker display
    Eyelink('Command', 'record_status_message "TRIAL %d/%d"', i, numTrial); % record_status_message: 호스트 화면 하단 상태줄에 "TRIAL i/numTrial" 표시

% (1) 이번 트라이얼 자극 좌표/프레임 계산  
grid   = make_grid(dp.resolution, ss);
nSlots = size(grid,2);
if nSlots < setsize
    warning('Grid slots (%d) < setsize (%d). Using %d.', nSlots, setsize, nSlots);
    setsize = nSlots;
end
grididx              = randperm(nSlots, setsize);
locs                 = grid(:, grididx);
tloc                 = randperm(setsize, 1);
results.target(i, :) = locs(:, tloc)';
tframe               = round([locs(1, tloc) - ss, locs(2, tloc) - ss, ...
                              locs(1, tloc) + ss, locs(2, tloc) + ss]);

% (2) 드리프트 보정(모든 모드에서 실행) → (3) 호스트 오버레이
% === 자동 드리프트 코렉션 조건 추가 ===
if dummymode==0
    currentBlock = Final.block(i);
    firstTrialInBlock = NaN;
    doDrift = false;

    if ~isempty(blockFirstTrial) && currentBlock>=1 && currentBlock<=numel(blockFirstTrial)
        firstTrialInBlock = blockFirstTrial(currentBlock);
        if ~isnan(firstTrialInBlock) && firstTrialInBlock>0
            since = i - firstTrialInBlock;              % 0,1,2,...
            doDrift = (since==0) || (mod(since,3)==0);  % 블록 첫 trial 및 매 3 trial
        end
    end

% --- 진단 로그(원인 확인용) ---
fprintf('[DRIFT] trial=%d | block=%d | first=%d | since=%d | doDrift=%d\n', ...
    i, currentBlock, firstTrialInBlock, i-firstTrialInBlock, doDrift);

% === [임시 패치] 드리프트 코렉션 비활성화 → 다시 활성화 ===
if doDrift
    % (1) 드리프트 시작 메시지 로그
    Eyelink('Message','AUTO_DRIFT_BEGIN %d', i);

    % (2) 오프라인 모드 전환 + 잠깐 대기
    Eyelink('SetOfflineMode');
    WaitSecs(0.05);  % 0.1로 늘려도 됨

    % (3) 키 상태 초기화 (핵심: 이전 키 입력 여파 제거)
    Eyelink('Flushkeybuttons');   % Eyelink 키버퍼 비우기
    KbReleaseWait(-1);            % PTB 키보드에서 모든 키가 올라갈 때까지 대기

    % (4) Drift correction 실행
    status = EyelinkDoDriftCorrection(el, cx, cy, 1, 1);
    fprintf('[DRIFT] trial=%d | status=%d\n', i, status);

    % (5) 드리프트 종료 메시지 + 상태 남기기
    Eyelink('Message','AUTO_DRIFT_END %d STATUS %d', i, status);
end

end

% (드리프트 코렉션을 했든 안 했든) 호스트 오버레이는 항상 그림
Eyelink('Command','clear_screen %d',7);
Eyelink('Command','draw_cross %d %d 15', cx, cy);
Eyelink('Command','draw_box %d %d %d %d 15', ...
    tframe(1), tframe(2), tframe(3), tframe(4));

% (4) 이제부터 녹화 시작
Eyelink('StartRecording');
WaitSecs(0.2); % StartRecording 직후 200ms 버퍼 안정화

% (유지) 눈 사용 정보
eye_used = Eyelink('EyeAvailable');   % 0=left, 1=right, 2=both, -1=none
if eye_used == 2
    eye_used = 0;                      % 왼눈 기본
elseif eye_used == -1
    Eyelink('Message','EYE_UNAVAILABLE_FALLBACK');
end

% === 입력 문자열 가드 ===
thisSent = Final.sentence(i);
if isstring(thisSent), sChar = char(thisSent); else, sChar = thisSent; end
fprintf('Trial %d | is_catch=%d | len=%d | preview="%s"\n', ...
    i, Final.is_catch(i), length(sChar), sChar(1:min(end,80)));

rawS   = Final.sentence(i);
if ismissing(rawS), rawS = ""; end
sClean = cleanSentence(rawS);
lenS   = strlength(strtrim(sClean));
if lenS >= 1
    thisSent = sClean;
else
    fprintf('[WARN] empty-ish after cleanup (i=%d)\n', i);
    thisSent = " ";
end
thisSentU8 = uint8(unicode2native(char(thisSent), 'UTF-8'));

rightX = computeRightABCx_fromDot(w, char(thisSent), TEXT_LEFT_X, dp, 3, outerDeg);

% === 그리기 직전 디버그(여기!) ===
if i <= 3
    % 안전한 80자 미리보기
    sPrev = char(thisSent);
    if isempty(sPrev)
        sPrev = '(empty)';
    else
        sPrev = sPrev(1:min(80, numel(sPrev)));
    end

    fprintf('DBG[%d] catch=%d | rawLen=%d | cleanLen=%d | use="%s"\n', ...
        i, Final.is_catch(i), strlength(rawS), strlength(sClean), sPrev);

    cp = double(char(thisSent));
    fprintf('   codepoints(1..20)=['); fprintf('%d ', cp(1:min(end,20))); fprintf(']\n');
end

fprintf('DBG trial=%d | catch=%d | show="%s"\n', i, Final.is_catch(i), char(thisSent));

% === 왼쪽 ABC -> 1초 후 OFF -> 문장+오른쪽 ABC ===
% 5cm 오프셋 좌표 계산
leftX  = ABC_OFFSET_PX;
rightX = computeRightABCx_fromDot(w, char(thisSent), TEXT_LEFT_X, dp, 3, outerDeg);
yMid   = cy;

outerDeg = 0.84;  
crossDeg = 0.10;  
innerDeg = 0.10;

% === 커서 표시 정책: 실모드=숨김, 더미모드=표시 ===
if dummymode ~= 0
    ShowCursor('Arrow', w);   % dummy 모드에서만 보이게
else
    HideCursor;               % real 모드에서는 계속 숨김
end

kb.only({'escape','r'});                        % ESC/R만 허용
kb.debounce();                                  % 눌려 있던 키 완전 해제
t0 = GetSecs; tEnter = NaN;                     % (필요시) 타이머 리셋

%% [PATCH MAIN-LEFT-REAL-ABC] LEFT dwell (REAL + DUMMY 통합 교체본)
if dummymode == 0
    Eyelink('Message','LEFT_FIX_CHECK');
    Eyelink('Command','clear_screen %d',7);
    hostDrawABC(leftX, yMid, dp.ppd, 'outer',1.5, 'cross',0.18, 'inner',0.18);
    Eyelink('Command','draw_box %d %d %d %d %d', ...
        tframe(1), tframe(2), tframe(3), tframe(4), 15);

    holdSec  = 1.0;   % 연속 응시 1초
    winPx    = 75;    % 허용 반경(px)
    timeoutS = 10;    % 최대 10초 안에 못 잡으면 타임아웃

    % --- 키보드/큐 정리 ---
    try KbQueueStop(-1);  catch, end
    try KbQueueFlush(-1); catch, end
    RestrictKeysForKbCheck([]); 
    DisableKeysForKbCheck([]);

    % --- 사용 눈 결정 ---
    eye_used = Eyelink('EyeAvailable'); 
    if eye_used == 2, eye_used = 0; end
    if eye_used < 0, eye_used = 0; end

    % === 여기서 실제 왼쪽 ABC를 "한 번" 그려줌 ===
    Screen('FillRect', w, bgColor);
    draw_ABC_fixation(w, leftX, yMid, dp.ppd, ...
        'outer', outerDeg, 'cross', crossDeg, 'inner', innerDeg, ...
        'circlecolor', [0 0 0], 'crosscolor', bgColor);
    [~, abcLeftOnset] = Screen('Flip', w);
    Eyelink('Message','ABC_LEFT_ONSET_MAIN');
    fprintf('[DBG MAIN-LEFT-ABC] trial=%d | ABC_LEFT_ONSET_MAIN @%.3f\n', i, abcLeftOnset);

    % dwell 타이머 초기화
    t0     = GetSecs;   % 타임아웃 기준
    tEnter = NaN;       % 윈도우에 처음 들어온 시점

    while true
        % === 링크 상태 체크: 끊기면 EDF 정리 후 종료 ===
        checkLinkOrAbort('MAIN_LEFT_DWELL', sprintf('trial %d LEFT_DWELL', i), ...
            dirname, edfFile, dummymode, w, bgColor, pickedFont, pickedRenderer);

        % === (윈도우 가드) 윈도우가 닫혔으면 다시 열고 ABC 재그리기 ===
        if ~Screen('WindowKind', w)
            [w, rect, cx, cy] = ensureWindow(w, whichscreen, bgColor, 40, pickedFont, pickedRenderer);
            Screen('FillRect', w, bgColor);
            draw_ABC_fixation(w, leftX, yMid, dp.ppd, ...
                'outer', outerDeg, 'cross', crossDeg, 'inner', innerDeg, ...
                'circlecolor', [0 0 0], 'crosscolor', bgColor);
            Screen('Flip', w);
            t0     = GetSecs;
            tEnter = NaN;
        else
            rect = Screen('Rect', w);
        end

        % --- ESC / 재캘 체크 ---
        [down,~,kc] = KbCheck(-3);
        if down
            if kc(keyEsc)
                kb.clear();
                try Eyelink('StopRecording'); end           
                Eyelink('SetOfflineMode'); WaitSecs(0.5);
                Eyelink('CloseFile'); filetransfer(dirname, edfFile, dummymode);
                cleanup; return;
            end

            if kc(recal)
                Eyelink('Message','!!! RECAL (LEFT ABC) !!!');
                tPause = GetSecs;

                try Eyelink('StopRecording'); end 
                Eyelink('SetOfflineMode'); WaitSecs(0.1);
                EyelinkDoTrackerSetup(el);
                Eyelink('StartRecording'); WaitSecs(0.1);

                % 호스트 오버레이 복구
                Eyelink('Command','clear_screen %d',7);
                Eyelink('Command','draw_cross %d %d 15', cx, cy);
                Eyelink('Command','draw_box %d %d %d %d %d', ...
                    tframe(1), tframe(2), tframe(3), tframe(4), 15);

                % 참가자 화면 복구: 다시 왼쪽 ABC 표시
                Screen('FillRect', w, bgColor);
                draw_ABC_fixation(w, leftX, yMid, dp.ppd, ...
                    'outer', outerDeg, 'cross', crossDeg, 'inner', innerDeg, ...
                    'circlecolor', [0 0 0], 'crosscolor', bgColor);
                Screen('Flip', w);

                kb.only({'escape','r'});
                kb.debounce();
                tEnter = NaN;
                t0     = GetSecs;
                continue;
            end
        end

        % --- 시선 샘플 → dwell 판정 ---
        if Eyelink('NewFloatSampleAvailable') > 0
            evt = Eyelink('NewestFloatSample');
            x = evt.gx(eye_used+1); 
            y = evt.gy(eye_used+1);
            if ~isnan(x) && ~isnan(y) && x~=el.MISSING_DATA && y~=el.MISSING_DATA
                inWin = abs(x-leftX)<=winPx && abs(y-yMid)<=winPx;
                if inWin
                    if isnan(tEnter), tEnter = GetSecs; end
                    if (GetSecs - tEnter) >= holdSec
                        Eyelink('Message','LEFT_FIX_OK');
                        break;  % dwell 충족 → 루프 탈출
                    end
                else
                    tEnter = NaN;
                end
            end
        end

        % --- 타임아웃: ABC 다시 보여주고 dwell 재시작 ---
        if (GetSecs - t0) > timeoutS
            Eyelink('Message','LEFT_FIX_TIMEOUT');

            Screen('FillRect', w, bgColor);
            draw_ABC_fixation(w, leftX, yMid, dp.ppd, ...
                'outer', outerDeg, 'cross', crossDeg, 'inner', innerDeg, ...
                'circlecolor', [0 0 0], 'crosscolor', bgColor);
            Screen('Flip', w);

            t0     = GetSecs;
            tEnter = NaN;
            continue;
        end

        WaitSecs(0.005);   % CPU 과부하 방지
    end

    kb.clear();   % LEFT_FIX_OK로 정상 종료 후 한 번만

else
    % === 더미모드: 마우스 “응시” 판정 ===
    RestrictKeysForKbCheck([]);
    while KbCheck(-3); WaitSecs(0.01); end

    holdSec  = 1.0;
    winPx    = 75;
    timeoutS = 10;
    t0       = GetSecs;
    tEnter   = NaN;

    % 왼쪽 ABC 초기 표시
    Screen('FillRect', w, bgColor);
    draw_ABC_fixation(w, leftX, yMid, dp.ppd, ...
        'outer', outerDeg, 'cross', crossDeg, 'inner', innerDeg, ...
        'circlecolor', [0 0 0], 'crosscolor', bgColor);
    Screen('Flip', w);

    while true
        % --- ESC / 재캘(더미에선 재캘은 단순 리셋 용도) ---
        [down,~,kc] = KbCheck(-3);
        if down && kc(keyEsc)
            fprintf('[ESC] LEFT-REAL pressed\n'); 
            RestrictKeysForKbCheck([]); cleanup; return;
        end

        if down && kc(recal)
            Screen('FillRect', w, bgColor);
            draw_ABC_fixation(w, leftX, yMid, dp.ppd, ...
                'outer', outerDeg, 'cross', crossDeg, 'inner', innerDeg, ...
                'circlecolor', [0 0 0], 'crosscolor', bgColor);
            Screen('Flip', w);
            t0 = GetSecs; tEnter = NaN; safeKbReleaseWait;
            continue;
        end

        % 윈도우 가드
        if ~Screen('WindowKind', w)
            [w, rect, cx, cy] = ensureWindow(w, whichscreen, bgColor, 40, pickedFont, pickedRenderer);
            Screen('FillRect', w, bgColor);
            draw_ABC_fixation(w, leftX, yMid, dp.ppd, ...
                'outer', outerDeg, 'cross', crossDeg, 'inner', innerDeg, ...
                'circlecolor', [0 0 0], 'crosscolor', bgColor);
            Screen('Flip', w);
            t0 = GetSecs; tEnter = NaN;
        end

        % 마우스 “응시” 판정
        [mx, my] = GetMouse;
        inWin = ~isnan(mx) && ~isnan(my) && ...
                abs(mx-leftX) <= winPx && abs(my-yMid) <= winPx;

        if inWin
            if isnan(tEnter), tEnter = GetSecs; end
            if (GetSecs - tEnter) >= holdSec
                break;
            end
        else
            tEnter = NaN;
        end

        if (GetSecs - t0) > timeoutS
            Screen('FillRect', w, bgColor);
            draw_ABC_fixation(w, leftX, yMid, dp.ppd, ...
                'outer', outerDeg, 'cross', crossDeg, 'inner', innerDeg, ...
                'circlecolor', [0 0 0], 'crosscolor', bgColor);
            Screen('Flip', w);
            t0 = GetSecs; tEnter = NaN;
        end

        WaitSecs(0.005);
    end
end

% === [NEW] 한 줄 문장 단어별 ROI 계산 + EDF IAREA 메시지 기록 ===
% (TRIALID 보낸 직후, SearchOnset 찍기 전에 넣으세요)
% 현재 trial 문장: thisSent  (이미 위에서 정제된 string)
sentChar = char(thisSent);
Screen('Preference','TextRenderer', pickedRenderer);
Screen('TextFont',  w, pickedFont);
curSz = Screen('TextSize', w); %#ok<NASGU>

% === [NEW] 한 줄 문장 단어별 ROI 계산 + EDF IAREA 메시지 기록 ===
[wordRects, words, lineRect] = makeWordROIsSingleLine(w, cx, cy, sentChar, TEXT_LEFT_X);
results.words{i,1}     = words;
results.wordRects{i,1} = wordRects;
results.lineRect(i, :) = lineRect;

% === IAREA 좌표 안전화(정수/화면내/유효폭 보장) ===
scrW = dp.resolution(1);
scrH = dp.resolution(2);
safeRect = @(r) deal( ...
    max(0, min(scrW-1, round(r(1)))), ...
    max(0, min(scrH-1, round(r(2)))), ...
    max(0, min(scrW-1, round(r(3)))), ...
    max(0, min(scrH-1, round(r(4)))) );

% 기존 IAREA 전송 루프를 아래처럼 수정
Eyelink('Message','!V CLEAR IAREAS');
for k = 1:numel(words)
    [L,T,R,B] = safeRect(wordRects(k,:));
    if R <= L, R = min(scrW-1, L+1); end
    if B <= T, B = min(scrH-1, T+1); end
    Eyelink('Message','!V IAREA RECTANGLE %d W%d %d %d %d %d', k, k, L, T, R, B);
end

% === 안전 체크: 화면 폭 넘침 경고 ===
try
    scrW = rect(3);
    if wordRects(end, 3) > scrW - 5
        fprintf('[WARN] Sentence may overflow horizontally at trial %d (R=%d > scrW=%d)\n', ...
                i, wordRects(end,3), scrW);
    end
catch
end

% === (선택) 표적/포스트표적 라벨 ===
% Final.targetWordIdx(i)가 1..numel(words)라면 Data Viewer에서 관심 단어 추적이 편해져요.
if isfield(Final, 'targetWordIdx') && ~isnan(Final.targetWordIdx(i))
    tIdx = Final.targetWordIdx(i);

    Eyelink('Message', 'TARGET_WORD_IDX %d', tIdx);

    if tIdx >= 1 && tIdx <= numel(words)
        Eyelink('Message', 'TARGET_WORD %d', tIdx);
        if (tIdx+1) <= numel(words)
            Eyelink('Message', 'POST_TARGET_WORD %d', tIdx+1);
        end
    end
end

% (선택) 표적/포스트표적 라벨 남기기 (실험 테이블에 인덱스가 있을 때)
% 예: Final.targetWordIdx(i) 가 1~numel(words) 라고 가정
% if isfield(Final,'targetWordIdx') && ~isnan(Final.targetWordIdx(i))
%     tIdx = Final.targetWordIdx(i);
%     Eyelink('Message','TARGET_WORD %d', tIdx);
%     if tIdx+1 <= numel(words)
%         Eyelink('Message','POST_TARGET_WORD %d', tIdx+1);
%     end
% end

% ... IAREA 메시지들 ...

% === 문장 + 오른쪽 ABC 실제 표시 (여기가 진짜 온셋) ===
Screen('FillRect', w, bgColor);
thisSentU8 = uint8(unicode2native(char(thisSent), 'UTF-8'));
DrawFormattedText(w, thisSentU8, TEXT_LEFT_X, 'center', bk, 0, [], [], 1);
draw_ABC_fixation(w, rightX, yMid, dp.ppd, ...
    'outer', outerDeg, 'cross', crossDeg, 'inner', innerDeg, ...
    'circlecolor', [0 0 0], 'crosscolor', bgColor);
[~, SearchOnset] = Screen('Flip', w);   % ← 앵커
% --- DUMMY 오른쪽 ABC도 호스트에 오버레이 표시 ---
if dummymode ~= 0 && Eyelink('IsConnected') == 1
    Eyelink('Command','clear_screen %d',7);
    hostDrawABC(rightX, yMid, dp.ppd, 'outer',1.5, 'cross',0.18, 'inner',0.18);
    Eyelink('Command','draw_box %d %d %d %d %d', tframe(1), tframe(2), tframe(3), tframe(4), 15);
end

% === Flip 직후에만 온셋 로그 ===
if Final.is_catch(i)
    Eyelink('Message','SENTENCE_ONSET_RIGHTABC CATCH');
    Eyelink('Message','CATCH_ONSET');
else
    Eyelink('Message','SENTENCE_ONSET_RIGHTABC %s_%s', ...
        char(Final.freq(i)), char(Final.valence(i)));
end

% === MAIN: 오른쪽 ABC dwell 게이트 (필수) ===
% (여기에 네가 이미 넣은 메인 오른쪽 dwell while-루프가 옴)
% ... RIGHT_FIX_CHECK / 시선 판정 루프(실모드) 또는 마우스 판정(더미) ...
% === MAIN: 오른쪽 ABC dwell 게이트 ===
phase_log('BEGIN','RIGHT_DWELL','trial',i);

% === 키 응답 처리 (1 / 2 / ESC / R=Recal) ===
% (여기부터 기존 응답 while 루프 시작)
holdSec  = 1.0;
winPx    = 75;
timeoutS = 10;

% === 커서 표시 정책: 실모드=숨김, 더미모드=표시 ===
if dummymode ~= 0
    ShowCursor('Arrow', w);   % dummy 모드일 때만 커서
else
    HideCursor;               % real 모드에서는 계속 숨김 유지
end

kb.only({'escape','r'});                 
kb.debounce();                           
t0 = GetSecs; tEnter = NaN;
dwellPauseOffsetSec = 0;

if dummymode == 0
    Eyelink('Message','RIGHT_FIX_CHECK');
    Eyelink('Command','clear_screen %d',7);
    hostDrawABC(rightX, yMid, dp.ppd, 'outer',1.5, 'cross',0.18, 'inner',0.18);
    Eyelink('Command','draw_box %d %d %d %d %d', tframe(1), tframe(2), tframe(3), tframe(4), 15);

    eye_used = Eyelink('EyeAvailable'); 
    if eye_used == 2, eye_used = 0; end
    if eye_used < 0, eye_used = 0; end

    try KbQueueStop(-1);  catch, end
    try KbQueueFlush(-1); catch, end
    RestrictKeysForKbCheck([]); 
    DisableKeysForKbCheck([]);

      while true
        % === 링크 상태 체크: 끊기면 EDF 정리 후 종료 ===
        checkLinkOrAbort('MAIN_RIGHT_DWELL', sprintf('trial %d RIGHT_DWELL', i), ...
            dirname, edfFile, dummymode, w, bgColor, pickedFont, pickedRenderer);

        % --- ESC / 재캘 ---
        if ~Screen('WindowKind', w)
        [w, rect, cx, cy] = ensureWindow(w, whichscreen, bgColor, 40, pickedFont, pickedRenderer);
        Screen('FillRect', w, bgColor); Screen('Flip', w);
    else
        rect = Screen('Rect', w);
    end

        [down,~,kc] = KbCheck(-3);
        if down
            if kc(keyEsc)
                kb.clear();
                Eyelink('SetOfflineMode'); WaitSecs(0.5);
                Eyelink('CloseFile'); filetransfer(dirname, edfFile, dummymode);
                phase_log('END','RIGHT_DWELL', 'trial', i, 'abort', 'ESC');
                cleanup; return;
            elseif kc(recal)
                Eyelink('Message','RIGHT_FIX_MANUAL_RECAL');
                tPause = GetSecs;
                try Eyelink('StopRecording'); end
                Eyelink('SetOfflineMode'); WaitSecs(0.1);
                EyelinkDoTrackerSetup(el);
                Eyelink('StartRecording'); WaitSecs(0.1);
                eye_used = Eyelink('EyeAvailable'); 
                if eye_used == 2, eye_used = 0; end
                if eye_used < 0, eye_used = 0; end

    % 참가자 화면 복구
                Screen('FillRect', w, bgColor);
                DrawFormattedText(w, thisSentU8, TEXT_LEFT_X, 'center', bk, 0, [], [], 1);
                draw_ABC_fixation(w, rightX, yMid, dp.ppd, ...
                    'outer', outerDeg, 'cross', crossDeg, 'inner', innerDeg, ...
                    'circlecolor', [0 0 0], 'crosscolor', bgColor);
                Screen('Flip', w);

                % 호스트 오버레이 복구
                Eyelink('Command','clear_screen %d',7);
                hostDrawABC(rightX, yMid, dp.ppd, 'outer',1.5, 'cross',0.18, 'inner',0.18);
                Eyelink('Command','draw_box %d %d %d %d %d', tframe(1), tframe(2), tframe(3), tframe(4), 15);

                % 보정/재설정
                dwellPauseOffsetSec = dwellPauseOffsetSec + (GetSecs - tPause);
                kb.only({'escape','r'}); kb.debounce();
                tEnter = NaN; t0 = GetSecs;
                continue;
            end
        end

        % --- 시선 샘플 → dwell 판정 ---
        if Eyelink('NewFloatSampleAvailable') > 0
            evt = Eyelink('NewestFloatSample');
            x = evt.gx(eye_used+1); y = evt.gy(eye_used+1);
            if ~isnan(x) && ~isnan(y) && x~=el.MISSING_DATA && y~=el.MISSING_DATA
                inWin = abs(x-rightX)<=winPx && abs(y-yMid)<=winPx;
                if inWin
                    if isnan(tEnter), tEnter = GetSecs; end
                    if (GetSecs - tEnter) >= holdSec
                        rightDwellSec = max(0, (GetSecs - tEnter) - dwellPauseOffsetSec);
                        Eyelink('Message','RIGHT_FIX_OK');
                        Eyelink('Message','!V TRIAL_VAR RightDwellMs %d', round(rightDwellSec*1000));
                        results.rightDwell(i,1) = rightDwellSec;
                        break;
                    end
                else
                    tEnter = NaN;
                end
            end
        end

         % --- 타임아웃 → 단순 리셋(재캘 없음) ---
         if (GetSecs - t0) > timeoutS
             Eyelink('Message','RIGHT_FIX_TIMEOUT');

             % 호스트/참가자 화면만 복구
             Eyelink('Command','clear_screen %d',7);
             hostDrawABC(rightX, yMid, dp.ppd, 'outer',1.5, 'cross',0.18, 'inner',0.18);
             Eyelink('Command','draw_box %d %d %d %d %d', ...
                 tframe(1), tframe(2), tframe(3), tframe(4), 15);

             Screen('FillRect', w, bgColor);
             DrawFormattedText(w, thisSentU8, TEXT_LEFT_X, 'center', bk, 0, [], [], 1);
             draw_ABC_fixation(w, rightX, yMid, dp.ppd, ...
                 'outer', outerDeg, 'cross', crossDeg, 'inner', innerDeg, ...
                 'circlecolor', [0 0 0], 'crosscolor', bgColor);
             Screen('Flip', w);

             % dwell 재시작
             t0     = GetSecs; 
             tEnter = NaN;
             continue;
         end
        WaitSecs(0.005);
    end
     kb.clear();  % ← 루프 정상 종료 후 한 번만
else
    
 % === 더미모드: 마우스 “응시” 판정 ===
    holdSec = 1.0; winPx = 75; timeoutS = 10;
    t0 = GetSecs; tEnter = NaN;
    while true
        [down,~,kc] = KbCheck(-3);
        if down
            if kc(keyEsc)
                RestrictKeysForKbCheck([]); phase_log('END','RIGHT_DWELL', 'trial', i, 'abort', 'ESC'); cleanup; return;
            elseif kc(recal)
                Screen('FillRect', w, bgColor);
                DrawFormattedText(w, thisSentU8, TEXT_LEFT_X, 'center', bk, 0, [], [], 1);
                draw_ABC_fixation(w, rightX, yMid, dp.ppd, ...
                    'outer', outerDeg, 'cross', crossDeg, 'inner', innerDeg, ...
                    'circlecolor', [0 0 0], 'crosscolor', bgColor);
                Screen('Flip', w);
                t0 = GetSecs; tEnter = NaN; safeKbReleaseWait; continue;
            end
        end

% [ANCHOR M-RIGHT-DUMMY] GetMouse 직전
    % === (가드) PTB 윈도우 핸들 복구 ===
    if ~Screen('WindowKind', w)
        fprintf('[SAFE] PTB window handle lost before GetMouse — reopening...\n');
        [w, rect, cx, cy] = ensureWindow(w, whichscreen, bgColor, 40, pickedFont, pickedRenderer);
    end  

        % === (기존) 마우스 기반 “응시” 판정 ===
        [mx,my] = GetMouse;
        inWin = ~isnan(mx) && ~isnan(my) && abs(mx-rightX)<=winPx && abs(my-yMid)<=winPx;
        if inWin
            if isnan(tEnter), tEnter = GetSecs; end
            if (GetSecs - tEnter) >= holdSec
                rightDwellSec = (GetSecs - tEnter);
                Eyelink('Message','RIGHT_FIX_OK DUMMY');
                Eyelink('Message','!V TRIAL_VAR RightDwellMs %d', round(rightDwellSec*1000));
                results.rightDwell(i,1) = rightDwellSec;
                break;
            end
        else
            tEnter = NaN;
        end

        if (GetSecs - t0) > timeoutS
            Screen('FillRect', w, bgColor);
            DrawFormattedText(w, thisSentU8, TEXT_LEFT_X, 'center', bk, 0, [], [], 1);
            draw_ABC_fixation(w, rightX, yMid, dp.ppd, ...
                'outer', outerDeg, 'cross', crossDeg, 'inner', innerDeg, ...
                'circlecolor', [0 0 0], 'crosscolor', bgColor);
            Screen('Flip', w);
            t0 = GetSecs; tEnter = NaN;
        end
        WaitSecs(0.005);
    end
end

phase_log('END','RIGHT_DWELL','ok',i);


% (한 번만 정의돼 있으면 생략 가능. 안전하게 여기서도 보강)
keyA    = KbName('1');
keyL    = KbName('2');
keyEsc  = KbName('ESCAPE');
recal   = KbName('r');     % ← 재캘 키

% === PROMPT 화면 띄우고, 여기부터 1/2/ESC/R 받기 (Main) ===
Screen('FillRect', w, bgColor);
Screen('TextSize', w, 40);
promptTxt = ['방금 제시된 문장이 의미적으로 자연스러웠나요?\n' ...
             '자연스러웠다: 1, 어색했다: 2'];
DrawFormattedText(w, u8(promptTxt), 'center','center', bk, 70, [], [], 1.4);
[~, PromptOnset] = Screen('Flip', w);
Eyelink('Message','PROMPT_ONSET');

% 프롬프트 직전 초기화/관심키 설정 (queues 전용)
kbPad.clear();  kbKb.clear();
kbPad.debounce(); kbKb.debounce();
kbPad.only([RESP1_KEYS, RESP2_KEYS]);                 % keypad: 1/2
kbKb.only([RESP1_KEYS, RESP2_KEYS, KBD_ADMIN]);       % keyboard: 1/2 + ESC/R
pauseOffsetSec = 0;
                        
name = ""; 
tResp = NaN;

while name==""
    % keypad: 1/2만
    [nP, tP] = kbPad.wait([RESP1_KEYS, RESP2_KEYS], 0.02);
    % keyboard: 1/2 + ESC/R
    [nK, tK] = kbKb.wait([RESP1_KEYS, RESP2_KEYS, KBD_ADMIN], 0.02);

    % 1) ESC: 즉시 종료(키보드)
    if (nK~="") && strcmpi(nK,'escape')
        if dummymode==0
            Eyelink('SetOfflineMode'); WaitSecs(0.5);
            Eyelink('CloseFile'); filetransfer(dirname, edfFile, dummymode);
        end
        RestrictKeysForKbCheck([]); cleanup; return;
    end
     
      % 2) 재캘(R): 장비 재설정 → 화면 복구 → 시간 보정 → 다시 대기
    if (nK~="") && strcmpi(nK,'r') && dummymode==0
        tPauseStart = GetSecs;
        Eyelink('Message','!!! RECALIBRATION !!!');
        try Eyelink('StopRecording'); end
        Eyelink('SetOfflineMode');  WaitSecs(0.1);
        EyelinkDoTrackerSetup(el);
        Eyelink('StartRecording');  WaitSecs(0.1);

        % (선택) 호스트 오버레이 복구
        Eyelink('Command','clear_screen %d',7);
        Eyelink('Command','draw_cross %d %d 15', cx, cy);
        Eyelink('Command','draw_box %d %d %d %d %d', tframe(1), tframe(2), tframe(3), tframe(4), 15);

        % 참가자 화면 복구: 문장 + 오른쪽 ABC
        Screen('FillRect', w, bgColor);
        DrawFormattedText(w, thisSentU8, TEXT_LEFT_X, 'center', bk, 0, [], [], 1);
        draw_ABC_fixation(w, rightX, yMid, dp.ppd, ...
            'outer', outerDeg, 'cross', crossDeg, 'inner', innerDeg, ...
            'circlecolor', [0 0 0], 'crosscolor', bgColor);
        Screen('Flip', w);

         pauseOffsetSec = pauseOffsetSec + (GetSecs - tPauseStart);

        % 큐 재설정 후 재대기
        kbPad.clear();  kbKb.clear();
        kbPad.debounce(); kbKb.debounce();
        kbPad.only([RESP1_KEYS, RESP2_KEYS]);
        kbKb.only([RESP1_KEYS, RESP2_KEYS, KBD_ADMIN]);
        continue;
    end

  % 3) 응답키(키패드 우선, 없으면 키보드)
    if name==""
        if nP~=""
            name  = nP;  tResp = tP;   % keypad: 1/2
        elseif (nK~="") && ~any(strcmpi(char(nK), KBD_ADMIN))  % ESC/R 제외
            name  = nK;  tResp = tK;   % keyboard: 1/2
        end
    end  

% === 응답 시점 디버그 로그 + RT 계산 (루프 바깥에서 딱 한 번) ===
if ~isnan(tResp)
    RT = (tResp - PromptOnset) - pauseOffsetSec;
else
    RT = NaN;
end

fprintf('[DBG] END TRIAL %d: resp=%s | RT=%.3f s @%.3f\n', ...
        i, char(name), RT, tResp);

% === 4) 1 / 2 응답 처리 (기존 코드 계속) ===
n = char(string(name));

if exist('RESP1_KEYS','var') && ~isempty(RESP1_KEYS)
    is1 = any(strcmpi(n, RESP1_KEYS));
else
    is1 = strcmpi(n,'1');
end
if exist('RESP2_KEYS','var') && ~isempty(RESP2_KEYS)
    is2 = any(strcmpi(n, RESP2_KEYS));
else
    is2 = strcmpi(n,'2');
end

results.RT(i,1) = RT;

if is1 || is2
    % 시간 소스: tResp → GetSecs 폴백 (본실험은 PromptOnset 사용)
    if exist('tResp','var') && ~isempty(tResp) && isfinite(tResp)
        tUse = tResp;
    else
        tUse = GetSecs;
    end

    rtMs    = round( (tUse - PromptOnset - pauseOffsetSec) * 1000 );
    isCatch = Final.is_catch(i) ~= 0;

    if is1
        resp = '1'; acc = ~isCatch;
    else
        resp = '2'; acc =  isCatch;
    end

        results.RT(i,1)      = rtMs/1000;
        results.resp(i,1)    = string(resp);
        results.acc(i,1)     = acc;
        results.isCatch(i,1) = isCatch;

        Eyelink('Message', sprintf('RESPONSE %s %d ACC=%d IS_CATCH=%d', ...
            char(resp), int32(rtMs), int32(acc), int32(isCatch)));

        % 피드백
        if acc
            fbText='맞았습니다'; fbColor=[0 200 0]; fbTag='CORRECT';
        else
            fbText='틀렸습니다'; fbColor=[230 0 0]; fbTag='INCORRECT';
        end
        Eyelink('Message', ['FEEDBACK_ON ' fbTag]);
        Screen('FillRect', w, bgColor);
        Screen('TextSize', w, 50);
        DrawFormattedText(w, u8(fbText), 'center','center', fbColor, 0, [], [], 1.2);
        Screen('Flip', w);
        WaitSecs(1.0);
        Screen('FillRect', w, bgColor);
Screen('Flip', w);
        Eyelink('Message', 'FEEDBACK_OFF');

        safeKbReleaseWait;       % 중복 입력 방지
        break;               % ← 루프 탈출 (여기서 끝!)
    end
end

RestrictKeysForKbCheck([]);  % 해제

    % Stop recording eye movements at the end of each trial
    WaitSecs(0.1); % Add 100 msec of data to catch final events before stopping, 마지막 이벤트/샘플이 버퍼에 다 들어오도록 100ms 대기
    Eyelink('StopRecording'); % 이번 트라이얼의 시선 기록 종료

    % Add variables to EDF
    Eyelink('Message', '!V TRIAL_VAR TargetLocation %d %d', round(results.target(i, 1)), round(results.target(i, 2))); % TargetLocation → 타깃의 (x,y) 픽셀 좌표
    Eyelink('Message', '!V TRIAL_VAR ReactionTime %d', fix(results.RT(i)*1000)); % ms, round toward zero; ReactionTime → 반응시간을 ms로 변환. fix는 소수점 버림

    % Write TRIAL_RESULT message to EDF file: marks the end of a trial for DataViewer (paired with TRIALID)
    Eyelink('Message', 'TRIAL_RESULT 0'); % TRIAL_RESULT 0 = 트라이얼 정상 종료
    WaitSecs(0.01); % Allow some time before ending the trial

% 체크포인트 저장 O A(10 trial마다)
if mod(i,10)==0
    try
        save(fullfile(dirname, sprintf('checkpoint_%03d.mat', i)), 'results', 'dp', 'Final', 'i');

    end
end

        % --- 블록 사이 휴식(15초) ---
    if i < numTrial
        curBlock  = Final.block(i);
        nextBlock = Final.block(i+1);

        % 블록 경계에만 휴식 적용
        if nextBlock ~= curBlock
            % 로그
            Eyelink('Message','BLOCK_END %d', curBlock);
            Eyelink('Message','BREAK_START %d', 15);

            % 1) "쉬는 시간입니다" 화면 (15초) — ESC 허용
            Screen('FillRect', w, bgColor);
            DrawFormattedText(w, u8('쉬는 시간입니다.\n\n휴식 시간: 15초'), ...
                  'center','center', bk, 60, [], [], 1.4);
            Screen('Flip', w);

            kb.clear();          % 잔여 입력/큐 정리
            kb.debounce();       % 눌린 키가 남았다면 떼질 때까지 대기
            kb.only({'escape'});   % 휴식 중엔 ESC만 허용

            tEnd = GetSecs + 15;
            while GetSecs < tEnd
                [name] = kb.wait('escape', 0.05);  % 50 ms 타임슬라이스, non-blocking
                if strcmpi(name,'escape')
                    kb.clear();
                    if dummymode==0
                        Eyelink('SetOfflineMode'); WaitSecs(0.5);
                        Eyelink('CloseFile'); filetransfer(dirname, edfFile, dummymode);
                    end
                    cleanup; return;
                end
            end
            kb.clear();
            Eyelink('Message','BREAK_OVER');

            % 2) "Enter 버튼을 누르면 다음 블록 시작" — Enter/ESC 모두 허용
            Screen('FillRect', w, bgColor);
            DrawFormattedText(w, u8('15초가 지났습니다.\nEnter 버튼을 누르면 다음 블록을 시작합니다.'), ...
                              'center','center', bk, 60, [], [], 1.4);
            Screen('Flip', w);

            kb.debounce();                           % 눌려 있던 키 해제
kb.only({'space','return','escape'});     % enter/ESC만 허용
[name] = kb.wait({'space','return','escape'}, inf); % 블로킹 대기
kb.clear();                              % 제한/큐 해제

if strcmpi(name,'escape')
    if dummymode==0
        Eyelink('SetOfflineMode'); WaitSecs(0.5);
        Eyelink('CloseFile'); filetransfer(dirname, edfFile, dummymode);
    end
    cleanup; return;
end

% 여기까지 왔으면 enter → 다음 블록으로
Eyelink('Message','BLOCK_START %d', nextBlock);
Screen('FillRect', w, bgColor); Screen('Flip', w);  % (선택) 화면 정리
        end % <-- close 'if nextBlock ~= curBlock'
    end % <-- close 'if i < numTrial'

% --- 마지막 trial 처리 ---
if i == numTrial
    % 종료 메시지
    [w, ~, ~, ~] = ensureWindow(w, whichscreen, bgColor, 50, pickedFont, pickedRenderer);
    Screen('TextStyle', w, 0);
    msg = u8('실험이 모두 종료되었습니다. 실험자를 불러주세요.');
    DrawFormattedText(w, msg, 'center', 'center', [0 0 0]);
    Screen('Flip', w);

    % --- 자동 저장 (Flip 직후 즉시) ---
    try
        save(fullfile(dirname, 'EL_demo.mat'), 'results', 'dp', 'Final');
        fprintf('[SAVE] EL_demo.mat written to %s (auto after final flip)\n', dirname);
    catch ME_SAVE
        warning('AUTO SAVE ERROR → %s', getReport(ME_SAVE,'basic'));
    end

    % --- (실모드일 때만) EDF 닫고 수신 ---
    if exist('dummymode','var') && dummymode==0 && Eyelink('IsConnected')==1
        Eyelink('SetOfflineMode'); WaitSecs(0.3);
        try, Eyelink('CloseFile'); catch, end
        WaitSecs(0.3);
        try
            filetransfer(dirname, edfFile, dummymode);  % 로컬 함수 사용
        catch errEDF
            warning('EDF transfer failed: %s', getReport(errEDF,'basic'));
        end
    end

    DisableKeysForKbCheck([]);   % ← 이미 넣은 줄
try KbQueueRelease(-1); end  % ← 선택: 잔여 큐 있으면 해제
    
kb.debounce();            % 혹시 눌린 키 있으면 다 떼질 때까지
kb.only({'escape'});
kb.wait('escape', inf);
kb.clear();
cleanup;
return;                                     % 함수 종료
end
end % <-- close 'for i = 1:numTrial'

catch ME     % 비정상 종료라도 중간까지의 결과는 남기기 (가능할 때만)     
    try         
        if exist('results','var') && exist('dp','var') && exist('dirname','var')             
            save(fullfile(dirname,'EL_demo_crash.mat'), 'results', 'dp');         end
    end
    cleanup;
    rethrow(ME);
end   % try–catch 종료
end  % visualsearch_EL_demo  % ← 메인 함수 닫기(필수)
    
%% ===== Local functions =====
function cleanup
    % PTB/EL 상태를 가능한 한 조용히 원복
    try Screen('CloseAll'); end
    try Eyelink('Shutdown'); end
    try ShowCursor; end
end

function fix = infixationWindow(mx, my, rx, ry)
    % check if eye position falls within 75 pixel window 
    % (rx, ry): center of the ROI, (rx, ry) 중심으로 가로~세로 ±75px 사각 창 안에 들어오면 고정으로 간주
    fix = abs(rx - mx) < 75 && abs(ry - my) < 75; 
end

function filetransfer(dirname, edfFile, dummymode)
    if dummymode == 1
        fprintf('No EDF file saved in Dummy mode \n');  
        return;  % 더미 모드면 실제 EDF 없음 → 바로 종료
    end

    try
        fprintf('Receiving data file ''%s''\n', edfFile);
    
        status = Eyelink('ReceiveFile', edfFile, dirname, 1); % status: 전송된 바이트 수(>0이면 성공적으로 받음)
        
        if status > 0
            fprintf('\tEDF file size: %.1f KB\n', status/1024); % Divide file size by 1024 to convert bytes to KB
        end

        if 2==exist(edfFile, 'file')
            fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, dirname);
        end
    
    catch
        fprintf('Problem receiving data file ''%s''\n', edfFile);
    end
end

%% supporting functions for experiment
function draw_fixation(w, cx, cy, col)
    Screen('DrawLine', w, col, cx-20, cy, cx+20, cy, 5); % 가로 줄(길이 40px, 두께 5)
    Screen('DrawLine', w, col, cx, cy-20, cx, cy+20, 5); % 세로 줄
end

function [grid] = make_grid(rect, ss) % input = dp.resolution

xsize = rect(1);
ysize = rect(2);

slotsize = (ss*2)*1.5; % adjust the value depending on the display resolution, 자극 지름(2*ss)의 1.5배 간격 → 서로 겹치지 않게

numx = floor(xsize / slotsize) - 2; % 가장자리 여백을 위해 2칸 감산
numy = floor(ysize / slotsize) - 2;

xmargin = (xsize - (numx * slotsize)) / 2; % 좌우 중앙 정렬
ymargin = (ysize - (numy * slotsize)) / 2; % 상하 중앙 정렬

xx = (0 : slotsize : (slotsize * (numx - 1))) + (slotsize/2) + xmargin; % centres, 각 칸의 x-중심
yy = (0 : slotsize : (slotsize * (numy - 1))) + (slotsize/2) + ymargin; % 각 칸의 y-중심

ytemp = zeros(1, numx*numy);   % 최종 길이만큼 미리 공간 확보
for i = 1:numy
    idx = (i-1)*numx + (1:numx);
    ytemp(idx) = repmat(yy(i), 1, numx);
end

grid = [repmat(xx, 1, numy); ytemp]; % all possible combinations of x, y coordiates, 2×(numx*numy): 모든 (x,y) 조합

end

function out = cleanSentence(inStr)
    if ismissing(inStr), s = ''; else, s = char(inStr); end

    % 제어/DEL/NBSP → 공백
    s = regexprep(s, '[\x00-\x1F\x7F\xA0]', ' ');
    % 제로폭/BOM 제거
    s = regexprep(s, '[\u200B-\u200F\uFEFF]', '');
    % 서로게이트/PUA/비문자 제거
    s = regexprep(s, '[\uD800-\uDFFF\uE000-\uF8FF\uFDD0-\uFDEF\uFFFE\uFFFF]', '');

    % ← 정규식 화이트리스트 '삭제'
    % 유니코드 속성(문자/숫자/공백/구두점)만 유지
    keep = isstrprop(s,'alpha') | isstrprop(s,'digit') | isstrprop(s,'wspace') | isstrprop(s,'punct');
    s(~keep) = ' ';

    s = regexprep(s, '\s+', ' ');
    out = string(strtrim(s));
end

function [w, rect, cx, cy] = ensureWindow(w, whichscreen, bgColor, fontsize, pickedFont, pickedRenderer)
    if nargin < 4 || isempty(fontsize), fontsize = 40; end

    if ~exist('w','var') || isempty(w) || Screen('WindowKind', w) <= 0
        [w, rect] = Screen('OpenWindow', whichscreen, bgColor);
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    else
        rect = Screen('Rect', w);
    end

    [cx, cy] = RectCenter(rect);

    % === 텍스트 상태: 선택된 조합만 사용 ===
    Screen('Preference','TextEncodingLocale','UTF-8');
    if ~isnan(pickedRenderer)
        Screen('Preference','TextRenderer', pickedRenderer);
    end
    Screen('TextStyle', w, 0);
    Screen('TextFont',  w, pickedFont);
    Screen('TextSize',  w, fontsize);
end

function draw_ABC_fixation(w, cx, cy, ppd, varargin)
% ABC fixation (bulls-eye + crosshair)
% Default sizes (deg): outer 0.6°, inner 0.2°, cross thickness 0.2°
% Colors: circle=black, cross=white

p = inputParser;
addParameter(p, 'outer', 0.6);                 % deg (diameter)
addParameter(p, 'inner', 0.2);                 % deg (diameter)
addParameter(p, 'cross', []);                  % deg (thickness); default = inner
addParameter(p, 'circlecolor', [0 0 0]);       % RGB
addParameter(p, 'crosscolor',  [255 255 255]); % RGB
parse(p, varargin{:}); 

outerDeg = p.Results.outer;
innerDeg = p.Results.inner;
crossDeg = p.Results.cross;
if isempty(crossDeg), crossDeg = innerDeg; end

% 픽셀 단위 계산 (안정성 위해 정수로 반올림)
ro = round((outerDeg * ppd) / 2);   % outer radius (px)
ri = round((innerDeg * ppd) / 2);   % inner radius (px)
lw = max(1, round(crossDeg * ppd)); % cross thickness (px, >=1)

% 좌표도 정수로
cx = round(cx); cy = round(cy);

% 1) 바깥 원(검정)
Screen('FillOval', w, p.Results.circlecolor, [cx-ro, cy-ro, cx+ro, cy+ro]);

% 2) 십자선(사각형 두 개로)
half = floor(lw/2);
hRect = [cx-ro, cy-half, cx+ro, cy+half]; % 가로 막대
vRect = [cx-half, cy-ro, cx+half, cy+ro]; % 세로 막대
Screen('FillRect', w, p.Results.crosscolor, hRect);
Screen('FillRect', w, p.Results.crosscolor, vRect);

% 3) 안쪽 원(검정)
Screen('FillOval', w, p.Results.circlecolor, [cx-ri, cy-ri, cx+ri, cy+ri]);
end

function [wordRects, words, lineRect] = makeWordROIsSingleLine(w, cx, cy, sentChar, leftX)
% 단일 줄 문장(sentChar)을 화면 중심(cx,cy)에 표시한다고 가정하고
% 각 단어의 [L T R B] 박스, 단어 리스트, 전체 라인 박스를 반환.
%
% DrawFormattedText와 동일한 폰트/렌더러/사이즈가 이미 설정돼 있어야 함!
% (호출 전: Screen('Preference','TextRenderer', pickedRenderer);
%          Screen('TextFont', w, pickedFont);
%          Screen('TextSize', w, <그때 쓰는 값>))

    % UTF-8 헬퍼
    u8 = @(s) uint8(unicode2native(s,'UTF-8'));
    renderer = Screen('Preference','TextRenderer');
    useUTF8  = (renderer == 1);

    % 0) 전체 문장 bounds 측정
    if useUTF8
        tbAll = Screen('TextBounds', w, u8(sentChar));
    else
        tbAll = Screen('TextBounds', w, sentChar);
    end
    fullW = tbAll(3) - tbAll(1);
    fullH = tbAll(4) - tbAll(2);

    % === [ANCHOR ROI-FALLBACK] leftX 기본값 (하위호환) ===
    if nargin < 5 || isempty(leftX)
        leftX = cx - fullW/2;   % 과거(center) 동작과 동일
    end

    % DrawFormattedText의 'center','center'와 동일하게 좌상단 기준 계산
    startX = round(leftX);      % ← 호출부에서 TEXT_LEFT_X를 넘김
    startY = round(cy - fullH/2);
    lineRect = [startX, startY, startX + fullW, startY + fullH];

    % 1) 단어 분할 (공백 기준; 멀티 스페이스는 하나씩 취급)
    words = split(string(sentChar));  % string형 배열
    words = words(:).';               % 1×N

    % 2) 공백 폭 측정
    if useUTF8
        tbSp = Screen('TextBounds', w, u8(' '));
    else
        tbSp = Screen('TextBounds', w, ' ');
    end
    spaceW = tbSp(3) - tbSp(1);

    % 3) 각 단어 폭 측정 → 순차 배치
    N = numel(words);
    wordRects = zeros(N, 4);
    x = startX;
    for k = 1:N
        tok = char(words(k));
        if useUTF8
            tb = Screen('TextBounds', w, u8(tok));
        else
            tb = Screen('TextBounds', w, tok);
        end
        wtok = tb(3) - tb(1);

        L = round(x);
        T = startY;
        R = round(x + wtok);
        B = startY + fullH;

        wordRects(k, :) = [L T R B];

        % 다음 단어 시작점(마지막 단어 뒤에는 공백 추가 없음)
        x = x + wtok;
        if k < N
            x = x + spaceW;
        end
    end
    wordRects = round(wordRects);     % 좌표 정수화
    lineRect  = round(lineRect);      % (권장) 라인 박스도 정수화
end

function safeKbReleaseWait(deviceNumber)
    if nargin < 1 || isempty(deviceNumber) || deviceNumber < 0
        deviceNumber = -1;  % all devices
    end
    try
        % wait until all keys are released and PTB window is in focus
        while KbCheck(deviceNumber)
            WaitSecs(0.05);
        end

        % explicit short idle to stabilize focus handover
        WaitSecs(0.1);

        % instead of KbReleaseWait, use a manual release-wait loop
        while true
            [down,~,~] = KbCheck(deviceNumber);
            if ~down, break; end
            WaitSecs(0.01);
        end

    catch ME
        if contains(ME.message, 'AbortedByUser') || ...
           contains(ME.message, 'terminated by user')
            fprintf('[SAFE] KbReleaseWait aborted by user/focus change — ignored.\n');
        else
            rethrow(ME);
        end
    end
end

function [mx, my, w, rect, cx, cy] = safeGetMouse(w, whichscreen, bgColor, pickedFont, pickedRenderer, fontsize)
    if nargin < 6 || isempty(fontsize), fontsize = 40; end

    % 창 유효성 점검 & 필요 시 복구
    if ~Screen('WindowKind', w)
        fprintf('[SAFE] PTB window handle lost. Reopening...\n');
        [w, rect, cx, cy] = ensureWindow(w, whichscreen, bgColor, fontsize, pickedFont, pickedRenderer);
        Screen('FillRect', w, bgColor);
        Screen('Flip', w);
        WaitSecs(0.02);
    else
        rect = Screen('Rect', w);
        [cx, cy] = RectCenter(rect);
    end

    % === 여기! 윈도우 핸들 대신 스크린 번호 사용 ===
    try
        [mx, my] = GetMouse(whichscreen);   % <-- 핵심 수정
    catch
        % 폴백: 인자 없이(기본 스크린) 시도
        [mx, my] = GetMouse;
    end
end

function hostDrawABC(cx, cy, ppd, varargin)
% 호스트(Operator Display)에 ABC(bulls-eye + crosshair) 오버레이를 그림
% 사용 예: hostDrawABC(leftX, yMid, dp.ppd, 'outer',1.5, 'cross',0.18, 'inner',0.18);

p = inputParser;
addParameter(p, 'outer', 1.5);   % deg (바깥 원 지름)
addParameter(p, 'inner', 0.18);  % deg (안쪽 점 지름)
addParameter(p, 'cross', 0.18);  % deg (십자 두께/사이즈 기준치)
parse(p, varargin{:});
outerDeg = p.Results.outer;
innerDeg = p.Results.inner;
crossDeg = p.Results.cross;

% 픽셀 환산(반지름/두께 등 정수화)
ro = round((outerDeg * ppd) / 2);           % outer radius
ri = round((innerDeg * ppd) / 2);           % inner radius
cw = max(1, round(crossDeg * ppd));         % cross half-length/굵기 기준치

% 색상/굵기 파라미터: EyeLink 명령의 마지막 인자는 색/굵기 인덱스용으로 15를 관례대로 사용
col = 15;

% 안전: 좌표 정수화
cx = round(cx); cy = round(cy);

% 바탕을 유지하고 덧그리기만: 필요 시 clear_screen으로 지우고 다시 그릴 것
% 바깥 원 (Eyelink는 draw_circle 미지원 → 사각 테두리로 대체)
Eyelink('Command','draw_box %d %d %d %d %d', cx-ro, cy-ro, cx+ro, cy+ro, col);

% 십자(가로/세로 선)
Eyelink('Command','draw_line %d %d %d %d %d', cx-ro, cy, cx+ro, cy, col);
Eyelink('Command','draw_line %d %d %d %d %d', cx, cy-ro, cx, cy+ro, col);

% 안쪽 원도 동일하게 박스로 대체
Eyelink('Command','draw_box %d %d %d %d %d', cx-ri, cy-ri, cx+ri, cy+ri, col);
end

function debounce_pollonly(dev)
    if nargin<1 || isempty(dev), dev = -1; end
    quiet = 0.12;            % 120 ms 무입력 구간
    t0 = GetSecs;
    while true
        [down,~,~] = KbCheck(dev);
        if ~down
            if (GetSecs - t0) >= quiet, break; end
        else
            t0 = GetSecs;    % 키 눌리면 타이머 리셋
        end
        WaitSecs(0.01);
    end
    WaitSecs(0.03);          % 짧은 안정화
end

function devIdx = pick_keypad_device(timeoutSec)
% 외장 키패드 자동 탐지:
%  1) 후보 장치(키보드류) 스캔
%  2) 각 장치에 임시 키큐 생성 후 동시에 스타트
%  3) "외장 키패드 아무 키" 눌러달라고 안내
%  4) 이벤트가 들어온 장치를 devIdx로 선택
%
% 실패 시 -1 반환(전체 장치), 이후 이름 기반 매칭으로라도 동작

if nargin<1 || isempty(timeoutSec), timeoutSec = 7; end

try
    devices = PsychHID('Devices');
catch
    fprintf('[KEYPAD] PsychHID Devices() 실패 → -1\n');
    devIdx = -1; return;
end

% --- 1) 키보드류 후보만 추리기 ---
cands = [];
for k = 1:numel(devices)
    d = devices(k);
    % PTB 표준: usage/usageName가 Keyboard 또는 Keypad 계열
    isKbd = false;
    if isfield(d,'usageName') && ~isempty(d.usageName)
        nm = lower(string(d.usageName));
        isKbd = contains(nm, "keyboard") || contains(nm, "keypad");
    elseif isfield(d,'usage') && ~isempty(d.usage)
        % HID Usage 6 = Keyboard
        isKbd = (d.usage == 6);
    end
    if isKbd
        cands(end+1) = k; %#ok<AGROW>
    end
end

if isempty(cands)
    fprintf('[KEYPAD] 키보드류 HID 후보가 없습니다 → -1\n');
    devIdx = -1; return;
end

% --- 2) 각 후보 장치에 임시 큐 생성 & 시작 ---
queues = nan(size(cands));
nKeys = 256;  % 대부분 환경에서 충분. (KbMaxNumberOfKeys 대체)
keyFlags = zeros(1, nKeys); % double
keyFlags(:) = 1;            % 모든 키 이벤트 허용(장치 구분 목적)

for i = 1:numel(cands)
    di = devices(cands(i)).index;
    try
        % 혹시 기존 큐가 잡고 있으면 해제 시도
        try KbQueueRelease(di); catch, end
        KbQueueCreate(di, keyFlags);
        KbQueueStart(di);
        queues(i) = di;
        fprintf('[KEYPAD] 후보 %d: "%s" (index=%d)\n', i, devices(cands(i)).product, di);
    catch ME
        fprintf('[KEYPAD] 큐 생성 실패: "%s" (idx=%d) → %s\n', devices(cands(i)).product, di, ME.message);
    end
end

fprintf('\n[KEYPAD] 외장 키패드의 **아무 키**를 눌러주세요 (%.1f초 이내)…\n', timeoutSec);

% --- 3) 동시 폴링으로 이벤트가 들어온 첫 장치 선택 ---
t0 = GetSecs;
hitIdx = [];
while (GetSecs - t0) < timeoutSec && isempty(hitIdx)
    for i = 1:numel(cands)
        di = queues(i);
        if ~isnan(di)
            try
                [pressed, firstPress] = KbQueueCheck(di);
            catch
                % 장치 제거/에러 → 스킵
                continue;
            end
            if pressed && any(firstPress>0)
                hitIdx = i; break;
            end
        end
    end
    WaitSecs(0.01);
end

% --- 4) 큐 정리 & 결과 반환 ---
for i = 1:numel(cands)
    di = queues(i);
    if ~isnan(di)
        try KbQueueStop(di); catch, end
        try KbQueueRelease(di); catch, end
    end
end

if ~isempty(hitIdx)
    picked = devices(cands(hitIdx));
    devIdx = picked.index;
    fprintf('[KEYPAD] 선택: "%s" (index=%d)\n\n', picked.product, devIdx);
else
    fprintf('[KEYPAD] 시간 내 입력 없음 → -1로 폴백\n\n');
    devIdx = -1;
end
end

function [expired] = watchdog(t0, hardSec, tag)
% hardSec 초가 지나면 true 반환하고, 화면/콘솔에 로그 남김
% tag: 현재 페이즈 이름(예: 'LEFT_DWELL', 'PROMPT_WAIT' 등)
    if GetSecs - t0 > hardSec
        expired = true;
        try
            fprintf('[WATCHDOG] %s hard-timeout (>%gs) → force break\n', tag, hardSec);
        end
    else
        expired = false;
    end
end


