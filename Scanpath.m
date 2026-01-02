%% TRIALID Scanpath + ROI + 단어라벨 + fixation(scatter) + start/end
%    + COORD모드표기 + word1 inROI QC
%    + SENTENCE_ONSET ~ PROMPT_ONSET(미포함) fixation만 보기
%    + prePROMPT sample: gx/gy가 정규화(0~1)면 픽셀로 변환해서 표시(아니면 px로 간주)
%
% NOTE: ABC(LEFT/RIGHT) + gate 관련 코드는 완전 제거

trialNumWanted = 66;

% ===== SENTENCE_ONSET ~ PROMPT_ONSET(미포함)만 보기 =====
onlySentenceBeforePrompt = true;   % true: SENTENCE_ONSET <= fix.onset < PROMPT_ONSET

% 색
scanBlue    = [0 0.4470 0.7410];
startYellow = [0.9290 0.6940 0.1250];
endOrange   = [0.8500 0.3250 0.0980];

% fixation scatter 옵션
dotSizeMode  = "dur";   % "fixed" 또는 "dur"
dotSizeFixed = 40;
dotSizeMin   = 20;
dotSizeMax   = 120;

%% 0) 기본 전제 체크
assert(exist("subj","var")==1 && isfield(subj,"trial") && isfield(subj,"event") && isfield(subj.event,"fix"), ...
    "subj / subj.trial / subj.event.fix 필요");
assert(exist("results","var")==1, "results 필요");

%% 1) TRIALID -> subj trial index (tSubj)
trialIDstr = sprintf("TRIALID %d", trialNumWanted);
tSubj = find(strcmp({subj.trial.id}, trialIDstr), 1);
assert(~isempty(tSubj), "%s not found in subj.trial", trialIDstr);

% ★ 핵심: results 인덱스는 tSubj 기준으로 접근
tw = tSubj;

%% 1-1) ROI(rects) / words(wlist): results{tw}
assert(isfield(results,'wordRects') && tw>=1 && tw<=numel(results.wordRects) && ~isempty(results.wordRects{tw}), ...
    "results.wordRects{%d} missing/empty (tw=tSubj)", tw);

rects  = results.wordRects{tw};   % [L T R B]
nWords = size(rects,1);

wlist = string.empty(0,1);
if isfield(results,'words') && tw<=numel(results.words) && ~isempty(results.words{tw})
    wlist = string(results.words{tw});
end

if ~isempty(wlist)
    fprintf('[SENT] %s | tw=%d | %s\n', trialIDstr, tw, strjoin(wlist," "));
else
    fprintf('[SENT] %s | tw=%d | (wlist empty)\n', trialIDstr, tw);
end

%% 2) fixations (subj.trial(tSubj)에서 추출)
useDurClean = false;
if isfield(subj.trial(tSubj),'fixIdxDurClean') && ~isempty(subj.trial(tSubj).fixIdxDurClean)
    fixIdx = subj.trial(tSubj).fixIdxDurClean(:);
    useDurClean = true;
elseif isfield(subj.trial(tSubj),'fixIdx') && ~isempty(subj.trial(tSubj).fixIdx)
    fixIdx = subj.trial(tSubj).fixIdx(:);
else
    error("No fixation indices for %s.", trialIDstr);
end

fixIdx = fixIdx(fixIdx>=1 & fixIdx<=numel(subj.event.fix));
fix = subj.event.fix(fixIdx);

% validity 필터(있으면)
if isfield(fix,'isValid')
    keep = [fix.isValid]';
    fix    = fix(keep);
    fixIdx = fixIdx(keep);
end
assert(~isempty(fix), "No valid fixations left for %s.", trialIDstr);

% 시간순 정렬
if isfield(fix,'onset')
    [~, srt] = sort([fix.onset]');
    fix    = fix(srt);
    fixIdx = fixIdx(srt);
end

% duration(점 크기)
if isfield(fix,'dur')
    dur = [fix.dur]';
else
    dur = nan(numel(fix),1);
end

srcName = "fixIdx";
if useDurClean, srcName = "fixIdxDurClean"; end
fprintf('[FIXSRC] %s (n=%d)\n', srcName, numel(fix));

%% 2-1) Raw vs Corr 모드 선택 + dx/dy(median) 계산 (필터 전에 dx/dy 산출용)
hasCorrFields = isfield(fix,'xCorr') && isfield(fix,'yCorr') && isfield(fix,'x') && isfield(fix,'y');

useCorr = false;
if hasCorrFields
    useCorr = all(isfinite([fix.xCorr]')) && all(isfinite([fix.yCorr]'));
end

if useCorr
    coordMode = "CORR";
    dxMed_all = median([fix.xCorr]' - [fix.x]', 'omitnan');
    dyMed_all = median([fix.yCorr]' - [fix.y]', 'omitnan');
else
    coordMode = "RAW";
    dxMed_all = NaN; dyMed_all = NaN;
end

%% ===== SENTENCE_ONSET ~ PROMPT_ONSET(미포함) fixation만 남기기 =====
t = tSubj;
tOn = NaN; tPrompt = NaN;

if onlySentenceBeforePrompt
    assert(isfield(subj.trial(t),'sentenceOnset') && isfinite(subj.trial(t).sentenceOnset), ...
        'sentenceOnset missing/NaN for this trial');
    tOn = subj.trial(t).sentenceOnset;

    msgIdxT = subj.trial(t).msgIdx(:);
    txtM = string({subj.event.msg(msgIdxT).text})';
    timM = [subj.event.msg(msgIdxT).time]';
    pIdx = find(contains(txtM,"PROMPT_ONSET",'IgnoreCase',true), 1, 'first');
    assert(~isempty(pIdx), 'PROMPT_ONSET not found in msgs');
    tPrompt = timM(pIdx);

    inWin = [fix.onset]' >= tOn & [fix.onset]' < tPrompt;

    fix    = fix(inWin);
    fixIdx = fixIdx(inWin);
    dur    = dur(inWin);

    assert(~isempty(fix), 'No fixation in SENTENCE_ONSET~PROMPT_ONSET window');

    fprintf('[WIN] kept %d fix | rel=[%.0f %.0f] ms | window=SENTENCE_ONSET~PROMPT_ONSET(excl)\n', ...
        numel(fix), min([fix.onset]'-tOn), max([fix.onset]'-tOn));
end
%% ================================================================

%% 2-2) 최종 좌표 (필터 후) + dx/dy(median) 재계산(표시용)
if useCorr
    x = [fix.xCorr]';  y = [fix.yCorr]';
    dxMed = median([fix.xCorr]' - [fix.x]', 'omitnan');
    dyMed = median([fix.yCorr]' - [fix.y]', 'omitnan');
else
    x = [fix.x]';      y = [fix.y]';
    dxMed = NaN; dyMed = NaN;
end
fprintf('[COORD] %s | dxMed=%.2f dyMed=%.2f (tw=%d)\n', coordMode, dxMed, dyMed, tw);

%% 3) word1 inROI QC (가능할 때만)
word1InROI = NaN; word1N = 0;
if isfield(fix,'word') && nWords>=1
    wFix = [fix.word]';
    isW1 = (wFix==1);
    word1N = sum(isW1);

    L1 = rects(1,1); T1 = rects(1,2); R1 = rects(1,3); B1 = rects(1,4);
    if word1N > 0
        word1InROI = mean( x(isW1)>=L1 & x(isW1)<=R1 & y(isW1)>=T1 & y(isW1)<=B1 );
    end
end
if isfinite(word1InROI)
    fprintf('[QC] word1 inROI = %d/%d (%.1f%%)\n', round(word1InROI*word1N), word1N, 100*word1InROI);
else
    fprintf('[QC] word1 inROI = NaN (no word mapping or no word1 fix)\n');
end

%% 4) plot
figure('Color','w'); hold on;

% ROI boxes
hROI = plot(nan,nan,'s','DisplayName','word ROI');
for i = 1:nWords
    L = rects(i,1); T = rects(i,2); R = rects(i,3); B = rects(i,4);
    rectangle('Position',[L T (R-L) (B-T)], 'LineWidth', 0.8);
end

% word labels
if ~isempty(wlist)
    minW = min(rects(:,3)-rects(:,1));
    if minW < 35
        fsz = 7;
    elseif minW < 55
        fsz = 8;
    else
        fsz = 10;
    end

    for i = 1:min(nWords,numel(wlist))
        L = rects(i,1); T = rects(i,2); R = rects(i,3); B = rects(i,4);
        text((L+R)/2, (T+B)/2, wlist(i), ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'Interpreter','none','FontSize',fsz);
    end
end

% scanpath
hScan = plot(x,y,'-','LineWidth',1,'DisplayName','scanpath');
set(hScan,'Color',scanBlue);

% fixation scatter
if strcmp(dotSizeMode,"dur") && all(~isnan(dur)) && numel(dur)>1
    d0 = max(dur,1);
    s  = dotSizeMin + (d0-min(d0))./max(eps,(max(d0)-min(d0))) * (dotSizeMax-dotSizeMin);
else
    s  = dotSizeFixed * ones(size(x));
end
hFix = scatter(x,y,s,'MarkerEdgeColor',scanBlue,'DisplayName','fixation locations');

% fixation order numbers
for i = 1:numel(x)
    text(x(i), y(i), sprintf(' %d', i), 'FontSize', 10, 'VerticalAlignment','bottom');
end

% start/end
hStart = plot(x(1),y(1),'d','MarkerSize',10,'LineWidth',2,'MarkerEdgeColor',startYellow,'DisplayName','start fixation');
hEnd   = plot(x(end),y(end),'s','MarkerSize',10,'LineWidth',2,'MarkerEdgeColor',endOrange,'DisplayName','end fixation');

% axes
set(gca,'YDir','reverse'); axis equal; grid on;
xlabel('X (px)'); ylabel('Y (px)');

% 제목: 모드 + drift + QC + 윈도우 정보
if isfinite(word1InROI)
    qcStr = sprintf('w1inROI=%.0f%%(n=%d)', 100*word1InROI, word1N);
else
    qcStr = 'w1inROI=NaN';
end

winStr = "allFix";
if onlySentenceBeforePrompt, winStr = "SENT~prePROMPT"; end

title(sprintf('Scanpath: %s (subjIdx=%d, tw=%d, nFix=%d, %s, dx=%.1f, dy=%.1f, %s, %s)', ...
    string(subj.trial(tSubj).id), tSubj, tw, numel(x), coordMode, dxMed, dyMed, qcStr, winStr));

%% ===== prePROMPT sample 표시: norm->px + (CORR면) dx/dy 적용 =====
showPromptSample = false;
if showPromptSample

    % screenW/H: subj.displayRect 우선, 없으면 results.resolution 우선 (마지막 fallback은 NaN 유지)
    screenW = NaN; screenH = NaN;
    if isfield(subj,'displayRect') && numel(subj.displayRect)==4 && all(isfinite(subj.displayRect))
        screenW = subj.displayRect(3) - subj.displayRect(1) + 1;
        screenH = subj.displayRect(4) - subj.displayRect(2) + 1;
    elseif isfield(results,'dp') && isfield(results.dp,'resolution') && numel(results.dp.resolution)>=2
        screenW = results.dp.resolution(1);
        screenH = results.dp.resolution(2);
    elseif isfield(results,'resolution') && numel(results.resolution)>=2
        screenW = results.resolution(1);
        screenH = results.resolution(2);
    end

    % PROMPT_ONSET 시간 (이미 위에서 찾았으면 재사용)
    if ~isfinite(tPrompt)
        msgIdxT = subj.trial(t).msgIdx(:);
        txtM = string({subj.event.msg(msgIdxT).text})';
        timM = [subj.event.msg(msgIdxT).time]';
        pIdx = find(contains(txtM,"PROMPT_ONSET",'IgnoreCase',true), 1, 'first');
        if ~isempty(pIdx), tPrompt = timM(pIdx); end
    end

    if isfinite(tPrompt)
        sIdx = subj.trial(t).sampleIdx(:);
        sIdx = sIdx(sIdx>=1 & sIdx<=numel(subj.sample.time));
        st = subj.sample.time(sIdx);

        near = st >= (tPrompt-50) & st < tPrompt;

        gx = subj.sample.gx(sIdx);
        gy = subj.sample.gy(sIdx);

        if any(near)
            xRaw = median(gx(near),'omitnan');
            yRaw = median(gy(near),'omitnan');

            % 정규화 판별: 0~1 범위면 norm 후보, 단 screenW/H가 있어야 px 변환 가능
            isNorm01 = isfinite(xRaw) && isfinite(yRaw) && xRaw>=0 && xRaw<=1 && yRaw>=0 && yRaw<=1;

            if isNorm01 && isfinite(screenW) && isfinite(screenH)
                xP = xRaw * screenW;
                yP = yRaw * screenH;
                srcTag = "norm01->px";
            else
                xP = xRaw;
                yP = yRaw;
                srcTag = "pxOrUnknown";
            end

            if coordMode=="CORR" && isfinite(dxMed_all) && isfinite(dyMed_all)
                xP = xP + dxMed_all;
                yP = yP + dyMed_all;
                srcTag = srcTag + "+dxdy";
            end

            if isfinite(xP) && isfinite(yP)
                hPre = plot(xP,yP,'p','MarkerSize',12,'LineWidth',2,'DisplayName','pre-PROMPT sample');
                text(xP,yP,' prePROMPT','VerticalAlignment','bottom','Interpreter','none');
                uistack(hPre,'top');
                fprintf('[PREPROMPT] %s | median x=%.3f y=%.3f -> plot x=%.1f y=%.1f | n=%d\n', ...
                    srcTag, xRaw, yRaw, xP, yP, sum(near));
            else
                fprintf('[PREPROMPT] invalid sample -> skip\n');
            end
        else
            fprintf('[PREPROMPT] no samples in window\n');
        end
    else
        fprintf('[PREPROMPT] PROMPT_ONSET not found -> skip\n');
    end
end
%% =======================================================================

legend([hROI hScan hFix hStart hEnd],'Location','bestoutside');
hold off;

