function d = cohens_d(x, y)
    x = x(~isnan(x));
    y = y(~isnan(y));

    nx = numel(x);
    ny = numel(y);
    sx = std(x);
    sy = std(y);

    sp = sqrt(((nx-1)*sx^2 + (ny-1)*sy^2) / (nx+ny-2));  % pooled SD
    d  = (mean(x) - mean(y)) / sp;
end

%% 분석 결과
% Cohen''s d (goPast, N word): 
% N-P=0.218 → “작은 효과 +”
% N-U=0.331 → “작고-중간 사이 정도”
% P-U=0.129→ “아주 작은 효과”

% 부정(N) vs 긍정(P/중립 U):
% → go-past time에서 일관되게 더 느린 방향의 차이는 있고,
% → 크기는 “작지만 무시할 정도는 아닌” 수준
% 긍정 vs 중립 (P vs U)는 거의 차이가 없음
% 패턴은 "부정 단어만 특별히 더 오래 붙잡는다”는 자동 경계 가설 방향과 잘 맞음


% Cohen''s d (regPathDur, N+1): 
% N-P=0.040 → 매우 작은 효과(negligible effect)
% N-U=0.105 → small
% P-U=0.071→ medium
% → N, P, U 간 spillover(regPathDur)에서 큰 차이가 없음

% N 단어에서의 초기 처리(goPast)는 valence 차이가 뚜렷하지만,
% N+1 spillover 단계에서는 valence 효과가 거의 소멸되는 패턴.
% eye-tracking reading 연구에서 꽤 흔한 구조

% 자동 경계 가설 관점에서 설명하면
% Automatic Vigilance Hypothesis는:
% 부정 단어(N)에 빠르게 주의가 붙고
% 주의를 떼는 데 시간(cost)이 더 든다
% spillover(N+1)는 개별 연구마다 양상이 다름

% Target(N):goPast 효과 있음 → 부정 단어에서 초기 처리 비용 증가
% Spillover(N+1): regPathDur 효과 거의 없음 → 부정 단어 효과가 N+1로 carry-over되지 않음
% → 초기 포착 강하지만, spillover는 제한적