function key = makeSentenceKey(s)
    if ismissing(s) || strlength(string(s))==0
        key = "";
        return;
    end
    s = string(s);

    % 보이지 않는 문자/제어문자 제거 + 공백 통일
    s = replace(s, char(160), ' ');      % NBSP
    s = regexprep(s, '[\t\r\n]+', ' ');  % 탭/줄바꿈
    s = regexprep(s, '\s+', ' ');        % 연속 공백
    s = strtrim(s);

    % 스마트 따옴표/하이픈 통일
    s = replace(s, ["“","”","„","″"], '"');
    s = replace(s, ["‘","’","‚","′"], "'");
    s = replace(s, ["–","—","−"], "-");

    % 종결부호 흡수(원하면 주석처리 가능)
    s = regexprep(s, '[\.\!\?…]+$', '');
    s = strtrim(s);

    key = lower(s);
end
