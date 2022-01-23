%%
 % Project Title: GNSS-R SDR
 % Author       : John Bagshaw
 % Contact      : jotshaw@yorku.ca
 % Supervisor   : Prof.Sunil Bisnath
 % Institution  : York University, Canada.
%%

function print_string(str)
%%% Script to print string on console given the module ID

global printDebugInfo;
global swName;

if printDebugInfo
    % https://www.mathworks.com/matlabcentral/answers/9907-
    % is-there-a-function-that-retrieves-the-filename-of-the-current-script#answer_13613
    S = dbstack();
    callerPath = which(S(2).name);
    callerLineNum = S(2).line;
    
    callerPathSpltd = regexp(callerPath, '\', 'split');
    sD = 1;
    while ~isequal(callerPathSpltd{sD}, 'mat')
        sD = sD + 1;
    end
    
    callerName = callerPathSpltd{end};
    
    callerRelPath = '~\';
    while sD < length(callerPathSpltd)
        callerRelPath = [callerRelPath, callerPathSpltd{sD}, '\'];
        sD = sD + 1;
    end
            
    callerInfo = ['path: ', callerRelPath, ',  file: ', callerName, ',  line: ', num2str(callerLineNum)];
    fprintf("%s : %-20s :> <strong>%s</strong> (%s)\n", swName, getAllFuncNames(S), str, callerInfo);
else
    S = dbstack();
    fprintf("%s : %-20s :> <strong>%s</strong> \n", swName, getAllFuncNames(S), str);
    
end
end

function recurScrptName = getAllFuncNames(S)

    sD = length(S);
    recurScrptName = S(sD).name;
    sD = sD - 1;
    while sD > 1
        recurScrptName = [recurScrptName, '->', S(sD).name];
        sD = sD - 1;
    end
end

