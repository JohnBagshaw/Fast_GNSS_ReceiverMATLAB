%%
 % Project Title: GNSS-R SDR
 % Author       : John Bagshaw
 % Contact      : jotshaw@yorku.ca
 % Supervisor   : Prof.Sunil Bisnath
 % Institution  : York University, Canada.
%%

%%%================== GNSS SDR ===========================%%%
%                                                       
% @brief: Software Defined GNSS Receiver Software Model 
% @date: 11-19-2021                                     
%%%=======================================================%%%

%%% Clean up the environment first
clear; close all; clc; clear init_settings;

%%% Formats
format ('compact');
format ('long', 'g');

%%% Global variables
global printDebugInfo;
global swName;
swName = 'JohnBagshawGNSS_SDR_v0.1';
printDebugInfo = 0;

%%% Add all folders and sub folders in current path.
addpath(genpath('.'));

%%% Print startup
print_string("===============================================================");
print_string(['Welcome! MATLAB Reference Software-Defined Radio: ', swName]);
print_string("===============================================================");

% Do parameter configuration.
sdrParams = config_sdr_params();


% Do multiple data stream processing (offline)
% TODO: support real-time.

while sdrParams.stateParams.numFilesProcessed < ...
    sdrParams.stateParams.numFilesToProcess

    print_string('-----------------------------------------------------------------');
    print_string(['Data processing started for file: ',...
        sdrParams.stateParams.fileNames{sdrParams.stateParams.numFilesProcessed+1}]);
    print_string('-----------------------------------------------------------------');


    % Reset frame counter
    sdrParams.stateParams.currFrameNum = 0;
        
    %%% pre-processing per input data file.
    [ppData, rxData] = pre_process(sdrParams);
    
    %%% Data frame processing
    while sdrParams.stateParams.currFrameNum < ...
        sdrParams.stateParams.numTotalFrames

        %%% Processing
        % Process function is based on processing of data frames per
        % currently processing file.
        processResults = process(sdrParams, ppData, rxData);
        
        %%% post-processing
        % Post processing is done for data frames with one time final
        % global post processing managed internally.
        postProcessResults = post_process(sdrParams, rxData, processResults);
        
        sdrParams.stateParams.currFrameNum = ...
            sdrParams.stateParams.currFrameNum + 1;        
    end
    
    
    %%% Plotting the results.
    plot_acquisition_results(sdrParams, postProcessResults);

    
    %%% Save results to configured file format.
    save_acquisition_resutls(sdrParams, postProcessResults);

    
    print_string('-----------------------------------------------------------------');
    print_string(['Data processing completed for file: ',...
        sdrParams.stateParams.fileNames{sdrParams.stateParams.numFilesProcessed+1}]);
    print_string('-----------------------------------------------------------------');
    
    % Next file.
    sdrParams.stateParams.numFilesProcessed = ...
        sdrParams.stateParams.numFilesProcessed + 1;

end
print_string("============================");
print_string("Program Completed. Good Bye!");
print_string("============================");

