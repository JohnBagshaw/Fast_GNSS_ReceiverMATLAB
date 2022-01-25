function save_acquisition_resutls(sdrParams, postProcessResults)
%%%

%%% Extract coherent frame data.
currFileNum = sdrParams.stateParams.numFilesProcessed + 1;
fileName    = sdrParams.stateParams.fileNames{currFileNum};
dataParams  = sdrParams.dataFileParamsList{currFileNum};

channelId = 1;
if dataParams.selectedChannel ~= -1
    numChannels = 1;
    channelId = dataParams.selectedChannel;
else
    numChannels = dataParams.totalChannels;
    channelId = 1:dataParams.totalChannels;
end
                 
acqAlgoList = sdrParams.sysParams.acqAlgosList;
numAcqAlgos = length(acqAlgoList);

print_string('====Displaying Acquisition Results: Start====');
for chIdx=1:numChannels
    for acqAlgoIdx=1:numAcqAlgos
        print_string(sprintf('(File: %s, Channel: %d, Algorithm: %s)', ...
            fileName, chIdx, acqAlgoList{acqAlgoIdx}));
        
        acqResults = postProcessResults.acqResults{chIdx, acqAlgoIdx}.acqResults;        
         if ~isempty(acqResults)
             
             chAlgoAcqResults = acqResults.chAlgoAcqResults;

             %%% Save
             sFileName = sprintf('file_%s_ch_%d_frame_%d_algo_%s_%s.mat', ...
                 fileName(1:end-4),...
                 chIdx,...
                 sdrParams.stateParams.currFrameNum+1, ...
                 acqAlgoList{acqAlgoIdx},...
                 datestr(now, 'yyyy-mm-dd-HH-MM-SS'));

             % clean dir
             
             if currFileNum == 1 && ...
                chIdx == 1 && ...
                sdrParams.stateParams.currFrameNum == 0 &&...
                acqAlgoIdx == 1
                 
                 
                 % Check to make sure that folder actually exists.  Warn user if it doesn't.
                 if ~isdir(sdrParams.stateParams.dataPathOut)
                     errorMessage = sprintf('Error: The following folder does not exist:\n%s', ...
                         sdrParams.stateParams.dataPathOut);
                     uiwait(warndlg(errorMessage));
                     return;
                 end
                 % Get a list of all files in the folder with the desired file name pattern.
                 filePattern = fullfile(sdrParams.stateParams.dataPathOut, '*.*');
                 theFiles = dir(filePattern);
                 for k = 3 : length(theFiles)
                     baseFileName = theFiles(k).name;
                     fullFileName = fullfile(sdrParams.stateParams.dataPathOut, baseFileName);
                     print_string(sprintf('Now deleting %s', fullFileName));
                     delete(fullFileName);
                 end
             end
             % save
             save([sdrParams.stateParams.dataPathOut, sFileName], 'acqResults');
             

             
             %%% Dump to console
             for prn=1:length(chAlgoAcqResults)
                 print_string(sprintf('(PRN: %3d, Metric: %8s, Code: %5d, Frequency: %s)', ...
                 chAlgoAcqResults{prn}.satellitePrn,...
                 num2str(chAlgoAcqResults{prn}.peakMetric),...
                 chAlgoAcqResults{prn}.codeDelay,...
                 num2str(chAlgoAcqResults{prn}.dopplerShiftHz)));
             end   
        else
            print_string('No satellites found.');
        end
    end
end
print_string('====Displaying Acquisition Results: End====');
end
