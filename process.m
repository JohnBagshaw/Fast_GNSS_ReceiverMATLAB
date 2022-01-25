function [ processResults] = process( sdrParams, ppData, rxData )
% SUMMARY Function to process coherent processing interval data.
% This function takes in input Raw data, preprocessed common signals
% and processes current frame (coherent processing interval for current
% input data file. process_results returns cell of length C where C is
% number of channels and each cell contains L cells each corresponding to 
% different algorithm and each each contains struct containing processing
% results of delay doppler map and acquisition stats.



%%% Extract coherent frame data.
currFileNum = sdrParams.stateParams.numFilesProcessed + 1;
dataParams = sdrParams.dataFileParamsList{currFileNum};
numSamplesPerMs = dataParams.samplingFreqHz * 1e-3;
currFrameNum    = sdrParams.stateParams.currFrameNum + 1;
numTotalFrames  = sdrParams.stateParams.numTotalFrames;
numSamplesPerFrame = sdrParams.sysParams.coherentProcessingTimeMS * numSamplesPerMs;
frameDataIndex = (currFrameNum-1)*numSamplesPerFrame + 1: ...
             currFrameNum * numSamplesPerFrame ;
numChannels = length(rxData);
acqAlgoList = sdrParams.sysParams.acqAlgosList;
numAcqAlgos = length(acqAlgoList);
processResults = cell(numChannels, numAcqAlgos);
for chIdx=1:numChannels
    
    % Rx data for current frame
    % used for coherent integration 
    % within the function.
    frameData = rxData{chIdx}(frameDataIndex);
        
    % Iterate over all algorithms and process the data.
    % These are core algorithms for acquisition.
    for acqAlgoIdx = 1:numAcqAlgos
        
        
        %%% Acquisition

        print_string(sprintf('Processing for frame: %d/%d, channel: %d/%d, algorithm: %s', ...
        sdrParams.stateParams.currFrameNum+1, ...
        sdrParams.stateParams.numTotalFrames, ...
        chIdx, ...
        numChannels, ...
        acqAlgoList{acqAlgoIdx}...
        ));

     
        clear( acqAlgoList{acqAlgoIdx} ); % Clear previous calls  
        
        processResults{chIdx, acqAlgoIdx} = feval(acqAlgoList{acqAlgoIdx}, ...       % Call acqusition function 
                                            sdrParams, ...
                                            ppData{acqAlgoIdx}, ...
                                            frameData);
                                        
                                        
        %%% Tracking
        

                                            
    end
end    
    
    
end


