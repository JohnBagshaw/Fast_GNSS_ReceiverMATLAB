function [ postProcessResuts ] = post_process(sdrParams, rxData, ...
    processing_results...
)
% SUMMARY Post processes correlation for 1 frame of coherent processing interval.
% This function receives processing results for configured number of
% satellites for 1 coherent processing interval and converts correlation
% delay doppler profile to code delay and doppler shift parameters as well
% as refines them through advanced signal processing techniques. It stores
% results for one incoheren processing interval frame and then outputs them
% either to plotting functions or dumps to files. 

postProcessResuts = struct();

currFileNum         = sdrParams.stateParams.numFilesProcessed+1;
dataFileParams      = sdrParams.dataFileParamsList{currFileNum};
numChannels         = dataFileParams.totalChannels;
selectedChannel     = dataFileParams.selectedChannel;
acqAlgos            = sdrParams.sysParams.acqAlgosList;
numAcqAlgos         = length(acqAlgos);

if selectedChannel ~= -1
    numChannels = 1;
end

%%% Define buffer for incoherent integration for all data channels 
% and all acquisition algorithms types.

persistent ddMapBuffer;
if isempty(ddMapBuffer)
    ddMapBuffer = cell(numChannels, numAcqAlgos);
end

% Define data buffer here that are static
acqResults = cell(numChannels, numAcqAlgos);


for chIdx = 1:numChannels
    for algoIdx = 1:numAcqAlgos 
                
        %%% Add the delay doppler map for multiple incoherent integration
        % intervals.
        
        if isempty(ddMapBuffer) || ...
           sdrParams.stateParams.currFrameNum == 0
            ddMapBuffer{chIdx, algoIdx}.ddMapBuffer = processing_results{chIdx, algoIdx}.ddMap;
        else
            ddMapBuffer{chIdx, algoIdx}.ddMapBuffer = ddMapBuffer{chIdx, algoIdx}.ddMapBuffer + ...
                                                      processing_results{chIdx, algoIdx}.ddMap;
        end
        
        
        print_string(sprintf('Post processing for frame: %d/%d, channel: %d/%d, algorithm: %s', ...
        sdrParams.stateParams.currFrameNum+1, ...
        sdrParams.stateParams.numTotalFrames, ...
        chIdx, ...
        numChannels, ...
        acqAlgos{algoIdx}...
        ));

        
        %%% Do post-processing for each algorithm.
        switch(acqAlgos{algoIdx})
            case 'norm_acq_parcode'
                chAlgoAcqResults = post_process_norm_acq_parcode(...
                    sdrParams, ...
                    dataFileParams, ...
                    processing_results{chIdx, algoIdx}, ...
                    rxData{chIdx}, ...
                    ddMapBuffer{chIdx, algoIdx}.ddMapBuffer...
                    );
            case  'weak_acq_dbzp'
                chAlgoAcqResults = post_process_weak_acq_dbzp( ... 
                    sdrParams, ...
                    dataFileParams, ...
                    processing_results{chIdx, algoIdx}, ...
                    rxData{chIdx}, ...
                    ddMapBuffer{chIdx, algoIdx}.ddMapBuffer...                    
                    );
            otherwise
                error('Acquisition algorithm not identified.');
        end
        
        acqResults{chIdx, algoIdx}.acqResults = chAlgoAcqResults;
    end
    
    % Save results as well as spit out.

end
  

%%% visualizations of different algorithms
% plot_acquisition(acqResults);
    postProcessResuts.acqResults = acqResults;
    
end
