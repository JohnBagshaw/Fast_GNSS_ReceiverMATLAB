function acqResults = norm_acq_parcode(sdrParams, ppData, rxFrameData)
%%% This function runs acquisition algorithm known as 
% 'Parallel Code Search Acquisition'. The minimum processing
% block is for 1ms and this function takes care to give 
% results for coherentIntegrationInterval results. Results is a
% Nd x N matrix containing complex correlation matrix for Nd doppler bins
% and N samples per code. It also contains relevant stats used for 
% optimizations for examples number of samples averaged etc.
% acqResults.ddMap;
% acqResults.stats;

optimizeOption = ppData.optimizeOption ;
tic;

if optimizeOption == 1
    
% Parameters
averFactor = ppData.averFactor;
prnList = sdrParams.sysParams.acqSatelliteList;
numPrns = length(prnList);
numBlocks = sdrParams.sysParams.coherentProcessingTimeMS;
currFile = sdrParams.stateParams.numFilesProcessed+1;
numCodeSamples = sdrParams.dataFileParamsList{currFile}.samplingFreqHz * 1e-3;
numDopplerSamples = floor(sdrParams.sysParams.acqDopplerBwKhz * 1e3 / ...
                          sdrParams.sysParams.acqDopplerResHz) + 1;

rxDataStartIdx = 1;

% Define buffer for correlation matrix
delayDopplerCohMat = zeros(numPrns, numDopplerSamples, numCodeSamples/averFactor);
circCorrOutMat = zeros(numPrns, numDopplerSamples, numCodeSamples/averFactor);

% Iterate over each block to do the processing.
for blIdx=1:numBlocks    
    
    
    % Multiply carrier signal with block data.
    blockData = rxFrameData(rxDataStartIdx:rxDataStartIdx+numCodeSamples-1) .*...
                ppData.dopplerFreqExp;

    % Take average of consecutive samples.
    blockDataN = zeros(numDopplerSamples, numCodeSamples/averFactor);
    for dbin=1:numDopplerSamples
        blockVec = blockData(dbin, :);
        blockVec = reshape(blockVec, averFactor, numCodeSamples/averFactor);
        blockDataN(dbin, :) = mean(blockVec, 1);        
    end
    if 0
        blockData = permute(reshape(blockData, numDopplerSamples, averFactor, ...
                    numCodeSamples/averFactor), [1,3,2]);
        blockData = mean(blockData, 3);
    end

    
    % Do the correaltion for 32 PRn iterations with each 
    % circular correlation using FFT
    for prn=1:length(prnList)
        blockDataModFd   = fft(blockDataN, [], 2);
        circCorrOutMatFd = blockDataModFd .* ppData.caCodesTable(prn,:);
        circCorrOutMat(prn,:,:) = ifft(circCorrOutMatFd, [], 2);
    end
                
    % Add result to correlation matrix
    delayDopplerCohMat = delayDopplerCohMat + abs(circCorrOutMat).^2;
    rxDataStartIdx = rxDataStartIdx + numCodeSamples;
 

end
acqResults.ddMap = delayDopplerCohMat;
acqResults.averFactor = averFactor;
acqResults.dopplerResHz = ppData.dopplerResHz;

elseif optimizeOption == 2
    
    
% Parameters

averFactor = ppData.averFactor;
prnList = sdrParams.sysParams.acqSatelliteList;
numPrns = length(prnList);
numBlocks = sdrParams.sysParams.coherentProcessingTimeMS;
currFile = sdrParams.stateParams.numFilesProcessed+1;
numCodeSamples = sdrParams.dataFileParamsList{currFile}.samplingFreqHz * 1e-3;
numDopplerSamples = ppData.numDopplerSamples;
rxDataStartIdx = 1;

% Define buffer for correlation matrix
delayDopplerCohMat = zeros(numPrns, numDopplerSamples, numCodeSamples/averFactor);
circCorrOutMat = zeros(numPrns, numDopplerSamples, numCodeSamples/averFactor);

% Iterate over each block to do the processing.
for blIdx=1:numBlocks    
    
    blockData = rxFrameData(rxDataStartIdx:rxDataStartIdx+numCodeSamples-1);
    
    % Multiply carrier signal with block data.
    blockData = blockData .* ppData.dopplerFreqInitExp;
    blockData  = reshape(blockData, averFactor, numCodeSamples/averFactor);
    blockData = mean(blockData, 1);
    blockData = fft(blockData);
     
    blockDataN = zeros(numDopplerSamples, numCodeSamples/averFactor);
    if 1
        shiftFactor = ppData.shiftFactor;
        sh=0;
        for dbin=1:numDopplerSamples
            blockDataN(dbin, :) = circshift(blockData, sh);
            sh = floor(shiftFactor * dbin);
        end
    else
%         r=blockData;
%         c=circshift(fliplr(blockData), 1);
%         c = c(1:numDopplerSamples);
%         blockDataN = toeplitz(c, r);        

    numCols = numCodeSamples/averFactor;
    blockDataN = repmat(blockData, numDopplerSamples, 1);
    blockDataN = ifft(fft(blockDataN,[],2) .* ...
                 exp(-2i*pi/numCols*(0:numDopplerSamples-1)'*...
                 (0:numCols-1)) ,[],2);

    end
    
    
    % Do the correaltion for 32 PRn iterations with each 
    % circular correlation using FFT
    for prn=1:length(prnList)
        circCorrOutMatFd = blockDataN .* ppData.caCodesTable(prn,:);
        circCorrOutMat(prn,:,:)   = ifft(circCorrOutMatFd, [], 2);
    end
    % Add result to correlation matrix
    delayDopplerCohMat = delayDopplerCohMat + abs(circCorrOutMat).^2;
    rxDataStartIdx = rxDataStartIdx + numCodeSamples;
end

acqResults.ddMap = delayDopplerCohMat;
acqResults.averFactor = averFactor;
acqResults.numCodeSamples = numCodeSamples;
acqResults.numDopplerBins = numDopplerSamples;
acqResults.dopplerResHz = ppData.dopplerResHz;
end
end

