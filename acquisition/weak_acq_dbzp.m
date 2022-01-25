function acqResults = weak_acq_dbzp(sdrParams, ppData, rxFrameData)
%%% This function runs acquisition algorithm known as 
% 'Double block zero padding weak acquisition'. The minimum processing
% block is for 1ms and this function takes care to give 
% results for coherentIntegrationInterval results. Results is a
% Nd x N matrix containing complex correlation matrix for Nd doppler bins
% and N samples per coherent processing interval divided by averaging factor.
% It also contains relevant stats used for otimizations for examples number 
% of samples averaged etc.
% acqResults.ddMap;
% acqResults.stats;


%%% Parameters

caCodesTable       = ppData.caCodesTable;
dopplerFreqExp     = ppData.dopplerFreqExp;
numBlocks          = ppData.numBlocks;
numSamplesPerBlock = ppData.numSamplesPerBlock;
averFactor         = ppData.averFactor;
prnList            = sdrParams.sysParams.acqSatelliteList;
currFile           = sdrParams.stateParams.numFilesProcessed+1;
samplingFreqHz     = sdrParams.dataFileParamsList{currFile}.samplingFreqHz;
numPrns            = length(prnList);
numDopplerSamples  = numBlocks;
numCohIntMs        = sdrParams.sysParams.coherentProcessingTimeMS;

dopplerDftInterpFactor = 4;
numDopplerFftBins      = dopplerDftInterpFactor*numDopplerSamples;

% Define buffer for correlation matrix
circCorrPartial = zeros(numSamplesPerBlock*numDopplerSamples, numDopplerSamples);
circCorrOutMat  = zeros(numPrns, numSamplesPerBlock*numDopplerSamples/numCohIntMs, numDopplerFftBins);


%%% Convert signal to baseband.
dataBB = rxFrameData .* dopplerFreqExp;


%%% Average the input signal
dataBB = reshape(dataBB, averFactor, length(dataBB)/averFactor);
dataBB = mean(dataBB, 1);

%%% Convert signal to DBZP format matix.
dataBB = reshape(dataBB, numSamplesPerBlock, numBlocks);
dataBB = [dataBB; circshift(dataBB, -1, 2)];

% dataBB = repmat(dataBB, 1, numBlocks);
% r=1:numBlocks;
% c=circshift(fliplr(1:numBlocks), 1);
% matIdx = toeplitz(c, r)';   
% permIdx = matIdx(:);
% caCodeDbzpPerm = squeeze(caCodesTable(prnIdx, :, permIdx));
% parCorrOutput = ifft(fft(dataBB) .* conj(fft(caCodeDbzpPerm)));
% parCorrOutput = parCorrOutput(1:numSamplesPerBlock, :);
% parCorrOutput = mat2cell(parCorrOutput, [numSamplesPerBlock], [53*ones(1, 53)])';
% parCorrOutput = cell2mat(parCorrOutput);

%%% Perform search for all listed PRN numbers ..
for prnIdx = 1:numPrns
    
    %%% Generate local PRN C/A code

    %%% generate code for current PRN
    caCode = squeeze(caCodesTable(prnIdx, :, :));
    
    %%% Iterate over each permuation
    for permIdx = 1:numBlocks
        
        %%% Generate Permuted DBZP matrix
        permMap = circshift(1:numBlocks, permIdx-1, 2);
        caCodeDbzpedPerm = caCode(:, permMap);
        
        %%% Partial correlation using FFT
        parCorrOutput = ifft(fft(dataBB) .* conj(fft(caCodeDbzpedPerm)));
        
        %%% Put partial correlation results in the buffer
        
        circCorrPartial((permIdx-1)*numSamplesPerBlock+1:permIdx*numSamplesPerBlock, :) = ...
            parCorrOutput(1:numSamplesPerBlock, :);
    end
    
    %%% FFT across the columsn (doppler bins)
    
    dftCircCorrPartial = abs(fftshift(fft(circCorrPartial, numDopplerFftBins, 2), 2)).^2;
    dftCircCorrPartialCellArr = mat2cell(dftCircCorrPartial, ...
        size(dftCircCorrPartial, 1)/numCohIntMs*ones(1, numCohIntMs), size(dftCircCorrPartial, 2));
    
    dftCircCorrPartial = zeros(size(dftCircCorrPartialCellArr{1}));
    for m=1:numCohIntMs
        dftCircCorrPartial = dftCircCorrPartial + dftCircCorrPartialCellArr{m};
    end    
    circCorrOutMat(prnIdx, :, :) = dftCircCorrPartial;    
end

% Pack results in the output.
acqResults.ddMap = circCorrOutMat;
acqResults.averFactor = averFactor;
acqResults.numSamplesPerBlock = numSamplesPerBlock;
acqResults.numBlocks = numBlocks/numCohIntMs;
acqResults.numDopplerBins = numDopplerSamples;
acqResults.dopplerResHz = (samplingFreqHz/averFactor/numSamplesPerBlock) / numDopplerFftBins;
acqResults.numDopplerFftBins = numDopplerFftBins;

end


