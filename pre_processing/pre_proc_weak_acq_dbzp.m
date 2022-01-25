function pprocSignals = pre_proc_weak_acq_dbzp(sdrParams, caCode)

%%% Parameters
fileNum         = sdrParams.stateParams.numFilesProcessed + 1;
samplingFreqHz  = sdrParams.dataFileParamsList{fileNum}.samplingFreqHz;
intermFreqHz    = sdrParams.dataFileParamsList{fileNum}.intermFreqHz;
chipRateHz      = sdrParams.sysParams.caCodeChipRateHz;
acqDopplerBwKhz = sdrParams.sysParams.acqDopplerBwKhz;
cpiTimeMs       = sdrParams.sysParams.coherentProcessingTimeMS;

samplesPerCpi = samplingFreqHz * cpiTimeMs *1e-3;
samplesPerMs  = samplingFreqHz * 1e-3;
numChipsPerMs = chipRateHz * 1e-3;
numCaCodeFolds = cpiTimeMs;
prnList = sdrParams.sysParams.acqSatelliteList;

if floor(samplesPerCpi) ~= samplesPerCpi
    error('Number of samples per coherent processing interval can''t have decimal point.');
end

% Prepare per algorithm per block common signals

%%% Determine maximum extent to which averaging can be done
averFactor = [1];
temAverFactor = 1;
while samplingFreqHz/temAverFactor > sdrParams.sysParams.minSamplingFreqHz
    temAverFactor=temAverFactor+1;
    goodAverFactor = samplesPerCpi/temAverFactor;
    if floor(goodAverFactor) == goodAverFactor
        averFactor = [averFactor, temAverFactor];
    end
end
% choose highest dividing factor to maximize number of input samples
% for one averaged output.
averFactor = averFactor(end);


%%% Find out integer number of samples per block and doppler bins

numSamplesPerCpi = samplesPerCpi / averFactor;
%number of doppler bins should be atleast 2/t_cpi
minDopplerBins = round(acqDopplerBwKhz * cpiTimeMs);
candidateNumDopplerBins = [];
for ndb=minDopplerBins:minDopplerBins*2
    if floor(numSamplesPerCpi/ndb) == (numSamplesPerCpi/ndb)
        candidateNumDopplerBins = [candidateNumDopplerBins, ndb];
    end
end
if isempty(candidateNumDopplerBins)
    error('Unable to find integer multiple of blockSize for DBZP.');
end
numDopplerBins = min(candidateNumDopplerBins);
numSamplesPerBlock = numSamplesPerCpi/numDopplerBins;



%%% generate code for current PRN
% Convert CA code to DBZP format.
    
%%% Map to samples of 1ms
caCodeSampleMap = floor((0:samplesPerMs/averFactor-1) * (numChipsPerMs*averFactor/samplesPerMs)) + 1;
caCodeMapped = caCode(prnList, caCodeSampleMap);
caCodeMapped = repmat(caCodeMapped, 1, numCaCodeFolds); 

numBlocks = samplesPerCpi/averFactor/numSamplesPerBlock;
caCodeDbzpTable = zeros(length(prnList), 2*numSamplesPerBlock, numBlocks);
for prnIdx=1:length(prnList)
    caCodePerPrn = reshape(caCodeMapped(prnIdx, :), numSamplesPerBlock, numBlocks);
    caCodePerPrn = [caCodePerPrn; zeros(numSamplesPerBlock, numBlocks)];
    caCodeDbzpTable(prnIdx, :, :) = caCodePerPrn;
end


%%% Prepare baseband doppler modulated matrix of data

dopplerFreqExp  = exp(2i * pi * (intermFreqHz/samplingFreqHz) * (0:samplesPerCpi-1));


%%% Save preprocessing common signals and parameters to used in the 
% processing of DBZP algorithm.

pprocSignals.caCodesTable       = caCodeDbzpTable;
pprocSignals.dopplerFreqExp     = dopplerFreqExp;
pprocSignals.numBlocks          = numBlocks;
pprocSignals.numSamplesPerBlock = numSamplesPerBlock;
pprocSignals.averFactor         = averFactor;

end