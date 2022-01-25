function pprocSignals = pre_proc_norm_acq_parcode(sdrParams, caCode)

optimizeOption = 2; % 2 for double efficiency
if optimizeOption ~= 1 && optimizeOption ~= 2
    error('Only two optimization options are available for parallel code search algorithm.');
end

%%% Parameters
fileNum = sdrParams.stateParams.numFilesProcessed + 1;
samplingFreqHz = sdrParams.dataFileParamsList{fileNum}.samplingFreqHz;
intermFreqHz = sdrParams.dataFileParamsList{fileNum}.intermFreqHz;
chipRateHz = sdrParams.sysParams.caCodeChipRateHz;
acqDopplerBwKhz = sdrParams.sysParams.acqDopplerBwKhz;
acqDopplerResHz = sdrParams.sysParams.acqDopplerResHz;

samplesPerMs = samplingFreqHz * 1e-3;
numChipsPerMs = chipRateHz * 1e-3;
prnList = sdrParams.sysParams.acqSatelliteList;

if floor(samplesPerMs) ~= samplesPerMs
    error('Number of samples per millisecond can''t have decimal point.');
end

% Prepare per algorithm per block common signals

%%% Determine maximum extent to which averaging can be done
averFactor = [1];
temAverFactor = 1;
while samplingFreqHz/temAverFactor > sdrParams.sysParams.minSamplingFreqHz
    temAverFactor=temAverFactor+1;
    goodAverFactor = samplesPerMs/temAverFactor;
    if floor(goodAverFactor) == goodAverFactor
        averFactor = [averFactor, temAverFactor];
    end
end

averFactor = averFactor(end);

%%% each PRN C/A code sequence frequency transform here.
caCodeMappingInd = floor((0:samplesPerMs/averFactor-1) * ...
    (chipRateHz *averFactor/ samplingFreqHz)) + 1;

caCodeMappingInd(caCodeMappingInd == 0) = 1;
caCodeMappingInd(caCodeMappingInd > numChipsPerMs) = numChipsPerMs;
caCodesTable = caCode(:, caCodeMappingInd);
caCodesTable = conj(fft(caCodesTable, [], 2));


% Prepare baseband doppler modulated matrix of data
if optimizeOption == 1
    numDopplerSamples = floor(acqDopplerBwKhz * 1e3 / acqDopplerResHz) + 1;
    intermFreqVec   = intermFreqHz;
    dopplerFreqVec  = acqDopplerResHz*(0:numDopplerSamples-1);
    dopplerFreqVec  = dopplerFreqVec - median(dopplerFreqVec); % make it two sided
    iFfreqVec = intermFreqVec + dopplerFreqVec;
    dopplerFreqExp  = exp(2i * pi * iFfreqVec' * (0:samplesPerMs-1)...
                          /(samplingFreqHz));    
                      
    pprocSignals.dopplerFreqExp = dopplerFreqExp;
    pprocSignals.dopplerResHz = acqDopplerResHz;

elseif optimizeOption == 2
    
    % In this method, frequency shift is integer in terms of 1/1ms=1khz
    acqDopplerResHz = 1e3; % 1/0.001 for time duration of 1 millisecond
    numDopplerSamples = floor(acqDopplerBwKhz * 1e3 / acqDopplerResHz) + 1;    
    iFfreqVec = intermFreqHz - acqDopplerBwKhz * 0.5e3;
    dopplerFreqInitExp = exp(2i * pi * iFfreqVec * (0:samplesPerMs-1)...
        /(samplingFreqHz));
    dopplerFreqDeltaExp = exp(2i * pi * acqDopplerResHz * (0:samplesPerMs-1)...
        /(samplingFreqHz));
    
    pprocSignals.shiftFactor = acqDopplerResHz * samplesPerMs / ...
                                                 samplingFreqHz;
    pprocSignals.dopplerFreqInitExp = dopplerFreqInitExp;
    pprocSignals.dopplerFreqDeltaExp = dopplerFreqDeltaExp;
    pprocSignals.numDopplerSamples = numDopplerSamples;
    pprocSignals.dopplerResHz = acqDopplerResHz;
        
end

pprocSignals.caCodesTable = caCodesTable;
pprocSignals.averFactor = averFactor;
pprocSignals.optimizeOption = optimizeOption;


end