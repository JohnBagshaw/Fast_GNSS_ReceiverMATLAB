function acqResults = post_process_norm_acq_parcode(...
    sdrParams, ...
    dataFileParams, ...
    processing_results, ...
    rxData, ...
    delayDopplerMapAcq...
)

%%% Get Parameters from results.
averFactor         = processing_results.averFactor;
numCodeSamples     = processing_results.numCodeSamples;
numDopplerSamples  = processing_results.numDopplerBins;
dopplerResHz       = processing_results.dopplerResHz;

satelliteList       = sdrParams.sysParams.acqSatelliteList;
maxNumAcqSatellites = sdrParams.sysParams.numAcqSatellites;
numPrns             = length(satelliteList);

chipRateHz          = sdrParams.sysParams.caCodeChipRateHz;
numChipsPerMs       = floor(chipRateHz * 1e-3);

if size(delayDopplerMapAcq, 1) ~= length(satelliteList)
    error('The number of returned correlation maps does not equal number of satellites.');
end

%%% Find coarse results

peakValueList        = zeros(1, numPrns);
peakMetricList       = zeros(1, numPrns);
peakLocList          = zeros(1, numPrns);

excludeBinRangePeakAtZero = [10:floor(numCodeSamples/averFactor)-10];

% Check if a satellite should be acquired
% if its peak metric is greater than threshold.

for prnIdx = 1:numPrns
    delayDopplerMapAcqPrn = squeeze(delayDopplerMapAcq(prnIdx, :, :));
    [peakValue, peakLoc] = max(delayDopplerMapAcqPrn(:));
    
    if length(peakValue) > 1
        peakValue = peakValue(1);
        peakLoc = peakLoc(1);
    end
    
    % Calculate peak metric for all satellite signals
    dopplerShiftInt = mod(peakLoc-1, numDopplerSamples)+1;
    codeDelay       = floor(peakLoc/numDopplerSamples);
    excludeBinRange = mod(excludeBinRangePeakAtZero + codeDelay+1,...
        numCodeSamples/averFactor)+1;
    nonPeakCorrProfile = delayDopplerMapAcqPrn(dopplerShiftInt, excludeBinRange);
    secondMaxPeak = max(nonPeakCorrProfile);
    peakMetric = peakValue/secondMaxPeak;
    
    % Fill in the list of peak Metric, peak value and locations.
    peakLocList(:, prnIdx)    = peakLoc;
    peakValueList(:, prnIdx)  = peakValue;
    peakMetricList(:, prnIdx) = peakMetric;    
end

% Identify acquired and non-acquired satellites.
acqPrnIdx = (peakMetricList>sdrParams.sysParams.acqThreshold);
nacqPrnIdx = ~acqPrnIdx;

aPrnList           = satelliteList(acqPrnIdx);
aprnLocList        = peakLocList(acqPrnIdx);
aprnPeakMetricList = peakMetricList(acqPrnIdx);

nonAcqPeakMetricList = peakMetricList(nacqPrnIdx);
naprnList            = satelliteList(nacqPrnIdx);

acqPrnCount = length(aPrnList);

%%% Refine results for acquired PRNs
acqResults = struct([]);

if acqPrnCount
    
    % if acquired satellites are greater than required, choose strongest.
    if acqPrnCount > maxNumAcqSatellites
 
        % Sort the acquired satellite signals according to peak strength
        [sortedAprnPeakMetricList, ~] = ...
            sort(aprnPeakMetricList, 'descend');
        
        acqPrnIdx = ismember(aprnPeakMetricList, ...
            sortedAprnPeakMetricList(1:maxNumAcqSatellites));
        nacqPrnIdx = ~acqPrnIdx;
        
        % Choose maxNumAcqSatellites
        
        nonAcqPeakMetricList = [nonAcqPeakMetricList, aprnPeakMetricList(nacqPrnIdx)];
        naprnList            = [naprnList, aPrnList(nacqPrnIdx)];

        aPrnList           = aPrnList(acqPrnIdx);
        aprnLocList        = aprnLocList(acqPrnIdx);
        aprnPeakMetricList = aprnPeakMetricList(acqPrnIdx);

        
        acqPrnCount = maxNumAcqSatellites;
    end    
    
    print_string([num2str(acqPrnCount), ' satellites acquired for ', ...
        num2str(sdrParams.stateParams.currFrameNum+1),...
        '/', num2str(sdrParams.stateParams.numTotalFrames), '.']);
    
    % Resize if size is not ok.
    chAlgoAcqResults = cell(1, acqPrnCount);
    for aprnIdx = 1:acqPrnCount
        
        prnIdx = aPrnList(aprnIdx);
        location = aprnLocList(aprnIdx);
        
        % Coarse code and doppler delay.
        
        dopplerShiftInt = mod(location-1, numDopplerSamples)+1;
        codeDelay = floor(location/numDopplerSamples);
        dopplerSpectrum = squeeze(delayDopplerMapAcq(prnIdx, 1:numDopplerSamples, codeDelay+1));
        
        
        % First perform fine frequency estimation using quadratic
        % interpolation
        
        maxValArrIdx = mod(dopplerShiftInt-2:dopplerShiftInt, numDopplerSamples) + 1;
        maxValArr = sqrt(dopplerSpectrum(maxValArrIdx));
        doppleShiftError = 0.5*(maxValArr(1) - maxValArr(3)) / ...
            (maxValArr(1) - 2*maxValArr(2)+maxValArr(3));
        dopplerShiftInt = dopplerShiftInt - floor(numDopplerSamples/2)-1;
        dopplerShiftHz = (dopplerShiftInt + doppleShiftError) * dopplerResHz;
        dopplerShiftAxis = ((0:numDopplerSamples-1) - floor(numDopplerSamples/2))* dopplerResHz;
        
        % Refine code delay results.
        ifFreqEst = dataFileParams.intermFreqHz + dopplerShiftHz;
        
        dopplerShiftHz = 0;
        codeDelayError = 0;
        dataIdx = (codeDelay*averFactor+1):(codeDelay*averFactor+1)+numCodeSamples-1;
        
        if dataIdx(end) <= length(rxData)
            
            dopplerFreqExp  = exp(2i * pi * ifFreqEst * ...
                (0:numCodeSamples-1)/(dataFileParams.samplingFreqHz));

            rxDataMs = rxData(dataIdx) .* dopplerFreqExp;
            
            caCodeMappingInd = floor((0:numCodeSamples-1)*(chipRateHz / dataFileParams.samplingFreqHz)) + 1;
            caCodeMappingInd(caCodeMappingInd == 0) = 1;
            caCodeMappingInd(caCodeMappingInd > numChipsPerMs) = numChipsPerMs;
            caCode = gen_ca_code(sdrParams.stateParams.dataPathIn, prnIdx);
            caCode = caCode(caCodeMappingInd);
            
            [corrOut, lags] = xcorr(rxDataMs, caCode, ceil(averFactor/2));
            [~, maxIdx] = max(abs(corrOut));
            if length(maxIdx) > 1
                maxIdx = maxIdx(1);
            end
            codeDelayError = lags(maxIdx);
                        
            % Refine doppler frequency agains.
            
            rxDataMs = rxData(dataIdx) .* caCode;
            rxDataMs = rxDataMs .* dopplerFreqExp;
            
            acqDopplerBW = sdrParams.sysParams.acqDopplerBwKhz * 1e3;
            fftNumPts = 8*2^(nextpow2(length(rxDataMs)));
            deltaF = dataFileParams.samplingFreqHz/fftNumPts;
            pbins = ceil(0.5 * acqDopplerBW / deltaF);
            faxis = 0.5*fftNumPts-pbins:pbins+0.5*fftNumPts;
            fftFreqBins = (faxis - floor(0.5*fftNumPts)-1) * deltaF;
            
            fftxc = abs(fftshift(fft(rxDataMs, fftNumPts)));
            fftxc = fftxc(faxis);
            [~, fftMaxIndex] = max(fftxc);
            
            maxValArrIdx = fftMaxIndex-1:fftMaxIndex+1;
            maxValArr = sqrt(fftxc(maxValArrIdx));
            
            % Apply quadratic interpolation again.
            
            doppleShiftError1 = 0.5*(maxValArr(1) - maxValArr(3)) / ...
                (maxValArr(1) - 2*maxValArr(2)+maxValArr(3));
            dopplerShiftHz  = fftFreqBins(fftMaxIndex) + doppleShiftError1*deltaF;
            
        end
        
        codeDelay = codeDelayError + codeDelay * averFactor + 1;
        ifFreqEst   = ifFreqEst - dopplerShiftHz;
        
        % Save results per PRN
        chAlgoAcqResults{aprnIdx}.dopplerShiftHz = dopplerShiftHz;
        chAlgoAcqResults{aprnIdx}.codeDelay      = codeDelay;
        chAlgoAcqResults{aprnIdx}.peakMetric     = aprnPeakMetricList(aprnIdx);
        chAlgoAcqResults{aprnIdx}.dopplerShiftHz = ifFreqEst;
        chAlgoAcqResults{aprnIdx}.satellitePrn   = prnIdx;
        chAlgoAcqResults{aprnIdx}.ddm            = squeeze(delayDopplerMapAcq(prnIdx, :, :));
        chAlgoAcqResults{aprnIdx}.dopplerShiftAxis=dopplerShiftAxis;
        chAlgoAcqResults{aprnIdx}.codeDelayAxis  = 0:averFactor:numCodeSamples-1;

    end
    
    acqResults(1).chAlgoAcqResults = chAlgoAcqResults;
    acqResults(1).nacqStats.nacqPrnList = naprnList;
    acqResults(1).nacqStats.nacqPrnPeakMetric = nonAcqPeakMetricList;
end
end