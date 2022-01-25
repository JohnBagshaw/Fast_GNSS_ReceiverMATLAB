%function acqResults = delay_doppler_map(longSignal, settings)

if isequal(settings.acquisitionType, 'weakAcquisition')
    
    %% Process delay doppler map through weak signal acquisition
    
    longSignal = data;
    
    samplesPerCode = round(settings.samplingFreq / ...
        (settings.codeFreqBasis / settings.codeLength));

    % Create two 1msec vectors of data to correlate with and one with zero DC
    signal0DC = longSignal - mean(longSignal); 

    % Find sampling period
    ts = 1 / settings.samplingFreq;

    % Find phase points of the local carrier wave 
    phasePoints = (0 : (samplesPerCode-1)) * 2 * pi * ts;

    % Number of the frequency bins for the given acquisition band (500Hz steps)
    numberOfFrqBins = round(settings.acqSearchBand * 2) + 1;

    % Generate all C/A codes and sample them according to the sampling freq.
    caCodesTable = make_ca_table(settings);


    %--- Initialize acqResults ------------------------------------------------
    % Carrier frequencies of detected signals
    acqResults.carrFreq     = zeros(1, 32);
    % C/A code phases of detected signals
    acqResults.codePhase    = zeros(1, 32);
    % Correlation peak ratios of the detected signals
    acqResults.peakMetric   = zeros(1, 32);
    
    fprintf('(');
    
    Nfd        = settings.acqSearchBand / 0.5 + 1;
    Sblock     = floor((samplesPerCode * settings.PIT) / Nfd);
    Nint       = Sblock * Nfd;
    Nblocks    = floor(settings.totalSamples / Sblock) - Nfd;
    
    % Perform search for all listed PRN numbers ...
    PRN = 8;
    
    signalFound = 0;
    block = 0;
    blockStart = 1;
    
    while block < Nblocks
        
        % 1. Convert signal to baseband
        
        % 2. Define block length and rearrange the signal into the matrix
        signal1   = longSignal(blockStart:blockStart+Nint-1);
        signal2   = longSignal(blockStart+Sblock:blockStart+Nint+Sblock-1);
        signal    = signal1 + signal2;
        sigMatrix = reshape(signal, Sblock, Nfd).';
        
        % 3. Define local PN code and rearrange it also into the matrix
        caCodeMatrix = repmat(caCodesTable(PRN, :), 1, ceil(Nint/samplesPerCode));
        caCodeMatrix = caCodeMatrix(1:Nint);
        caCodeMatrix = reshape(caCodeMatrix, Sblock, Nfd)';
        
        % 4. Do the circular correlation for each of the row
        corrMatrix = zeros(Nfd, Sblock);
        for row = 1:Nfd
            corrMatrix(row, :) = ifft(fft(sigMatrix(row, :)) .* conj(fft(caCodeMatrix(row, :))));
        end
        
        % 5. Do the DFT for each column and check if not greater than threshold
        corrMatrix = abs(fft(corrMatrix)).^2;
        
        % 6. If not, move one block onto next data and start again
        [~, frequencyBinIndex] = max(max(corrMatrix, [], 2));
        [maxFftAbs, codePhase] = max(max(corrMatrix));
        
        samplesPerCodeChip = round(settings.samplingFreq / settings.codeFreqBasis);
        excludeRangeIndex1 = codePhase - samplesPerCodeChip;
        excludeRangeIndex2 = codePhase + samplesPerCodeChip;
        codePhaseRange     = [1:excludeRangeIndex1, excludeRangeIndex2 : Sblock];
        secondMaxFftAbs    = max(corrMatrix(frequencyBinIndex, codePhaseRange));
        
        peakTo2ndPeakRatio         = maxFftAbs/secondMaxFftAbs;
        acqResults.peakMetric(PRN) = peakTo2ndPeakRatio;
        
        % If the result is above threshold, then there is a signal ...
        if peakTo2ndPeakRatio > settings.acqThreshold
            
            % Save properties of the detected satellite signal
            acqResults.carrFreq(PRN)  = settings.IF + 1 / (settings.PIT * 1e-3 * frequencyBinIndex);
            acqResults.codePhase(PRN) = codePhase;
            
            signalFound=1;
            break;
        end
        
        block = block + 1;
        blockStart = blockStart + Sblock;
    end % for next block in the matrix
    
    if signalFound == 1
        fprintf('%d ', PRN);
    else
        fprintf('. ');
    end
    %=== Acquisition is over ==================================================
    fprintf(')\n');

    % In this case, correlation matrix is equal to results matrix for case
    % of normal signal acquisition
    results = corrMatrix;
    frqBins = settings.IF - (settings.acqSearchBand/2) * 1000 + 0.5e3 .* (0:(Nfd-1));
    
elseif isequal(settings.acquisitionType, 'normalAcquisition')
    
    %% Initialization DDM=========================================================
    longSignal = data;
    
    samplesPerCode = round(settings.samplingFreq / ...
        (settings.codeFreqBasis / settings.codeLength));
    
    % Create two 1msec vectors of data to correlate with and one with zero DC
    signal1 = longSignal(1 : samplesPerCode);
    signal2 = longSignal(samplesPerCode+1 : 2*samplesPerCode);
    signal0DC = longSignal - mean(longSignal);
    
    % Find sampling period
    ts = 1 / settings.samplingFreq;
    % Find phase points of the local carrier wave
    phasePoints = (0 : (samplesPerCode-1)) * 2 * pi * ts;
    % Number of the frequency bins for the given acquisition band (500Hz steps)
    numberOfFrqBins = round(settings.acqSearchBand * 2) + 1;
    % Generate all C/A codes and sample them according to the sampling freq.
    caCodesTable = make_ca_table(settings);
    %--- Initialize arrays to speed up the code -------------------------------
    % Search results of all frequency bins and code shifts (for one satellite)
    results     = zeros(numberOfFrqBins, samplesPerCode);
    % Carrier frequencies of the frequency bins
    frqBins     = zeros(1, numberOfFrqBins);
    
    %--- Initialize acqResults ------------------------------------------------
    % Carrier frequencies of detected signals
    acqResults.carrFreq     = zeros(1, 32);
    % C/A code phases of detected signals
    acqResults.codePhase    = zeros(1, 32);
    % Correlation peak ratios of the detected signals
    acqResults.peakMetric   = zeros(1, 32);
    
    fprintf('(');
    
    % Perform search for all listed PRN numbers ...
    %for PRN = settings.acqSatelliteList
    
    PRN=8;
    
    %% Correlate signals ======================================================
    %--- Perform DFT of C/A code ------------------------------------------
    caCodeFreqDom = conj(fft(caCodesTable(PRN, :)));
    
    FrequencySearch1 = abs(fft(caCodesTable(PRN, :).*signal1)).^2;
    FrequencySearch2 = abs(fft(caCodesTable(PRN, :).*signal2)).^2;
    
    %--- Make the correlation for whole frequency band (for all freq. bins)
    for frqBinIndex = 1:numberOfFrqBins
        
        %--- Generate carrier wave frequency grid (0.5kHz step) -----------
        frqBins(frqBinIndex) = settings.IF - ...
            (settings.acqSearchBand/2) * 1000 + ...
            0.5e3 * (frqBinIndex - 1);
        
        %--- Generate local sine and cosine -------------------------------
        sinCarr = sin(frqBins(frqBinIndex) * phasePoints);
        cosCarr = cos(frqBins(frqBinIndex) * phasePoints);
        
        %--- "Remove carrier" from the signal -----------------------------
        I1      = sinCarr .* signal1;
        Q1      = cosCarr .* signal1;
        I2      = sinCarr .* signal2;
        Q2      = cosCarr .* signal2;
        
        %--- Convert the baseband signal to frequency domain --------------
        IQfreqDom1 = fft(I1 + 1i*Q1);
        IQfreqDom2 = fft(I2 + 1i*Q2);
        
        %--- Multiplication in the frequency domain (correlation in time
        %domain)
        convCodeIQ1 = IQfreqDom1 .* caCodeFreqDom;
        convCodeIQ2 = IQfreqDom2 .* caCodeFreqDom;
        
        %--- Perform inverse DFT and store correlation results ------------
        acqRes1 = abs(ifft(convCodeIQ1)) .^ 2;
        acqRes2 = abs(ifft(convCodeIQ2)) .^ 2;
        
        %--- Check which msec had the greater power and save that, will
        %"blend" 1st and 2nd msec but will correct data bit issues
        if (max(acqRes1) > max(acqRes2))
            results(frqBinIndex, :) = acqRes1;
        else
            results(frqBinIndex, :) = acqRes2;
        end
        
    end % frqBinIndex = 1:numberOfFrqBins
    
    %% Look for correlation peaks in the results ==============================
    % Find the highest peak and compare it to the second highest peak
    % The second peak is chosen not closer than 1 chip to the highest peak
    
    %--- Find the correlation peak and the carrier frequency --------------
    [peakSize frequencyBinIndex] = max(max(results, [], 2));
    %--- Find code phase of the same correlation peak ---------------------
    [peakSize codePhase] = max(max(results));
    
    %--- Find 1 chip wide C/A code phase exclude range around the peak ----
    samplesPerCodeChip   = round(settings.samplingFreq / settings.codeFreqBasis);
    excludeRangeIndex1 = codePhase - samplesPerCodeChip;
    excludeRangeIndex2 = codePhase + samplesPerCodeChip;
    
    %--- Correct C/A code phase exclude range if the range includes array
    %boundaries
    if excludeRangeIndex1 < 2
        codePhaseRange = excludeRangeIndex2 : ...
            (samplesPerCode + excludeRangeIndex1);
        
    elseif excludeRangeIndex2 >= samplesPerCode
        codePhaseRange = (excludeRangeIndex2 - samplesPerCode) : ...
            excludeRangeIndex1;
    else
        codePhaseRange = [1:excludeRangeIndex1, ...
            excludeRangeIndex2 : samplesPerCode];
    end
    
    %--- Find the second highest correlation peak in the same freq. bin ---
    secondPeakSize = max(results(frequencyBinIndex, codePhaseRange));
    
    %--- Store result -----------------------------------------------------
    acqResults.peakMetric(PRN) = peakSize/secondPeakSize;
    
    % If the result is above threshold, then there is a signal ...
    if (peakSize/secondPeakSize) > settings.acqThreshold
        
        %% Fine resolution frequency search =======================================
        %--- Indicate PRN number of the detected signal -------------------
        fprintf('%02d ', PRN);
        
        %--- Generate 10msec long C/A codes sequence for given PRN --------
        caCode = gen_ca_code(PRN);
        codeValueIndex = floor((ts * (1:10*samplesPerCode)) / ...
            (1/settings.codeFreqBasis));
        longCaCode = caCode((rem(codeValueIndex, 1023) + 1));
        
        %--- Remove C/A code modulation from the original signal ----------
        % (Using detected C/A code phase)
        xCarrier = ...
            signal0DC(codePhase:(codePhase + 10*samplesPerCode-1)) ...
            .* longCaCode;
        
        %--- Find the next highest power of two and increase by 8x --------
        fftNumPts = 8*(2^(nextpow2(length(xCarrier))));
        %--- Compute the magnitude of the FFT, find maximum and the
        %associated carrier frequency
        fftxc = abs(fft(xCarrier, fftNumPts));
        
        uniqFftPts = ceil((fftNumPts + 1) / 2);
        [fftMax, fftMaxIndex] = max(fftxc(5 : uniqFftPts-5));
        
        fftFreqBins = (0 : uniqFftPts-1) * settings.samplingFreq/fftNumPts;
        
        %--- Save properties of the detected satellite signal -------------
        acqResults.carrFreq(PRN)  = fftFreqBins(fftMaxIndex);
        acqResults.codePhase(PRN) = codePhase;
        %% DDM Method
        Doppler_Freq =  acqResults.carrFreq(PRN) - settings.IF;
        ZeroDelay =    acqResults.codePhase(PRN);
    else
        %--- No signal with this PRN --------------------------------------
        fprintf('. ');
    end   % if (peakSize/secondPeakSize) > settings.acqThreshold
end


figure(1023)
subplot(2,1,1)
mesh(results)
xlabel('Code phase [chips]')
ylabel('Frequency [MHz]')
title('correlation Map')

subplot(2,1,2)
imagesc(results)
colorbar

figddm = figure(1024);
subplot(2,1,1)

dopplerAxis = frequencyBinIndex-20:frequencyBinIndex+20;
dopplerAxis(find(dopplerAxis<1)) = [];
if isequal(settings.acquisitionType, 'weakAcquisition')
    dopplerAxis(find(dopplerAxis>Nfd)) = [];
else
    dopplerAxis(find(dopplerAxis>numberOfFrqBins)) = [];
end

delayAxis = codePhase-50:codePhase+100;
if isequal(settings.acquisitionType, 'weakAcquisition')
    delayAxis(find(delayAxis > Sblock)) = [];
else
    delayAxis(find(delayAxis > samplesPerCode)) = [];
end

DopplerY = frqBins(dopplerAxis)-acqResults.carrFreq(PRN);
TotalDelaysec = ((1:samplesPerCode)-codePhase)*ts*1e9;
DelayX = TotalDelaysec(delayAxis);
[Xax,Yax]=meshgrid(DelayX,DopplerY);
DDM = results(dopplerAxis,delayAxis);


mesh(Xax,Yax,DDM)
xlabel('Delay axis [ns]')
ylabel('Doppler axis [Hz]')
title('Delay Doppler Map')
box on

subplot(2,1,2)
imagesc(DelayX,DopplerY,DDM)
colorbar
xlabel('Delay axis [ns]')
ylabel('Doppler axis [Hz]')
title('Delay Doppler Map')

%     filename = strsplit(settings.fileName,'/');
%Savefilename = strsplit(filename{end},'.bin');
%saveas(figddm,sprintf('%s.jpg',Savefilename{1}))
%end    % for PRN = satelliteList
%=== Acquisition is over ==================================================
fprintf(')\n');
