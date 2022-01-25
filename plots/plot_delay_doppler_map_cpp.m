
clear; %close all;
% For C++ 
fileName = ['NTLab_Bands_GPS_GLONASS_L12'];
load(['../GNSS_SDR_C++/',fileName,'_ddm_mat.mat']);
load(['../GNSS_SDR_C++/',fileName,'_cpp.mat']);

%For matlab
% load('results_sim0_18.mat');
% load('acqResults_sim0.mat');

PRN = 8;
frequencyBinIndex = frequencyBinIndex + 1;
codePhaseLast = codePhaseLast + 1;
results = reshape(results, samplesPerCode, numberOfFrqBins)';

ts = 1 / samplingFreq;
figure(1025)
    subplot(2,1,1)
    mesh(results)
    xlabel('Code phase [chips]')
    ylabel('Frequency [MHz]')
    title('correlation Map')
    
    subplot(2,1,2)
    imagesc(results)
    colorbar
    
    figddm = figure(1026);
    subplot(2,1,1)

    dopplerAxis = frequencyBinIndex-20:frequencyBinIndex+20;
    dopplerAxis(find(dopplerAxis < 1)) = [];
    dopplerAxis(find(dopplerAxis > numberOfFrqBins)) = [];
    
    delayAxis = codePhaseLast-50:codePhaseLast+100;
    delayAxis(find(delayAxis < 1)) = [];
    delayAxis(find(delayAxis > samplesPerCode)) = [];
    
    DopplerY = frqBins(dopplerAxis)-carrFreq(PRN);
    TotalDelaysec = ((1:samplesPerCode)-codePhaseLast)*ts*1e9;
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
    title('Delay Doppler Map');
