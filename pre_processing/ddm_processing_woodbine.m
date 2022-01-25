% clear; close all; 
format ('compact');
format ('long', 'g');
%--------------------------------------------------------------
addpath include             % The software receiver functions
addpath geoFunctions        % Position calculation related functions
%--------------------------------------------------------------
GPSL1freq = 1575.42*power(10,6);

Channel = 1;

switch Channel
    case 1
%         NLHCP = load('./../../Woodbine_47a_ch1_1.mat');
        SourceFilepath = sprintf('./../../Woodbine_47a_ch1_1.bin');
    case 2
        NRHCP = load('./../../Woodbine_47a_ch2_1.mat');
        SourceFilepath = sprintf('./../../Woodbine_47a_ch2_1.bin');
    case 3
        ZRHCP = load('./../../Woodbine_cf7_ch1_1.mat');
        SourceFilepath = sprintf('./../../Woodbine_cf7_ch1_1.bin');
end

[~,SourceFilename,ext] = fileparts(SourceFilepath);
settings = initSettings(SourceFilepath);
samplesPerCode = round(settings.samplingFreq/(settings.codeFreqBasis / settings.codeLength));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fid, message] = fopen(settings.fileName, 'r','ieee-le');
fseek(fid, settings.skipNumberOfBytes, 'bof');    
[Sample, count] = fread(fid,[1, samplesPerCode*1000*2], settings.dataType);
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I = Sample(1:2:end)/2048;
Q = Sample(2:2:end)/2048;
%--------------------------------------------------------------
ts = 1/settings.samplingFreq;   % Sampling period in sec
tc = 1/settings.codeFreqBasis;  % C/A chip period in sec
%--------------------------------------------------------------
DelayDopplerMap = zeros(11,21);
DDMcollector = cell(100,32);
SatelliteList = 1:32;
%--------------------------------------------------------------
Delayaxis = -5:1:5;
Doppleraxis = -5000:500:5000;
DDMPhase = (0 : (samplesPerCode-1)) * 2 * pi * ts;

k=2;
%for k=1:10

switch Channel
    case 1
        PRN = SatelliteList(NLHCP.acqResults(k).peakMetric>1.3);
        DopplerFreq = NLHCP.acqResults(k).carrFreq(PRN);
        Delay = NLHCP.acqResults(k).codePhase(PRN);        
    case 2
        PRN = SatelliteList(NRHCP.acqResults(k).peakMetric>1.3);
        DopplerFreq = NRHCP.acqResults(k).carrFreq(PRN);
        Delay = NRHCP.acqResults(k).codePhase(PRN);        
    case 3
        PRN = SatelliteList(ZRHCP.acqResults(k).peakMetric>1.3);
        DopplerFreq = ZRHCP.acqResults(k).carrFreq(PRN);
        Delay = ZRHCP.acqResults(k).codePhase(PRN);
end

    %--------------------------------------------------------------
    signal_i = I(samplesPerCode*(k-1)+1:samplesPerCode*k);
    signal_q = Q(samplesPerCode*(k-1)+1:samplesPerCode*k);
    RecSignal = signal_i+1i*signal_q;
    %RecSignal = signal_i-signal_q;
    %--------------------------------------------------------------
    tmpDDMcollector = cell(1,32);

    for ll = 1:length(PRN)
        CAcode = generateCAcode(PRN(ll));
        codeValueIndex = ceil((ts * (1:samplesPerCode)) / tc);
        codeValueIndex(end) = 1023;
        tmpCACodesTable= CAcode(codeValueIndex);

        for ii = 1:length(Delayaxis)
            for jj = 1:21
                CACodesTable = circshift(tmpCACodesTable,-(Delay(ll)-39*Delayaxis(ii)));
                %Y = sum(RecSignal.*CACodesTable.*(cos(DDMPhase.*(DopplerFreq(ll)-Doppleraxis(jj)))+1i*sin(DDMPhase.*(DopplerFreq(ll)-Doppleraxis(jj)))),'omitnan');
                YY = sum(RecSignal.*CACodesTable.*(exp(-1i*DDMPhase.*(GPSL1freq+Doppleraxis(jj)))) ,'omitnan');
                DelayDopplerMap(ii,jj) = (abs(YY))^2;
            end
        end

        tmpDDMcollector{PRN(ll)}=DelayDopplerMap;

    end

    DDMcollector(k,:) = tmpDDMcollector;
    
    fprintf('%d th acquistion results has processed \n',k);
%end

%save(sprintf('/Volumes/TOSHIBA/bladeRF/200223/1/Results/cf7_ch1/DDM_ZRHCP_1_sub.mat'),'DDMcollector');
    
max(max(DDMcollector{2, 3}))

% figure(1)
% imagesc(Doppleraxis,Delayaxis,DelayDopplerMap)
% xlabel('\DeltaDoppler Frequency')
% ylabel('Delay (chip)')
