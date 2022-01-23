%%
 % Project Title: GNSS-R SDR
 % Author       : John Bagshaw
 % Contact      : jotshaw@yorku.ca
 % Supervisor   : Prof.Sunil Bisnath
 % Institution  : York University, Canada.
%%

function sdrParams = config_sdr_params()
%%% This function is responsible for creating global
% and per file settings. This is run once per data
% file, stream and 

print_string("Initializing SDR parameters.");

%%% Receiver system parameters

% Users can specify most important parameters here.
% TODO: Move these important user settings to settings.txt input file.
sdrParams.sysParams.acqAlgosList = {'norm_acq_parcode', ...
                                    'weak_acq_dbzp',...
                                   %  'weak_acq_dbzp'   ,...
                                   % Add more... e.g.
                                   % 'weak_acquisition_hb',...
                                  };
    
sdrParams.sysParams.coherentProcessingTimeMS   = 1;
sdrParams.sysParams.incoherentProcessingTimeMS = 2;

if sdrParams.sysParams.coherentProcessingTimeMS > 10
    sdrParams.sysParams.coherentProcessingTimeMS = 10; % Cannot be greater than 10ms.
end

if sdrParams.sysParams.incoherentProcessingTimeMS < sdrParams.sysParams.coherentProcessingTimeMS
    error("Incoherent processing time must be greater than coherent processing time.");
end

% System parameters
sdrParams.sysParams.c             = 299792458; % The speed of light, [m/s]
sdrParams.sysParams.startOffset   = 68.025;    %[ms] Initial sign. travel time
sdrParams.sysParams.caCodeChipRateHz = 1.023e6;  % Code frequency basis [Hz]
sdrParams.sysParams.minSamplingFreqHz = sdrParams.sysParams.caCodeChipRateHz * 2;
sdrParams.sysParams.sampleInterpOrder = 2;

% Acquisition Parameters
sdrParams.sysParams.numberOfChannels  = 8;       % Number of channels to be used for signal processing
sdrParams.sysParams.skipNumberOfBytes = 0;       % Move the starting point of processing.
sdrParams.sysParams.skipAcquisition   = 0;       % Skips acquisition in the script postProcessing.m if set to 1
sdrParams.sysParams.acqSatelliteList  = 1:32;    % List of satellites to look for. Some satellites can be excluded to speed up acquisition [PRN numbers]
sdrParams.sysParams.acqDopplerBwKhz   = 50;      % Band around IF to search for satellite signal. Depends on max Doppler
sdrParams.sysParams.acqDopplerResHz   = 1000;    % Doppler resolution
sdrParams.sysParams.acqThreshold      = 2.5;     % Threshold for the signal presence decision rule
sdrParams.sysParams.numAcqSatellites  = 8;       % Keep record of strongest 8 satellites.

% Tracking Parameters
sdrParams.sysParams.msToProcess          = 36;  % [ms]
sdrParams.sysParams.dllDampingRatio      = 0.7; % Code tracking loop parameters
sdrParams.sysParams.dllNoiseBandwidth    = 2;   %[Hz]
sdrParams.sysParams.dllCorrelatorSpacing = 0.5; %[chips]
sdrParams.sysParams.pllDampingRatio      = 0.7; % Carrier tracking loop parameters
sdrParams.sysParams.pllNoiseBandwidth    = 25;  %[Hz]

%%% Navigation solution settings
sdrParams.sysParams.navSolPeriod       = 500; % Period for calculating pseudoranges and position [ms]
sdrParams.sysParams.elevationMask      = 10;  % Elevation mask to exclude signals from satellites at low elevation[degrees 0 - 90]
sdrParams.sysParams.useTropCorr        = 1;   % % Enable/dissable use of tropospheric correction [0 - Off, 1 - On]
sdrParams.sysParams.truePosition.E     = nan; % True position of the antenna in UTM system (if known). Otherwise enter
sdrParams.sysParams.truePosition.N     = nan; % all NaN's and mean position will be used as a reference .
sdrParams.sysParams.truePosition.U     = nan;

%%% Data parameters
sdrParams.stateParams.dataPathIn        = '..\data\data_in\';
sdrParams.stateParams.dataPathOut       = '..\data\data_out\';
dataInfoFileName  = 'data_file_list.txt';
if ~(exist([sdrParams.stateParams.dataPathIn, ...
        dataInfoFileName], 'file') == 2)
    error(['Data list file is not present. Please define a ',...
        'data_file_list.txt file in ~\..\data folder']);
end

fileNameStr = fileread([sdrParams.stateParams.dataPathIn, ...
    dataInfoFileName]); % Read data file names form data list file.
sdrParams.stateParams.fileNames = regexp(fileNameStr, '\r\n|\r|\n', 'split');

for fileName=sdrParams.stateParams.fileNames
    
    % Check if file exists in the path.
    if exist([sdrParams.stateParams.dataPathIn, fileName{1}], 'file') == 2
        
        %sdrParams.dataParams.fileDataParams;
        dataParams = [];
        % TODO: Convert formats to JSON.
        switch fileName{1}
            case 'Woodbine_47a'
                dataParams.intermFreqHz     = 10e6;    % Intermediate frequency
                dataParams.samplingFreqHz   = 40e6;    % Sampling frequency
                dataParams.dataType         = 'int16'; % Data type used to store one sample
                dataParams.isBasebandSignal = 1;
                dataParams.totalChannels    = 1;
                dataParams.selectedChannel  = 1;
                
            case 'GPS_and_GIOVE_A-NN-fs16_3676-if4_1304.bin'
                dataParams.intermFreqHz     = 4.1304e6;  % Intermediate frequency
                dataParams.samplingFreqHz   = 16.367e6; % Sampling frequency
                dataParams.dataType         = 'int8';    % Data type used to store one sample
                dataParams.isBasebandSignal = 0;
                dataParams.totalChannels    = 1;
                dataParams.selectedChannel  = 1;
                
            case 'GPSdata-DiscreteComponents-fs38_192-if9_55.bin'
                dataParams.intermFreqHz     = 9.548e6;  % Intermediate frequency
                dataParams.samplingFreqHz   = 38.192e6; % Sampling frequency
                dataParams.dataType         = 'int8';   % Data type used to store one sample
                dataParams.isBasebandSignal = 0;
                dataParams.totalChannels    = 1;
                dataParams.selectedChannel  = 1;
                
            case 'NTLab_Bands_GPS_GLONASS_L12.bin'
                dataParams.intermFreqHz     = (1590-1575.42)*1e6; % Intermediate frequency
                dataParams.samplingFreqHz   = 53e6;               % Sampling frequency
                dataParams.dataType         = 'ubit2';            % Data type used to store one sample
                dataParams.isBasebandSignal = 0;
                dataParams.totalChannels    = 4;
                dataParams.selectedChannel  = 4; % set to -1 for all channels
                
              case 'Woodbine_47a_1.bin'
                dataParams.intermFreqHz     = 9.548e6;  % Intermediate frequency
                dataParams.samplingFreqHz   = 38.192e6; % Sampling frequency
                dataParams.dataType         = 'int8';   % Data type used to store one sample
                dataParams.isBasebandSignal = 0;
                dataParams.totalChannels    = 1;
                dataParams.selectedChannel  = 1;
                
              case 'Woodbine_47a_2.bin'
                dataParams.intermFreqHz     = 9.548e6;  % Intermediate frequency
                dataParams.samplingFreqHz   = 38.192e6; % Sampling frequency
                dataParams.dataType         = 'int8';   % Data type used to store one sample
                dataParams.isBasebandSignal = 0;
                dataParams.totalChannels    = 1;
                dataParams.selectedChannel  = 1;
                
            otherwise
                error("Listed file is not recognized by SDR. If newer, please add to init_settings or json format.");
        end
        
        if isfield(sdrParams, 'dataFileParamsList')
            sdrParams.dataFileParamsList{end+1, :} = dataParams;
        else
            sdrParams.dataFileParamsList = {dataParams};
        end
        
    else
        error([fileName, ' is not present in ', ...
            sdrParams.sysParams.dataPathIn, '.']);
    end
end

%%% State Params
sdrParams.stateParams.numFilesProcessed = 0; % Set processed files to zero
sdrParams.stateParams.numFilesToProcess = length(sdrParams.stateParams.fileNames); % number of total data sets to process

sdrParams.stateParams.currFrameNum  = 0;
sdrParams.stateParams.numTotalFrames = sdrParams.sysParams.incoherentProcessingTimeMS / ...
                                       sdrParams.sysParams.coherentProcessingTimeMS;

if floor(sdrParams.stateParams.numTotalFrames) ~= sdrParams.stateParams.numTotalFrames
    error('Number of incoherent frames of data not integer');
end




%%% end.
