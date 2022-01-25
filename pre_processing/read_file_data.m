function [ outData ] = read_file_data( sdrParams )
%%%READ_FILE_DATA Read input file data and converts to suitable format.
%   Reads input file data and returns a cell of dimension
% of  1xM is number of channels in the data file.


% Extract relevant parameters
fileNum = sdrParams.stateParams.numFilesProcessed + 1;
fileName = sdrParams.stateParams.fileNames{fileNum};
dataType = sdrParams.dataFileParamsList{fileNum}.dataType;
samplingFreqHz = sdrParams.dataFileParamsList{fileNum}.samplingFreqHz;
totalChannels = sdrParams.dataFileParamsList{fileNum}.totalChannels;
selectedChannel = sdrParams.dataFileParamsList{fileNum}.selectedChannel;
skipNumberOfBytes = sdrParams.sysParams.skipNumberOfBytes;

% Define data buffer
outData = cell(1);

% Open file.
currFileFullName = [sdrParams.stateParams.dataPathIn, fileName];
[fid, message] = fopen(currFileFullName, 'rb', 'ieee-le');
if fid <= 0
    % Error while opening the data file.
    error('Unable to read file %s: %s.', stateParams.currFileName, message);
else
    
    % File opened successfully, read the contents.
    processDataDurMs  = max([sdrParams.sysParams.coherentProcessingTimeMS, ...
        sdrParams.sysParams.incoherentProcessingTimeMS]);
    processDataDurSec = processDataDurMs  * 1e-3;
    numSamples = floor(processDataDurSec * samplingFreqHz) * totalChannels;
    
    % Skip initial bytes
    fseek(fid, skipNumberOfBytes, 'bof');
    
    % Read from data file.
    rawData = fread(fid, numSamples, dataType)';
    
    % Check if sufficient data is present in the file.
    if length(rawData) < numSamples
        error('Configured processing interval data is not present in data file.');
    end
    
    % Convert from bits to symbols if necessary
    % This is file specific.
    switch fileName
        case 'NTLab_Bands_GPS_GLONASS_L12.bin'
            mOneIdx   = rawData == 1;
            mThreeIdx = rawData == 3;
            pOneIdx   = rawData == 0;
            pThreeIdx = rawData == 2;
            
            rawData(pOneIdx)   = 1;
            rawData(mOneIdx)   = -1;
            rawData(pThreeIdx) = 3;
            rawData(mThreeIdx) = -3;
            
            if selectedChannel == -1
                % Read all data from all channel
                for chIdx=1:totalChannels
                        print_string(['Reading input data for channel : ', ...
                            num2str(chIdx), '/', ... 
                            num2str(totalChannels)]);
                    outData{chIdx, :} = rawData(chIdx:totalChannels:end);
                end
            else
                % Read all data from selected channel
                print_string(['Reading input data for channel : ', ...
                    num2str(selectedChannel), '/', ... 
                    num2str(selectedChannel)]);
                outData{1} = rawData(selectedChannel:totalChannels:end);
            end
            
        case 'GPS_and_GIOVE_A-NN-fs16_3676-if4_1304.bin'
            
        otherwise
            print_string('File data is assumed to be in symbols not bits.');
            
    end
end
end
    

