function probe_data(varargin)


%%% Check the number of arguments 
if (nargin == 1)
    settings = deal(varargin{1});
    fileNameStr = settings.fileName;
elseif (nargin == 2)
    [fileNameStr, settings] = deal(varargin{1:2});
    if ~ischar(fileNameStr)
        error('File name must be a string');
    end
else
    error('Incorect number of arguments');
end
    
%%% Generate plot of raw data 
[fid, message] = fopen(fileNameStr, 'rb');
settings.fileID = fid;
if (fid > 0)
    % Move the starting point of processing. Can be used to start the
    % signal processing at any point in the data record (e.g. for long
    % records).
    %fseek(fid, settings.skipNumberOfBytes, 'bof');    
     fseek(fid,0, 'bof');    
    
    % Read 10ms of signal
    [tmpdata, count] = fread(fid, [1, 10*settings.samplesPerCode], settings.dataType);
    fclose(fid);
    
    % convert the signal to symbols from bits
    % TODO: this moves to pre-processing for
    
    channelNum = 4;
    numChannels = 4;
    tmpdata = int32(tmpdata(channelNum:numChannels:end));
    data  = double(1 + 2 * idivide(tmpdata, 2) .* (1 - 2 * mod(tmpdata, 2)));    
    
    if (count < 10*settings.samplesPerCode)
        % The file is to short
        error('Could not read enough data from the data file.');
    end
    
    %%% Initialization 
    figure(100);
    clf(100);
    timeScale = 0 : 1/settings.samplingFreq : 5e-3;
    
    %%% Time domain plot
    subplot(2, 2, 1); plot(1000 * timeScale(1:round(settings.samplesPerCode/50)), ...
                           data(1:round(settings.samplesPerCode/50)));
    axis tight; grid on; 
    title ('Time domain plot');
    xlabel('Time (ms)'); 
    ylabel('Amplitude');
    ylim([-5 5])
    
    %%% Frequency domain plot
    subplot(2,2,2); pwelch(data-mean(data), 16384, ...
                           1024, 2048, settings.samplingFreq/1e6,'centered');
    
    axis tight; grid on;
    title ('Frequency domain plot');
    xlabel('Frequency (MHz)'); ylabel('Magnitude');
    
    %%% Histogram
    subplot(2, 2, 3.5); hist(data, -128:128)
    dmax = max(abs(data)) + 1;
    axis tight;
    adata = axis;
    axis([-dmax dmax adata(3) adata(4)]);
    grid on;
    title ('Histogram'); 
    xlabel('Bin'); ylabel('Number in bin');
else
    %%% Error while opening the data file 
    error('Unable to read file %s: %s.', fileNameStr, message);
end % if (fid > 0)
