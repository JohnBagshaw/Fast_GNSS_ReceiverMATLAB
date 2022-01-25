% Script postProcessing.m processes the raw signal from the specified data
% file (in settings) operating on blocks of 37 seconds of data.
%
% First it runs acquisition code identifying the satellites in the file,
% then the code and carrier for each of the satellites are tracked, storing
% the 1msec accumulations.  After processing all satellites in the 37 sec
% data block, then postNavigation is called. It calculates pseudoranges
% and attempts a position solutions. At the end plots are made for that
% block of data.

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Darius Plausinaitis
% Written by Darius Plausinaitis, Dennis M. Akos
% Some ideas by Dennis M. Akos
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

%                         THE SCRIPT "RECIPE"
%
% The purpose of this script is to combine all parts of the software
% receiver.
%
% 1.1) Open the data file for the processing and seek to desired point.
%
% 2.1) Acquire satellites
%
% 3.1) Initialize channels (preRun.m).
% 3.2) Pass the channel structure and the file identifier to the tracking
% function. It will read and process the data. The tracking results are
% stored in the trackResults structure. The results can be accessed this
% way (the results are stored each millisecond):
% trackResults(channelNumber).XXX(fromMillisecond : toMillisecond), where
% XXX is a field name of the result (e.g. I_P, codePhase etc.)
%
% 4) Pass tracking results to the navigation solution function. It will
% decode navigation messages, find satellite positions, measure
% pseudoranges and find receiver position.
%
% 5) Plot the results.

%% Initialization =========================================================
disp ('Starting processing...');

[fid, message] = fopen(settings.fileName, 'rb', 'ieee-le');

%If success, then process the data
if (fid > 0)
    
    % Move the starting point of processing. Can be used to start the
    % signal processing at any point in the data record (e.g. good for long
    % records or for signal processing in blocks).
    fseek(fid, settings.skipNumberOfBytes, 'bof');
    % Do acquisition if it is not disabled in settings or if the variable
    % acqResults does not exist.
    
    if ((settings.skipAcquisition == 0) || ~exist('acqResults', 'var'))

        % Read data for acquisition.
        fseek(fid, settings.skipNumberOfBytes, 'bof');
        tmpdata = fread(fid, settings.dataExtractLen, settings.dataType)';

        tmpdata(tmpdata==1) = -1;
        tmpdata(tmpdata==3) = -3;
        tmpdata(tmpdata==0) = 1;
        tmpdata(tmpdata==2) = 3;

        L = floor((length(tmpdata)/4));

        Data4I = double(tmpdata(1:4:4*(L)));
        Data3I = double(tmpdata(2:4:4*(L)));
        Data2I = double(tmpdata(3:4:4*(L)));
        Data1I = double(tmpdata(4:4:4*(L)));

        data = Data1I;

        if settings.modulationRequired
            dataLen     = settings.dataExtractLen / 2;
            phasePoints = 2 * pi * settings.IF / settings.samplingFreq * (0 : dataLen-1);
            data        = data(1:2:end) .* cos(phasePoints) - data(2:2:end) .* sin(phasePoints);
        elseif settings.demodulationRequired
            phasePoints = (0 : settings.dataExtractLen/4-1) * 2 * pi * settings.IF / settings.samplingFreq;
            data        = data .* cos(phasePoints) + 1i * data .* sin(phasePoints);
        end
        
        disp ('   Acquiring satellites...');
        
        if isequal(settings.acquisitionType, 'normalAcquisition')
            acqResults = acquisition(data, settings);
        elseif isequal(settings.acquisitionType, 'weakAcquisition')
            acqResults = weak_acquisition(data, settings);
        elseif isequal(settings.acquisitionType, 'halfBitAcquisition')
            acqResults = half_bit_acquisition(data, settings, settings.dataExtractLen);
        end
        
    end
    
    plot_acquisition(acqResults);

%% Initialize channels and prepare for the run ============================
    % Start further processing only if a GNSS signal was acquired (the
    % field FREQUENCY will be set to 0 for all not acquired signals)
    if (any(acqResults.carrFreq))
        channel = pre_run(acqResults, settings);
        show_channel_status(channel, settings);
    else
        % No satellites to track, exit
        disp('No GNSS signals detected, signal processing finished.');
        trackResults = [];
        return;
    end

    % Close the data file
    fclose(fid);
else
    % Error while opening the data file.
    error('Unable to read file %s: %s.', settings.fileName, message);
end % if (fid > 0)
