function plot_acquisition_results(sdrParams, postProcessResults)
% SUMMAR for plotting different plots for comparison
% of performance, code profiling and accuracy of results.


%%% Extract coherent frame data.
currFileNum     = sdrParams.stateParams.numFilesProcessed + 1;
dataParams      = sdrParams.dataFileParamsList{currFileNum};

channelId = 1;
if dataParams.selectedChannel ~= -1
    numChannels = 1;
    channelId = dataParams.selectedChannel;
else
    numChannels = dataParams.totalChannels;
    channelId = 1:dataParams.totalChannels;
end
                 
acqAlgoList = sdrParams.sysParams.acqAlgosList;
numAcqAlgos = length(acqAlgoList);

for chIdx=1:numChannels
    for acqAlgoIdx=1:numAcqAlgos
        
        acqResults = postProcessResults.acqResults{chIdx, acqAlgoIdx}.acqResults;        
         if ~isempty(acqResults)
             
             print_string(['Plotting bar graph fpr ', num2str(chIdx),...
                             '/', num2str(numChannels), ' for algorithm ', acqAlgoList{acqAlgoIdx}]);

             % Bar plot for peak metric bar graph of acquired satellites
             plot_acq_peak_metric_bars(sdrParams, postProcessResults, chIdx, channelId, acqAlgoIdx);
             
             % DDM plots
             plot_acq_ddm_map(sdrParams, postProcessResults, chIdx, channelId, acqAlgoIdx);
             
         end
    end
end
end

function plot_acq_peak_metric_bars(sdrParams, postProcessResults, chIdx, channelId, algoIdx)

acqResults        = postProcessResults.acqResults{chIdx, algoIdx}.acqResults.chAlgoAcqResults;
nacqPrnPeakMetric = postProcessResults.acqResults{chIdx, algoIdx}.acqResults.nacqStats.nacqPrnPeakMetric;
nacqPrnList       = postProcessResults.acqResults{chIdx, algoIdx}.acqResults.nacqStats.nacqPrnList;
numAcqSatellites  = length(acqResults);
currFileNum       = sdrParams.stateParams.numFilesProcessed + 1;
currFrameNum      = sdrParams.stateParams.currFrameNum;
numTotalFrames    = sdrParams.stateParams.numTotalFrames;
acqAlgoList = sdrParams.sysParams.acqAlgosList;


% prepare fields
titleStr = sprintf('Acquisition Results. Total Frames %d\n', numTotalFrames);
channelStr = sprintf('Data channel: %d\n', channelId(chIdx));
acqAlgoStr = sprintf('Acquisition Algorithm: %s', acqAlgoList{algoIdx});

acqPeakMetric = [];
acqPeakMetricLoc = [];
for aprnIdx=1:numAcqSatellites
    acqPeakMetric = [acqPeakMetric, acqResults{aprnIdx}.peakMetric];
    acqPeakMetricLoc = [acqPeakMetricLoc, acqResults{aprnIdx}.satellitePrn];
end
peakMetric = [acqPeakMetric, nacqPrnPeakMetric];
peakMetricLoc = [acqPeakMetricLoc, nacqPrnList];

figWindowTitle = sprintf('Acquisition Bar Graph. File#%d, Data_channel#%d, Algo#%d, Frame#%d', ...
    currFileNum, chIdx, algoIdx, currFrameNum);
figure('NumberTitle', 'off', 'Name', figWindowTitle);

ylimit = max(acqPeakMetric)*1.5;
hAxes = newplot();
hAxes.NextPlot = 'replaceall';
bar(hAxes, peakMetricLoc, peakMetric);
title (hAxes, [titleStr, channelStr, acqAlgoStr], 'Interpreter', 'none');
xlabel(hAxes, 'PRN number (no bar - SV is not in the acquisition list)');
ylabel(hAxes, 'Acquisition Metric');
hAxes.YLim = [0, ylimit]; % Set ylimit to 20 for easier eyeball comparison of plots.

oldAxis = axis(hAxes);
axis  (hAxes, [0, 33, 0, oldAxis(4)]);
set   (hAxes, 'XMinorTick', 'on');
set   (hAxes, 'YGrid', 'on');

% Mark acquired signals
qualifiedSignals = peakMetric .* (peakMetric > sdrParams.sysParams.acqThreshold);
hold(hAxes, 'on');
bar (hAxes, peakMetricLoc, qualifiedSignals, 'FaceColor', [0 0.8 0]);

acquiredSignals = peakMetric .* ismember(peakMetric, acqPeakMetric);
bar (hAxes, peakMetricLoc, acquiredSignals, ...
    'FaceColor',[0 0.8 0],'EdgeColor',[0 .9 .9],'LineWidth',1.5);
legend(hAxes, ...
    'Not acquired signals', ...
    'Acquired signals', ...
    'Selected signals', ...
    'AutoUpdate', 'off');

axisPos = hAxes.Position;
leftBottomPointX = axisPos(1);
leftBottomPointY = axisPos(2);
rightTopX = axisPos(3)+axisPos(1);
rightTopY = axisPos(4)+axisPos(2);
plotWidth = rightTopX;
plotHeight = rightTopY;

x = [leftBottomPointX + plotWidth  * 0.8, leftBottomPointX + plotWidth  * 0.8 ];
y = [leftBottomPointY + plotHeight * 0.6, leftBottomPointY + ...
    (sdrParams.sysParams.acqThreshold/ylimit)*axisPos(4)];
annotation('textarrow',x,y,'String','Threshold ', 'Units','normalized', 'Color', 'r');
plot(0:33, sdrParams.sysParams.acqThreshold*ones(1, 34), 'r--');
hold(hAxes, 'off');
end


function plot_acq_ddm_map(sdrParams, postProcessResults, chIdx, channelId, algoIdx)

acqResults        = postProcessResults.acqResults{chIdx, algoIdx}.acqResults.chAlgoAcqResults;
numAcqSatellites  = length(acqResults);
currFileNum       = sdrParams.stateParams.numFilesProcessed + 1;
numTotalFrames    = sdrParams.stateParams.numTotalFrames;
acqAlgoList = sdrParams.sysParams.acqAlgosList;


% prepare fields

for aprnIdx=1:numAcqSatellites
    
    figWindowTitle = sprintf('Acquisition Delay Doppler Map. File#%d, Data_channel#%d, Algo#%d, Frames#%d', ...
    currFileNum, chIdx, algoIdx, numTotalFrames);

    figure('NumberTitle', 'off', 'Name', figWindowTitle);
    
    dopplerShiftAxis = acqResults{aprnIdx}.dopplerShiftAxis;
    codeDelayAxis = acqResults{aprnIdx}.codeDelayAxis;

    
    mesh(codeDelayAxis, dopplerShiftAxis, acqResults{aprnIdx}.ddm);


    titleStr = sprintf('Acquisition Delay Doppler Map\n ');
    channelStr = sprintf('Data channel: %d\n', channelId(chIdx));
    acqAlgoStr = sprintf('Acquisition Algorithm: %s\n', acqAlgoList{algoIdx});
    prnStr     = sprintf('PRN : %d', acqResults{aprnIdx}.satellitePrn);

    
    title ([titleStr, channelStr, acqAlgoStr, prnStr], 'Interpreter', 'none');
    zlabel('Correlation Magnitude');
    ylabel('Doppler Shift (Hz)');
    xlabel('Code Phase (samples)');

end

end