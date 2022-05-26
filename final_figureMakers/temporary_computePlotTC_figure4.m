%% New tuning curves
% remake the entire calculation, not based on single runs, but on the whole cluster.
% pipeline:
% for each cluster
% 1. store ext baseline points -> average -> calculate std (and store it)
% 2. store the fucking response points.
% 3. use a fucking mov average to find 5-points peak -> average (and store) 
% 
% 
% sizeMap = [60, 86];
% clust = remapIndexedPointsToRun_RoisMaps(Tparent, pixel2keep, sizeMap);
% ROIs = clust.ROIs; %these are still unsorted, so do not remove the last one!
% 
% trialTypes = {chunks.pips.amplitude_1, chunks.pips.amplitude_2, chunks.pips.amplitude_3, ...
%               chunks.tones.amplitude_1_carrPhase_1, chunks.tones.amplitude_2_carrPhase_1 }; %keep separate subsets
% saveNameFTData = fullfile(clusterFolder, 'FTdata.mat');
ts_dec = T(3).ts_dec;
saveName = fullfile(Folder2Save, 'R_chunks.mat');


freqTuning_pipsTones_extendedBaselines_extractTraces(trialTypes, aZ, saveNameFTData); % only run if not already.
movingwindowsize = 1;   %for moving average
responseWindow = 1;     %seconds after stim onset. Empty [] for entire default window
if ~isempty(responseWindow)
    responseWindowPoints = sum(ts_dec>=0 & ts_dec <=responseWindow);
else
    responseWindowPoints = [];
end

%this is currently hijacked in script:
fTC = freqTuning_pipsTones_extendedBaselines(trialTypes, ROIs, movingwindowsize, saveNameFTData, responseWindowPoints); % baseline+ext baseline: 102 points
% load('fTC_data_1secResponse.mat') 
% load('fTC_data_02_07secAfterOnset.mat')
% load('fTC_data_01_06secAfterOnset.mat')

% calcluate perc dff of peak response
for i = 1 : length(trialTypes)
    for k = 1 : size(ROIs,2)
        fTC(i).dffPeak(:,:,k) = ( (max(fTC(i).movAvgResponse(:,:,k)) - fTC(i).meanNoise(:,:,k)) ./ fTC(i).meanNoise(:,:,k) ).* 100;
    end
end

%% find max dff values for each cluster from figure 2A
maxAllStimuli = [];
negMaxAllStimuli = [];
m = matfile(saveName);
for i = 1 : 3
    switch i
        case 1
            data = m.R_chirps;
        case 2
            data = m.R_songs;
        case 3
            data = m.R_steps;
    end
    for k = 1 : size(ROIs,2)
        maxAllStimuli(k,i) = max(mean(data(klust(k).k,:)));
        negMaxAllStimuli(k,i) = min(mean(data(klust(k).k,:)));
    end
end
% sanity check\ %plot step responses
% [hfig, hax] = figureTracesI_PP_pixels( 11, 200, 30 );
% hax = flipud(hax); % 1 is bottom, as indendrogram
% for Ki = 1:11
%     axes(hax(Ki))
%     plot(mean(data(klust(Ki).k,:)))
% end
%OK
% incorporate peaks from pips and tones
for i = 1:5
    for k = 1 : size(ROIs,2)
        maxAllStimuli(k,i+3) = max(fTC(i).dffPeak(:,:,k));
        negMaxAllStimuli(k,i+3) = min(fTC(i).dffPeak(:,:,k));
    end
end

dffLimits = [min(negMaxAllStimuli,[], 2), max(maxAllStimuli, [],2)]; %this takes K2use

%% rescale by baseline noise
% recalculate a single baseline noise level for each cluster
for k = 1 : size(ROIs,2)
    basTemp = [];
    for i = 1:5
        basTemp = cat(1, basTemp, reshape(fTC(i).baselines(:,:,k), [],1,1));
    end
    fullBaselPipsTones(:, k) = basTemp;
end
noiseClusters = std(fullBaselPipsTones); %ok.
scaledNoiseClusters = ( (noiseClusters-min(noiseClusters)) ./ max((noiseClusters-min(noiseClusters))) ) / 2 ; %scaled between 0 (min noise) and 0.5 (max noise). [To be subtracted]
% scaledNoiseClusters =
%     0.0434    0.0339    0.0133    0.5000    0.1400    0.1466    0.0201    0.0395    0.2126    0.0610         0

%% single curve with 4Hz from tones and rest from pips, max amplitude only
[xvar, ia, ic] = unique([chunks.pips.carrierFreqLevels; chunks.tones.carrierFreqLevels(1)]);

[hfig, hax] = figureTracesI_PP_pixels( nK, 200, 50 );
hax = flipud(hax); % 1 is bottom, as indendrogram   %check
for Ki = 1:nK   %fix
    axes(hax(Ki)), hold on
    K2use = find(id_colors==Ki);
    
    % amplitude 3: pips
    yvar3 = [fTC(3).dffPeak(:,:,K2use), NaN(1,8)]; %dff values
    yvar3 = yvar3(ia);
    % add tone 4Hz
    yvar2 = fTC(5).dffPeak(:,:,K2use); %dff values
    % put together
    yvar = [yvar2(1), yvar3(2:end)];
    TCfinal(K2use,:) = yvar; %NOTE K2USE!!!!!
    plot(xvar, yvar, 'Color', cmap(K2use,:))
%     plot(xvar(1), yvar(1), 'x', 'Color', cmap(K2use,:), 'MarkerSize',3)
    plot([xvar(1), xvar(end)],[0 0], '-k');   
    
    hax(Ki).XLim = [xvar(1), xvar(end)];
    hax(Ki).YLim = dffLimits(K2use, :).*1;
    hax(Ki).Box = 'off';
    hax(Ki).XAxis.Visible = 'off';
    hax(Ki).YTick = [];
    hax(Ki).XAxis.Scale = 'log';
    hax(Ki).XMinorGrid = 'off';
    set(gca, 'XTick', xvar)
    set(gca, 'XGrid', 'on')
    hax(Ki).TickDir = 'out';
    hax(Ki).XLim = [0,xvar(end)];
    if Ki == 1
        hax(Ki).XTickLabel = xvar;
        hax(Ki).XTickLabelRotation = 90;
        hax(Ki).XAxis.FontSize = 4;
        hax(Ki).XAxis.Visible = 'on';
    else
        hax(Ki).XTickLabel = [];
        hax(Ki).XAxis.Visible = 'off';
    end
    hax(Ki).XMinorTick = 'off';
    drawnow
end
export_fig(fullfile(folder2figure4, 'tuningCurves_01_06sec_4HzAndPips_MaxAmplit_65%.eps'))

%% find the maximum at highest amplitude
[~, peakFreqIdx] = max(TCfinal, [], 2);
[hfig, hax] = figureTracesI_PP_pixels( nK, 100, 50 );
hax = flipud(hax); 

for Ki = 1:nK 
    K2use = find(id_colors==Ki);
    if peakFreqIdx(K2use) == 1 
        peakFreqIdx(K2use) = 2; % the closest
    end
    % amplitude 1: pips
    yvar1 = [fTC(1).dffPeak(:,:,K2use), NaN(1,8)]; %dff values
    yvar1 = yvar1(ia);
    yvar1 = yvar1(peakFreqIdx(K2use));
    
    % amplitude 2: pips
    yvar2 = [fTC(2).dffPeak(:,:,K2use), NaN(1,8)]; %dff values
    yvar2 = yvar2(ia);
    yvar2 = yvar2(peakFreqIdx(K2use));
    
    % amplitude 3: pips
    yvar3 = [fTC(3).dffPeak(:,:,K2use), NaN(1,8)]; %dff values
    yvar3 = yvar3(ia);
    yvar3 = yvar3(peakFreqIdx(K2use));
    
    % plot
    axes(hax(Ki)), hold on
    plot(1:3, [yvar1, yvar2, yvar3], 'Color', cmap(K2use,:))
    
    plot([0, 4],[0 0], '-k');
    hax(Ki).XLim = [0,4];
    hax(Ki).YLim = dffLimits(K2use, :).*1;
    hax(Ki).Box = 'off';
    hax(Ki).XAxis.Visible = 'off';
    hax(Ki).YTick = [];
%     hax(Ki).XAxis.Scale = 'lin';
    hax(Ki).XMinorGrid = 'off';
    set(gca, 'XTick', 1:3)
    set(gca, 'XGrid', 'on')
    hax(Ki).TickDir = 'out';
    hax(Ki).YAxis.Visible = 'off';
    if Ki == 1
        hax(Ki).XTickLabel = {'low', 'mid', 'high'};
        hax(Ki).XTickLabelRotation = 90;
        hax(Ki).XAxis.FontSize = 5;
        hax(Ki).XAxis.Visible = 'on';
    else
        hax(Ki).XTickLabel = [];
        hax(Ki).XAxis.Visible = 'off';
    end
    hax(Ki).XMinorTick = 'off';
    drawnow 
end
export_fig(fullfile(folder2figure4, 'tuningAmplitude_01_06sec_65%.eps'))

