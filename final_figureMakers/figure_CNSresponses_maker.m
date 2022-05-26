load_wedAmmc_downsampled;
load('/Users/galileo/Dropbox (HMS)/Data/chunks_stimfinal.mat', 'chunks')

iRun = 14;
aZ = matfile(datalist{iRun}); %174_run01
T = aZ.T;
fastStimMeta = aZ.fastStimMeta;
clear aZ
for z = 1:NaZ
    aZ{z} = matfile(datalist{z});
end
load(fullfile(Folder2Save, 'R_matrix_Downsampled_wedammc.mat'), 'pxKeep')
load(fullfile(Folder2Save, 'R_matrix_sortedTones_Downsampled_smallWindow.mat'), 'R', 'pixel2keep', 'iKeep')
R = R(pxKeep,:);

clusterFolder = fullfile(Folder2Save, '/rawDataLinkage/ward_ZSCpix');
cd(clusterFolder)
assert(size(R,1) == sum(pxKeep))

folder2figure = '/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/figure5n';
if ~exist(folder2figure, 'dir')
    mkdir(folder2figure)
end

load('klust_maxClust19.mat'), nK = length(unique(Tparent)); 

%% 0 . una tantum change cluster identity (again) and for all. Sorted from top to bottom axis -DONE!!
% load('/Users/galileo/Dropbox (HMS)/Data/cmap19_spectral.mat', 'cmap') % no need to change this anymore.
% % clustIdSorting = fliplr([18 17 19 10 13 9 12 11 15 14 8 16 6 5 4 3 2 1 7]);
% clustIdSorting = [1:12 14 13 15:19];
% 
% TparentNew = nan(size(Tparent));
% for k = 1:nK
%     TparentNew(Tparent==clustIdSorting(k)) = k;
% end
% for k = 1:nK
%     klustNew(k).k = klust(clustIdSorting(k)).k;
% end
% for k = 1:nK
%     dendrOrder(k) = find(clustIdSorting==k);  %correct without flippings
% end
% Tparent = TparentNew;
% klust = klustNew;
% clear TparentNew klustNew
% save('klust_maxClust19.mat', 'Tparent', 'klust', 'cmap', 'dendrOrder')

% warning('need to recalculate fTC if swapped clusters')

%% recalculate FT data with new k labels (faccio prima, evito errori)
clust = remapIndexedPointsToRun_RoisMaps(Tparent, pixel2keep, sizeMap);
ROIs = clust.ROIs; %these are still unsorted, so do not remove the last one!

trialTypes = {chunks.pips.amplitude_1, chunks.pips.amplitude_2, chunks.pips.amplitude_3, ...
              chunks.tones.amplitude_1_carrPhase_1, chunks.tones.amplitude_2_carrPhase_1 };
% %%
% saveNameFTData = fullfile(clusterFolder, 'FTdata.mat');
% ts_dec = T(3).ts_dec;
% saveName = fullfile(Folder2Save, 'R_chunks.mat');
% 
% 
% freqTuning_pipsTones_extendedBaselines_extractTraces(trialTypes, aZ, saveNameFTData); % only run if not already.
% movingwindowsize = 1;   %for moving average
% responseWindow = 1;     %seconds after stim onset. Empty [] for entire default window
% if ~isempty(responseWindow)
%     responseWindowPoints = sum(ts_dec>=0 & ts_dec <=responseWindow);
% else
%     responseWindowPoints = [];
% end
% 
% %this is currently hijacked in script:
% fTC = freqTuning_pipsTones_extendedBaselines(trialTypes, ROIs, movingwindowsize, saveNameFTData, responseWindowPoints); % baseline+ext baseline: 102 points
% % load('fTC_data_01_06secAfterOnset.mat')


%% 1. ONLY check peak order for now.
load(fullfile(clusterFolder,'fTC_data_01_06secAfterOnset.mat'))
[xvar, ia, ic] = unique([chunks.pips.carrierFreqLevels; chunks.tones.carrierFreqLevels(1)]);

for i = 1 : length(fTC)
    for k = 1 : nK
        fTC(i).dffPeak(:,:,k) = ( (max(fTC(i).movAvgResponse(:,:,k)) - fTC(i).meanNoise(:,:,k)) ./ fTC(i).meanNoise(:,:,k) ).* 100;
    end
end
%get maximum
for Ki = 1:nK 
    % amplitude 3: pips
    yvar3 = [fTC(3).dffPeak(:,:,Ki), NaN(1,8)]; %dff values
    yvar3 = yvar3(ia);
    % add tone 4Hz
    yvar2 = fTC(5).dffPeak(:,:,Ki); %dff values
    % put together
    yvar = [yvar2(1), yvar3(2:end)];
    TCfinal(Ki,:) = yvar;
end
[~, peakFreqIdx] = max(TCfinal, [], 2);
disp('should be sorted:')
peakFreqIdx %should be sorted

% ok

% now recomputing everything with new cluster sorting, to avoid any errors
% I DO NOT ACTUALLY NEED TO RECALCULATE R_chunks
% also expand tones range


%% try considering cdfs
% first, interpolate tuning curves at linearly-or-log spaced intervals.
onset = 4;
offset = 600; %cropping before, makes things worse
cutoff = 0.5;

x2 = logspace(log10(onset), log10(offset), 2000);
clear cdfs

figure; hold on
plot([onset,offset], [cutoff cutoff], ':', 'Color', [0.6 0.6 0.6])
for Ki = 1:nK
    yvar2 = interp1(xvar, TCfinal(Ki,:), x2 );
    cdfs(Ki,:) = cumsum(yvar2)./max(cumsum(yvar2));
    plot(x2, cdfs(Ki,:), 'Color', cmap(Ki,:));
%     xvar3 = 1:0.01:19;
%     yvar3 = interp1(1:19, cdfs(Ki,:), xvar3 );
    [~, halff] = find(cdfs(Ki,:)>=cutoff,1);
    halfpoint(Ki) = x2(halff);
    plot([halfpoint(Ki), halfpoint(Ki)], [0,cutoff], 'Color', cmap(Ki,:));
    drawnow;
end
[a,b] = sort(halfpoint);
b
h = gca;
h.XTick = a;
h.XTickLabel = b;
h.TickDir = 'out';

xlim([onset,offset])
ylabel('frequency tuning cdf')
xlabel('frequency, log-spaced between 4 and 600 Hz. Ticks represent cluster # after sorting accordingly - except 1, 2, 19')

% export_fig(fullfile(folder2figure, 'cdf_freqTuningCurved.eps'))



h.XTick =  [4, 8, 16, 25, 50, 100, 150, 225, 350, 600];
h.XTickLabel =  [4, 8, 16, 25, 50, 100, 150, 225, 350, 600];
h.XGrid = 'on';
xlabel('');

h.XAxis.Scale = 'log';
h.XLim = [4,600];

h.XAxis.MinorTick = 'off';
h.XMinorGrid = 'off';
% export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/cdf_freqTuningCurved_CNS_LOG.eps')

%% chunk's R
cropbaseline = - 1.7; %first point give errors
saveName = fullfile(Folder2Save, 'R_chunks.mat');
ylims = [-20, 180];

%% chunk 1: frequency sweeps
CN = 'chirps_down_up';
% chunks.(CN)
% stimDur = 12;
sortedStimNs = {62, 61}
sortedClassNs = zeros(1, length(sortedStimNs));
for i = 1: length(sortedStimNs)
    sortedClassNs(i) = aZ{1}.stimuli2Class(sortedStimNs{i}(1), 1);
end
% [R_chirps, ~] = makeSaveR_mf_sortedStimuli( aZ, sets, pixel2keep, sortedStimNs, sortedClassNs, cropbaseline );
% save(saveName, 'R_chirps') %will be updated
load(saveName, 'R_chirps')

% CHIRPS - plot
[stimTraces, ts_Stim_all, stN, clN, zeros_st_idxs] = makeSortedStimulusSet_protot(T, sortedStimNs, sortedClassNs, cropbaseline);
stimTraces = stimTraces.*3; %um

sortedStimNs = cat(2, sortedStimNs{:});
i = 1;
effectiveFirst  = T(sortedClassNs(i)).ts_dec( find(T(sortedClassNs(i)).ts_dec >= cropbaseline, 1, 'first') );
effectiveLast   = T(sortedClassNs(i)).ts_dec( find(T(sortedClassNs(i)).ts_dec <= T(sortedClassNs(i)).stimDur + 10, 1, 'last') );
iKeepResponse = T(sortedClassNs(i)).ts_dec >= effectiveFirst   & T(sortedClassNs(i)).ts_dec <= effectiveLast;
ts_dec = T(sortedClassNs(i)).ts_dec(iKeepResponse) ;
deltaEndTss = T(sortedClassNs(i)).ts(end) - T(sortedClassNs(i)).ts_dec(end);
deltaStartTss = -(cropbaseline - ts_dec(1));

ts_dec_all = ts_dec - ts_dec(1) + deltaStartTss; %now it is correct!!
zero_tr_idxs = find(ts_dec>=0,1);
stN_tr = sortedStimNs(i) *ones(length(ts_dec), 1);

for i = 2 : length(sortedStimNs) %assuming they ally belong to the same or analogous class
    zero_tr_idxs = cat(1, zero_tr_idxs, find(ts_dec>=0,1) + length(ts_dec_all));
    ts_dec_all = cat(1, ts_dec_all, ts_dec - ts_dec(1) + ts_dec_all(end) + deltaStartTss + deltaEndTss );
    stN_tr = cat(1, stN_tr, sortedStimNs(i)*ones(length(ts_dec), 1) );
end

% start plotting
[hfig, hax] = figureTracesI_PP_pixels( nK, 200, 50 );
for Ki = 1:nK
    axes(hax(Ki))
    hold off
    allR = R_chirps(Tparent==Ki,:);
    avg = mean(allR,1);
    plot([ts_dec_all(1); ts_dec_all; ts_dec_all(end); ts_dec_all(1)], [0, avg, 0, 0], 'Color',  cmap(Ki,:)), hold on, axis tight   
    hax(Ki).XLim =  [ts_Stim_all(1), ts_Stim_all(end)];
    hax(Ki).XAxis.Visible = 'off';
    hax(Ki).Box = 'off';
    drawnow
end
for Ki = 1:nK
    hax(Ki).YLim = ylims;
    hax(Ki).YTick = [];
    hax(Ki).YAxis.Visible = 'off';
end
YaxReference = 100; %floor(diff(ylims)/100)*100;
plot([0,0], [ylims(1), ylims(1) + YaxReference], '-k')
hax(Ki).YAxis.Label.String = YaxReference;
hax(Ki).YAxis.Label.Visible = 'on';
hax(Ki).YAxis.Label.FontSize = 9;
XaxReference = min([floor(hax(Ki).XLim(end)/2)*2, 10]);
plot([0,XaxReference], [ylims(1), ylims(1)], '-k')
hax(Ki).XTick = [];
hax(Ki).XAxis.Label.String = sprintf('%d / %.2f',XaxReference, hax(Ki).XLim(end));
hax(Ki).XAxis.Label.Visible = 'on';
hax(Ki).XAxis.Label.FontSize = 9;
pause

% export and save!
export_fig(fullfile(folder2figure, 'resp_chirps.eps'))

%% chunk 2: adaptation to long tones
CN = 'tones';
% chunks.(CN)
sortedStimNs = {[85 89 93 75 77 79 81 83], [86 90 94 76 78 80 82 84]};
% sortedStimNs = {76};
sortedClassNs = zeros(1, length(sortedStimNs));
for i = 1: length(sortedStimNs)
    sortedClassNs(i) = aZ{1}.stimuli2Class(sortedStimNs{i}(1), 1);
end
% % tic
% [R_tones, ~] = makeSaveR_mf_sortedStimuli( aZ, sets, pixel2keep, sortedStimNs, sortedClassNs, cropbaseline );
% % toc
% save(saveName, 'R_tones', '-append')
load(saveName, 'R_tones');
[stimTraces, ts_Stim_all, stN, clN, zeros_st_idxs] = makeSortedStimulusSet_protot(T, sortedStimNs, sortedClassNs, cropbaseline);
stimTraces = stimTraces.*3; %um
sortedStimNs = cat(2, sortedStimNs{:});

%% TONE(s) - plot
% timestamps!
i = 1;
effectiveFirst  = T(sortedClassNs(i)).ts_dec( find(T(sortedClassNs(i)).ts_dec >= cropbaseline, 1, 'first') );
effectiveLast   = T(sortedClassNs(i)).ts_dec( find(T(sortedClassNs(i)).ts_dec <= T(sortedClassNs(i)).stimDur + 10, 1, 'last') );
iKeepResponse = T(sortedClassNs(i)).ts_dec >= effectiveFirst   & T(sortedClassNs(i)).ts_dec <= effectiveLast;
ts_dec = T(sortedClassNs(i)).ts_dec(iKeepResponse) ;
deltaEndTss = T(sortedClassNs(i)).ts(end) - T(sortedClassNs(i)).ts_dec(end);
deltaStartTss = -(cropbaseline - ts_dec(1));

ts_dec_all = ts_dec - ts_dec(1) + deltaStartTss; %now it is correct!!
zero_tr_idxs = find(ts_dec>=0,1);
stN_tr = sortedStimNs(i) *ones(length(ts_dec), 1);

for i = 2 : length(sortedStimNs) %assuming they all belong to the same or analogous class
    zero_tr_idxs = cat(1, zero_tr_idxs, find(ts_dec>=0,1) + length(ts_dec_all));
    ts_dec_all = cat(1, ts_dec_all, ts_dec - ts_dec(1) + ts_dec_all(end) + deltaStartTss + deltaEndTss );
    stN_tr = cat(1, stN_tr, sortedStimNs(i)*ones(length(ts_dec), 1) );
end

%so far all the stimuli from the previos section are included (so both R
%and stimTraces are in the same framework. However, only selectStimuli_stN
%will be plotted down here.
% start plotting
sel_stN = [86 90 94 80];
i_St = ismember(stN, sel_stN);
i_Tr = ismember(stN_tr, sel_stN);

ts_Stim_all_sel = ts_Stim_all(1:sum(i_St));
% onset_idx_sel = find(ismember( ts_Stim_all_sel, ts_Stim_all(zeros_st_idxs)));
ts_Stim_all_sel = ts_Stim_all_sel - ts_Stim_all_sel(1);
ts_dec_all_sel = ts_dec_all(1:sum(i_Tr));
ts_dec_all_sel = ts_dec_all_sel-ts_dec_all_sel(1);


% % plot stimulus
% YLIM = 15;
% [hfig, hax] = figureTracesI_PP_pixels( 1, 200, 50 );
% axes(hax(end)), hold on
% plot(ts_Stim_all_sel, stimTraces(i_St), '-k')
% hax(end).YLim = [-YLIM YLIM]; %um max range, fixed across stimuli
% hax(end).XLim =  [ts_Stim_all_sel(1), ts_Stim_all_sel(end)];
% sts = unique(sel_stN, 'stable');
% st_onoff_ts = [];
% xticks = [];
% hax(end).XTick = [];
% for n = 1:length(sts)
%     xticks = cat(2, xticks, ts_Stim_all_sel( find(stN==sortedStimNs(n), 1) ) ); %trial divisor
%     st_onoff_ts = cat(2, st_onoff_ts, ts_Stim_all_sel(zeros_st_idxs(n))); 
%     st_onoff_ts = cat(2, st_onoff_ts, ts_Stim_all_sel(zeros_st_idxs(n))+ 2 );
%     st_on = ts_Stim_all_sel(zeros_st_idxs(n)); 
%     st_off =ts_Stim_all_sel(zeros_st_idxs(n))+ 2;
%     plot([st_on, st_off, st_off, st_on, st_on], [-YLIM -YLIM YLIM YLIM -YLIM], 'Color', [0.7 0.7 0.7])
% end
% hax(end).XTick = unique(xticks);
% hax(end).XMinorTick='off';
% hax(end).XGrid = 'on';
% hax(end).XAxis.Visible = 'off';
% hax(end).YAxis.Visible = 'off';
% 
% export_fig(fullfile(folder2figure, sprintf('stim_tones_%s_4stSubset_HighAmpl.eps',num2str(sel_stN))))



[hfig, hax] = figureTracesI_PP_pixels( nK, 200, 50 );
for Ki = 1:nK
    axes(hax(Ki))
    hold off

    allR = R_tones(Tparent==Ki,:);
    avg = mean(allR,1);
    plot([ts_dec_all_sel(1); ts_dec_all_sel; ts_dec_all_sel(end); ts_dec_all_sel(1)], [0, avg(i_Tr), 0, 0], 'Color',  cmap(Ki,:)), hold on, axis tight

    hax(Ki).XLim =  [ts_Stim_all_sel(1), ts_Stim_all_sel(end)];
    hax(Ki).XAxis.Visible = 'off';
    hax(Ki).Box = 'off';
    drawnow
end

for Ki = 1:nK
    hax(Ki).YLim = ylims;
    hax(Ki).YTick = [];
    hax(Ki).YAxis.Visible = 'off';
end
YaxReference = 100; %floor(diff(ylims)/100)*100;
plot([0,0], [ylims(1), ylims(1) + YaxReference], '-k')
hax(Ki).YAxis.Label.String = YaxReference;
hax(Ki).YAxis.Label.Visible = 'on';
hax(Ki).YAxis.Label.FontSize = 9;
XaxReference = min([floor(hax(Ki).XLim(end)/2)*2, 10]);
plot([0,XaxReference], [ylims(1), ylims(1)], '-k')
hax(Ki).XTick = [];
hax(Ki).XAxis.Label.String = sprintf('%d / %.2f',XaxReference, hax(Ki).XLim(end));
hax(Ki).XAxis.Label.Visible = 'on';
hax(Ki).XAxis.Label.FontSize = 9;
pause

export_fig(fullfile(folder2figure, sprintf('resp_tones_%s_4stSubset_HighAmpl.eps',num2str(sel_stN))))


%% chunk 3: songs
CN = 'courtships';
% chunks.(CN)
% stimDur = 12;
sortedStimNs = {74, 73}
sortedClassNs = zeros(1, length(sortedStimNs));
for i = 1: length(sortedStimNs)
    sortedClassNs(i) = aZ{1}.stimuli2Class(sortedStimNs{i}(1), 1);
end
% [R_songs, ~] = makeSaveR_mf_sortedStimuli( aZ, sets, pixel2keep, sortedStimNs, sortedClassNs, cropbaseline );
% save(saveName, 'R_songs', '-append') %will be updated
load(saveName, 'R_songs')
R_songs_remain = R_songs;

[stimTraces, ts_Stim_all, stN, clN, zeros_st_idxs] = makeSortedStimulusSet_protot(T, sortedStimNs, sortedClassNs, cropbaseline);
stimTraces = stimTraces.*3; %um
sortedStimNs = cat(2, sortedStimNs{:});

% SONGS - plot
for i = 1:2
    sel_stN = sortedStimNs(i);
    stimDur = T(sortedClassNs(i)).stimDur - 2.6;
    effectiveFirst  = T(sortedClassNs(i)).ts_dec( find(T(sortedClassNs(i)).ts_dec >= cropbaseline, 1, 'first') );
    effectiveLast   = T(sortedClassNs(i)).ts_dec( find(T(sortedClassNs(i)).ts_dec <= T(sortedClassNs(i)).stimDur + 10, 1, 'last') );
    iKeepResponse = T(sortedClassNs(i)).ts_dec >= effectiveFirst   & T(sortedClassNs(i)).ts_dec <= effectiveLast;
    ts_dec = T(sortedClassNs(i)).ts_dec(iKeepResponse) ;
    deltaEndTss = T(sortedClassNs(i)).ts(end) - T(sortedClassNs(i)).ts_dec(end);
    deltaStartTss = -(cropbaseline - ts_dec(1));
    
    ts_dec_all_sel = ts_dec - ts_dec(1) + deltaStartTss; %now it is correct!!
    zero_tr_idxs = find(ts_dec>=0,1);
    stN_tr = sortedStimNs(i) *ones(length(ts_dec), 1);
    i_St = ismember(stN, sel_stN);
    ts_Stim_all_sel = ts_Stim_all(i_St);
    % onset_idx_sel = find(ismember( ts_Stim_all_sel, ts_Stim_all(zeros_st_idxs)));
    ts_Stim_all_sel = ts_Stim_all_sel - ts_Stim_all_sel(1);
    R_songs = R_songs_remain(:,1:length(ts_dec));
    R_songs_remain = R_songs_remain(:,length(ts_dec)+1:end);
    
    [hfig, hax] = figureTracesI_PP_pixels( nK, 200, 50 );
    for Ki = 1:nK
        axes(hax(Ki))
        hold off
        allR = R_songs(Tparent==Ki,:);
        avg = mean(allR,1);
        plot([ts_dec_all_sel(1); ts_dec_all_sel; ts_dec_all_sel(end); ts_dec_all_sel(1)], [0, avg, 0, 0], 'Color',  cmap(Ki,:)), hold on, axis tight
        % plot([ts_dec_all_sel(1), ts_dec_all_sel(end)], [0,0], 'Color', colors(Ki,:) )
        
        hax(Ki).XLim =  [ts_Stim_all_sel(1), ts_Stim_all_sel(end)];
        hax(Ki).XAxis.Visible = 'off';
        hax(Ki).Box = 'off';
        drawnow
    end
  
    for Ki = 1:nK
        hax(Ki).YLim = ylims;
        hax(Ki).YTick = [];
        hax(Ki).YAxis.Visible = 'off';
    end
    YaxReference = 100; %floor(diff(ylims)/100)*100;
    plot([0,0], [ylims(1), ylims(1) + YaxReference], '-k')
    hax(Ki).YAxis.Label.String = YaxReference;
    hax(Ki).YAxis.Label.Visible = 'on';
    hax(Ki).YAxis.Label.FontSize = 9;
    XaxReference = min([floor(hax(Ki).XLim(end)/2)*2, 10]);
    plot([0,XaxReference], [ylims(1), ylims(1)], '-k')
    hax(Ki).XTick = [];
    hax(Ki).XAxis.Label.String = sprintf('%d / %.2f',XaxReference, hax(Ki).XLim(end));
    hax(Ki).XAxis.Label.Visible = 'on';
    hax(Ki).XAxis.Label.FontSize = 9;
    
    
    pause
    % export and save!
    export_fig(fullfile(folder2figure, sprintf('resp_songs_%s.eps',num2str(sel_stN))))
end

%% chunk 4: steps
CN = 'steps';
chunks.(CN)
% stimDur = 12;
sortedStimNs = {55:60}
sortedClassNs = zeros(1, length(sortedStimNs));
for i = 1: length(sortedStimNs)
    sortedClassNs(i) = aZ{1}.stimuli2Class(sortedStimNs{i}(1), 1);
end
load(saveName, 'R_steps')

% STEPS - plot
% timestamps!
[stimTraces, ts_Stim_all, stN, clN, zeros_st_idxs] = makeSortedStimulusSet_protot(T, sortedStimNs, sortedClassNs, cropbaseline);
stimTraces = stimTraces.*3; %um

sortedStimNs = cat(2, sortedStimNs{:});
i = 1;
effectiveFirst  = T(sortedClassNs(i)).ts_dec( find(T(sortedClassNs(i)).ts_dec >= cropbaseline, 1, 'first') );
effectiveLast   = T(sortedClassNs(i)).ts_dec( find(T(sortedClassNs(i)).ts_dec <= T(sortedClassNs(i)).stimDur + 10, 1, 'last') );
iKeepResponse = T(sortedClassNs(i)).ts_dec >= effectiveFirst   & T(sortedClassNs(i)).ts_dec <= effectiveLast;
ts_dec = T(sortedClassNs(i)).ts_dec(iKeepResponse) ;
deltaEndTss = T(sortedClassNs(i)).ts(end) - T(sortedClassNs(i)).ts_dec(end);
deltaStartTss = -(cropbaseline - ts_dec(1));

ts_dec_all = ts_dec - ts_dec(1) + deltaStartTss; %now it is correct!!
zero_tr_idxs = find(ts_dec>=0,1);
stN_tr = sortedStimNs(i) *ones(length(ts_dec), 1);

for i = 2 : length(sortedStimNs) %assuming they ally belong to the same or analogous class
    zero_tr_idxs = cat(1, zero_tr_idxs, find(ts_dec>=0,1) + length(ts_dec_all));
    ts_dec_all = cat(1, ts_dec_all, ts_dec - ts_dec(1) + ts_dec_all(end) + deltaStartTss + deltaEndTss );
    stN_tr = cat(1, stN_tr, sortedStimNs(i)*ones(length(ts_dec), 1) );
end


% start plotting
[hfig, hax] = figureTracesI_PP_pixels( nK, 200, 50 );
for Ki = 1:nK
    axes(hax(Ki))
    hold off
    allR = R_steps(Tparent==Ki,:);
    avg = mean(allR,1);
%     cmap(Ki,:)
    plot([ts_dec_all(1); ts_dec_all; ts_dec_all(end); ts_dec_all(1)], [0, avg, 0, 0], 'Color',  cmap(Ki,:)), hold on, axis tight
    
    hax(Ki).XLim =  [ts_Stim_all(1), ts_Stim_all(end)];
    hax(Ki).XAxis.Visible = 'off';
    hax(Ki).Box = 'off';
    drawnow
end

for Ki = 1:nK
    hax(Ki).YLim = ylims;
    hax(Ki).YTick = [];
    hax(Ki).YAxis.Visible = 'off';
end
YaxReference = 100; %floor(diff(ylims)/100)*100;

plot([0,0], [ylims(1), ylims(1) + YaxReference], '-k')
hax(Ki).YAxis.Label.String = YaxReference;
hax(Ki).YAxis.Label.Visible = 'on';
hax(Ki).YAxis.Label.FontSize = 9;
XaxReference = min([floor(hax(Ki).XLim(end)/2)*2, 10]);
plot([0,XaxReference], [ylims(1), ylims(1)], '-k')
hax(Ki).XTick = [];
hax(Ki).XAxis.Label.String = sprintf('%d / %.2f',XaxReference, hax(Ki).XLim(end));
hax(Ki).XAxis.Label.Visible = 'on';
hax(Ki).XAxis.Label.FontSize = 9;

pause
export_fig(fullfile(folder2figure, 'resp_steps.eps'))



%% New tuning curves

% find max dff values for each cluster
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
%   % 1 is bottom, as indendrogram
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

dffLimits = [min(negMaxAllStimuli,[], 2), max(maxAllStimuli, [],2)]; %this takes Ki

% %% rescale by baseline noise
% % recalculate a single baseline noise level for each cluster
% for k = 1 : size(ROIs,2)
%     basTemp = [];
%     for i = 1:5
%         basTemp = cat(1, basTemp, reshape(fTC(i).baselines(:,:,k), [],1,1));
%     end
%     fullBaselPipsTones(:, k) = basTemp;
% end
% noiseClusters = std(fullBaselPipsTones); %ok.
% scaledNoiseClusters = ( (noiseClusters-min(noiseClusters)) ./ max((noiseClusters-min(noiseClusters))) ) / 2 ; %scaled between 0 (min noise) and 0.5 (max noise). [To be subtracted]
% % scaledNoiseClusters =
% %     0.0434    0.0339    0.0133    0.5000    0.1400    0.1466    0.0201    0.0395    0.2126    0.0610         0

%% single curve with 4Hz from tones and rest from pips, max amplitude only
[xvar, ia, ic] = unique([chunks.pips.carrierFreqLevels; chunks.tones.carrierFreqLevels(1)]);

[hfig, hax] = figureTracesI_PP_pixels( nK, 200, 50 );
for Ki = 1:nK   %fix
    axes(hax(Ki)), hold on   
    % amplitude 3: pips
    yvar3 = [fTC(3).dffPeak(:,:,Ki), NaN(1,8)]; %dff values
    yvar3 = yvar3(ia);
    % add tone 4Hz
    yvar2 = fTC(5).dffPeak(:,:,Ki); %dff values
    % put together
    yvar = [yvar2(1), yvar3(2:end)];
    TCfinal(Ki,:) = yvar; %NOTE Ki!!!!!
    plot(xvar, yvar, 'Color', cmap(Ki,:))
%     plot(xvar(1), yvar(1), 'x', 'Color', cmap(Ki,:), 'MarkerSize',3)
    plot([xvar(1), xvar(end)],[0 0], '-k');   
    
    hax(Ki).XLim = [xvar(1), xvar(end)];
    hax(Ki).YLim = dffLimits(Ki, :).*1;
    hax(Ki).Box = 'off';
    hax(Ki).XAxis.Visible = 'off';
    hax(Ki).YTick = [];
    hax(Ki).XAxis.Scale = 'log';
    hax(Ki).XMinorGrid = 'off';
    set(gca, 'XTick', xvar)
    hax(Ki).XGrid = 'off'; %was on
    hax(Ki).TickDir = 'out';
    hax(Ki).XLim = [0,xvar(end)];
    if Ki == nK
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
pause
export_fig(fullfile(folder2figure, 'tuningCurves_01_06sec_4HzAndPips_MaxAmplit_100%.eps'))

%% find the maximum at highest amplitude
[~, peakFreqIdx] = max(TCfinal, [], 2);


[hfig, hax] = figureTracesI_PP_pixels( nK, 80, 50 );
for Ki = 1:nK 
    if peakFreqIdx(Ki) == 1 
        peakFreqIdx(Ki) = 2; % the closest
    end
    % amplitude 1: pips
    yvar1 = [fTC(1).dffPeak(:,:,Ki), NaN(1,8)]; %dff values
    yvar1 = yvar1(ia);
    yvar1 = yvar1(peakFreqIdx(Ki));
    
    % amplitude 2: pips
    yvar2 = [fTC(2).dffPeak(:,:,Ki), NaN(1,8)]; %dff values
    yvar2 = yvar2(ia);
    yvar2 = yvar2(peakFreqIdx(Ki));
    
    % amplitude 3: pips
    yvar3 = [fTC(3).dffPeak(:,:,Ki), NaN(1,8)]; %dff values
    yvar3 = yvar3(ia);
    yvar3 = yvar3(peakFreqIdx(Ki));
    
    
    
    amplCurves(Ki,1) = yvar1;
    amplCurves(Ki,2) = yvar2;
    amplCurves(Ki,3) = yvar3;
    
    
    % plot
    axes(hax(Ki)), hold on
    plot(1:3, [yvar1, yvar2, yvar3], 'Color', cmap(Ki,:))
    
    plot([0.8, 3.2],[0 0], '-k');
    hax(Ki).XLim = [0.8, 3.2];
    hax(Ki).YLim = dffLimits(Ki, :).*1;
    hax(Ki).Box = 'off';
    hax(Ki).XAxis.Visible = 'off';
    hax(Ki).YTick = [];
%     hax(Ki).XAxis.Scale = 'lin';
    hax(Ki).XMinorGrid = 'off';
    set(gca, 'XTick', 1:3)
    hax(Ki).XGrid = 'off';
    hax(Ki).TickDir = 'out';
    hax(Ki).YAxis.Visible = 'on';
    if Ki == nK
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


save('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/amplitudeBestCDFStim_CNS.mat', 'amplCurves')


pause
export_fig(fullfile(folder2figure, 'tuningAmplitude_01_06sec_100%.eps'))





%% now add cluster heatmaps - ventral vs dorsal planes
edit WEDvsAMMC_clusterdistribution.m


