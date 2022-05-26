%% clean stimuli (command)
load_wedAmmc_downsampled;
load('/Users/galileo/Dropbox (HMS)/Data/chunks_stimfinal.mat', 'chunks')
% load('/Users/galileo/Dropbox (HMS)/Data_Raw_Metadata_BackedUp/BackedUp_Metadata/fly171_run01_metadata_2016-12-20_163728.mat', 'stimuli')


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
treeFileName = 'output_ward_Rzsp_linkage_wedammc.mat';
cd(clusterFolder)
load('klust_maxClust19.mat'), nK = length(unique(Tparent)); %ne' i colori ne' i clusters sono ordinati, ma tra di loro sono concordanti. Usa id_colors per tonotopy order
assert(size(R,1) == sum(pxKeep))
% clusterFolder = fullfile(clusterFolder, sprintf('%dclusters', nK));
% cd(clusterFolder)
% load('alignedMatrix_Images_maxClust15.mat')
% colors = Col(2:end,:);
% assert(size(R,1) == sum(pxKeep))

load('/Users/galileo/Dropbox (HMS)/Data/cmap19_spectral.mat', 'cmap')
id_colors = [2 3 4 5 6 7 1 19 13 15 11 12 14 9 10 8 16 18 17];
cmap = cmap(id_colors, :);

save('klust_maxClust19.mat', 'id_colors', 'cmap', '-append');
% folder2stimuli = '/Users/galileo/Dropbox (HMS)/figures/figure2';
% folder2figure1 = '/Users/galileo/Dropbox (HMS)/figures/figure1';
folder2figure4 = '/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/figure4';
% folder2figure3 = '/Users/galileo/Dropbox (HMS)/figures/figure3';

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
%% CHIRPS - plot
% timestamps!
[stimTraces, ts_Stim_all, stN, clN, zeros_st_idxs] = makeSortedStimulusSet_protot(T, sortedStimNs, sortedClassNs, cropbaseline);
stimTraces = stimTraces.*3; %um
YLIM = 15;

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
[hfig, hax] = figureTracesI_PP_pixels( nK, 200, 30 );
hax = flipud(hax); % 1 is bottom, as indendrogram
% % plot stimulus
% axes(hax(end)), hold on
% % plot(ts_Stim_all, stimTraces, '-k')
% % hax(end).YLim = [-YLIM YLIM]; %um max range, fixed across stimuli
% % hax(end).XLim =  [ts_Stim_all(1), ts_Stim_all(end)];
% sts = unique(stN, 'stable');
% st_onoff_ts = [];
% xticks = [];
% % hax(end).XTick = [];
% for n = 1:length(sts)
%     xticks = cat(2, xticks, ts_Stim_all( find(stN==sortedStimNs(n), 1) ) ); %trial divisor
%     st_onoff_ts = cat(2, st_onoff_ts, ts_Stim_all(zeros_st_idxs(n))); 
%     st_onoff_ts = cat(2, st_onoff_ts, ts_Stim_all(zeros_st_idxs(n))+ T(sortedClassNs(n)).stimDur );
%     st_on = ts_Stim_all(zeros_st_idxs(n)); 
%     st_off =ts_Stim_all(zeros_st_idxs(n))+ T(sortedClassNs(n)).stimDur;
%     plot([st_on, st_off, st_off, st_on, st_on], [-YLIM -YLIM YLIM YLIM -YLIM], 'Color', [0.7 0.7 0.7])
% end
% % hax(end).XTick = unique(xticks);
% % hax(end).XMinorTick='off';
% % hax(end).XGrid = 'on';
% % hax(end).XAxis.Visible = 'off';
% % hax(end).YAxis.Visible = 'off';


for Ki = 1:nK
    axes(hax(Ki))
    hold off
%     K2use = find(id_colors==Ki);
K2use = newAxisOrder(Ki);
    allR = R_chirps(Tparent==K2use,:);
    avg = mean(allR,1);
    plot([ts_dec_all(1); ts_dec_all; ts_dec_all(end); ts_dec_all(1)], [0, avg, 0, 0], 'Color',  cmap(K2use,:)), hold on, axis tight
    
    hax(Ki).XLim =  [ts_Stim_all(1), ts_Stim_all(end)];
    hax(Ki).XAxis.Visible = 'off';
    hax(Ki).Box = 'off';
    
%     ylims(Ki,:) = hax(Ki).YLim;
    drawnow
end

% ylims = [floor(min(ylims(:,1))/5)*5, ceil(max(ylims(:,2))/5)*5];
% ylims = [-50, 400];
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
export_fig(fullfile(folder2figure4, 'resp_chirps.eps'))
% timestamps_chirps = [];
% timestamps_chirps.ts_dec_all = ts_dec_all;
% timestamps_chirps.stimTraces = stimTraces;
% timestamps_chirps.ts_Stim_all = ts_Stim_all;
% timestamps_chirps.stN = stN;
% timestamps_chirps.zeros_st_idxs = zeros_st_idxs;
% timestamps_chirps.stN_tr = stN_tr;
% save( saveName, 'R_chirps', 'timestamps_chirps', '-append')


% figure
% spectrogram(stimTraces,10000,500,10000,4e4,'yaxis');    % used to be 20000 - 500 - 20000. Trying to distort less the temporal dimension.
% ax = gca;
% ax.TickDir = 'out';
% ax.YLim = [0,0.6];
% ax.YTick = [0, 0.3, 0.6];
% ax.YTickLabel = 0:300:600;
% ax.YLabel.String = 'Frequency (Hz)';
% ax.XTick = [];
% ax.XLabel.String = '';
% ax.Box = 'off';
% np = 128;
% expn = 1.75;
% cmap = brewermap(np, 'Oranges');
% x = 1:np;
% y = x.^expn;
% y = ceil(y./max(y)*np);
% cmap2 = interp1(x, cmap, y);
% colormap(cmap2)
% ax.XLim =  [ts_Stim_all(1), ts_Stim_all(end)];
% export_fig(fullfile(folder2figure4, 'chirps_spectrogram.eps'))
% colorbar off
% ax.Box = 'on';
% ax.XTick = [];
% ax.YTick = [];
% ax.XLabel.Visible = 'off';
% ax.YLabel.Visible = 'off';
% ax.XAxis.Color = [1,1,1];
% ax.YAxis.Color = [1,1,1];
% export_fig(fullfile(folder2figure4, 'chirps_spectrogram.jpg'), '-m4')


%% chunk 2: adaptation to long tones
CN = 'tones';
% chunks.(CN)
sortedStimNs = {[85 89 93 75 77 79 81 83], [86 90 94 76 78 80 82 84]};
% sortedStimNs = {76};
sortedClassNs = zeros(1, length(sortedStimNs));
for i = 1: length(sortedStimNs)
    sortedClassNs(i) = aZ{1}.stimuli2Class(sortedStimNs{i}(1), 1);
end
% tic
[R_tones, ~] = makeSaveR_mf_sortedStimuli( aZ, sets, pixel2keep, sortedStimNs, sortedClassNs, cropbaseline );
% toc
save(saveName, 'R_tones', '-append')
% load(saveName, 'R_tones');
[stimTraces, ts_Stim_all, stN, clN, zeros_st_idxs] = makeSortedStimulusSet_protot(T, sortedStimNs, sortedClassNs, cropbaseline);
stimTraces = stimTraces.*3; %um
YLIM = 15;
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
sel_stN = 76;
i_St = ismember(stN, sel_stN);
i_Tr = ismember(stN_tr, sel_stN);

ts_Stim_all_sel = ts_Stim_all(i_St);
onset_idx_sel = find(ismember( ts_Stim_all_sel, ts_Stim_all(zeros_st_idxs)));
ts_Stim_all_sel = ts_Stim_all_sel - ts_Stim_all_sel(1);
ts_dec_all_sel = ts_dec_all(i_Tr);
ts_dec_all_sel = ts_dec_all_sel-ts_dec_all_sel(1);

[hfig, hax] = figureTracesI_PP_pixels( nK, 200, 30 );
hax = flipud(hax); % 1 is bottom, as indendrogram
% plot stimulus
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
%     st_onoff_ts = cat(2, st_onoff_ts, ts_Stim_all_sel(zeros_st_idxs(n))+ T(sortedClassNs(n)).stimDur );
%     st_on = ts_Stim_all_sel(zeros_st_idxs(n)); 
%     st_off =ts_Stim_all_sel(zeros_st_idxs(n))+ T(sortedClassNs(n)).stimDur;
%     plot([st_on, st_off, st_off, st_on, st_on], [-YLIM -YLIM YLIM YLIM -YLIM], 'Color', [0.7 0.7 0.7])
% end
% hax(end).XTick = unique(xticks);
% hax(end).XMinorTick='off';
% hax(end).XGrid = 'on';
% hax(end).XAxis.Visible = 'off';
% hax(end).YAxis.Visible = 'off';
% 

for Ki = 1:nK
    axes(hax(Ki))
    hold off
    K2use = find(id_colors==Ki);
    allR = R_tones(Tparent==K2use,:);
    avg = mean(allR,1);
    plot([ts_dec_all_sel(1); ts_dec_all_sel; ts_dec_all_sel(end); ts_dec_all_sel(1)], [0, avg(i_Tr), 0, 0], 'Color',  cmap(K2use,:)), hold on, axis tight
    % plot([ts_dec_all_sel(1), ts_dec_all_sel(end)], [0,0], 'Color', colors(Ki,:) )
    
    hax(Ki).XLim =  [ts_Stim_all_sel(1), ts_Stim_all_sel(end)];
    hax(Ki).XAxis.Visible = 'off';
    hax(Ki).Box = 'off';
    
%     ylims(Ki,:) = hax(Ki).YLim;
    drawnow
end

% ylims = [floor(min(ylims(:,1))/5)*5, ceil(max(ylims(:,2))/5)*5];

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
export_fig(fullfile(folder2figure4, sprintf('resp_tones_%s_st76_125HzHighAmpl.eps',num2str(sel_stN))))
% timestamps_tones = [];
% timestamps_tones.ts_dec_all = ts_dec_all;
% timestamps_tones.stimTraces = stimTraces;
% timestamps_tones.ts_Stim_all_sel = ts_Stim_all_sel;
% timestamps_tones.stN = stN;
% timestamps_tones.zeros_st_idxs = zeros_st_idxs;
% timestamps_tones.stN_tr = stN_tr;
% save( saveName, 'R_tones', 'timestamps_tones', '-append')



% 
% figure
% spectrogram(stimTraces(i_St),10000,500,10000,4e4,'yaxis');    % used to be 20000 - 500 - 20000. Trying to distort less the temporal dimension.
% ax = gca;
% ax.TickDir = 'out';
% ax.YLim = [0,0.6];
% ax.YTick = [0, 0.3, 0.6];
% ax.YTickLabel = 0:300:600;
% ax.YLabel.String = 'Frequency (Hz)';
% % ax.XTick = [];
% % ax.XLabel.String = '';
% ax.Box = 'off';
% np = 128;
% expn = 1.75; 
% cmap = brewermap(np, 'Oranges');
% x = 1:np;
% y = x.^expn;
% y = ceil(y./max(y)*np);
% cmap2 = interp1(x, cmap, y);
% colormap(cmap2)
% ax.XLim =  [ts_Stim_all_sel(1), ts_Stim_all_sel(end)];
% export_fig(fullfile(folder2figure4, 'tone_spectrogram.eps'))
% colorbar off
% ax.Box = 'on';
% ax.XTick = [];
% ax.YTick = [];
% ax.XLabel.Visible = 'off';
% ax.YLabel.Visible = 'off';
% ax.XAxis.Color = [1,1,1];
% ax.YAxis.Color = [1,1,1];
% export_fig(fullfile(folder2figure4, 'tone_spectrogram.jpg'), '-m4')




%% chunk 3: songs
CN = 'courtships';
% chunks.(CN)
% stimDur = 12;
sortedStimNs = {74, 73}
sortedClassNs = zeros(1, length(sortedStimNs));
for i = 1: length(sortedStimNs)
    sortedClassNs(i) = aZ{1}.stimuli2Class(sortedStimNs{i}(1), 1);
end
[R_songs, ~] = makeSaveR_mf_sortedStimuli( aZ, sets, pixel2keep, sortedStimNs, sortedClassNs, cropbaseline );
save(saveName, 'R_songs', '-append') %will be updated
% load(saveName, 'R_songs')
R_songs_remain = R_songs;

[stimTraces, ts_Stim_all, stN, clN, zeros_st_idxs] = makeSortedStimulusSet_protot(T, sortedStimNs, sortedClassNs, cropbaseline);
stimTraces = stimTraces.*3; %um
YLIM = 15;
sortedStimNs = cat(2, sortedStimNs{:});

%% SONGS - plot
% timestamps!

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
    
    [hfig, hax] = figureTracesI_PP_pixels( nK, 200, 30 );
    hax = flipud(hax); % 1 is bottom, as indendrogram
    % plot stimulus
%     axes(hax(end)), hold on
%     plot(ts_Stim_all_sel, stimTraces(i_St), '-k')
%     hax(end).YLim = [-YLIM YLIM]; %um max range, fixed across stimuli
%     hax(end).XLim =  [ts_Stim_all_sel(1), ts_Stim_all_sel(end)];
%     sts = unique(sel_stN, 'stable');
%     st_onoff_ts = [];
%     xticks = [];
%     hax(end).XTick = [];
%     for n = 1:length(sts)
%         xticks = cat(2, xticks, ts_Stim_all_sel( find(stN==sortedStimNs(n), 1) ) ); %trial divisor
%         st_onoff_ts = cat(2, st_onoff_ts, ts_Stim_all_sel(zeros_st_idxs(n)));
%         st_onoff_ts = cat(2, st_onoff_ts, ts_Stim_all_sel(zeros_st_idxs(n))+ T(sortedClassNs(n)).stimDur );
%         st_on = ts_Stim_all_sel(zeros_st_idxs(n));
%         st_off =ts_Stim_all_sel(zeros_st_idxs(n))+ stimDur;
%         plot([st_on, st_off, st_off, st_on, st_on], [-YLIM -YLIM YLIM YLIM -YLIM], 'Color', [0.7 0.7 0.7])
%     end
%     hax(end).XTick = unique(xticks);
%     hax(end).XMinorTick='off';
%     hax(end).XGrid = 'on';
%     hax(end).XAxis.Visible = 'off';
%     hax(end).YAxis.Visible = 'off';
%     
    
    for Ki = 1:nK
        axes(hax(Ki))
        hold off
        K2use = find(id_colors==Ki);
        allR = R_songs(Tparent==K2use,:);
        avg = mean(allR,1);
        plot([ts_dec_all_sel(1); ts_dec_all_sel; ts_dec_all_sel(end); ts_dec_all_sel(1)], [0, avg, 0, 0], 'Color',  cmap(K2use,:)), hold on, axis tight
        % plot([ts_dec_all_sel(1), ts_dec_all_sel(end)], [0,0], 'Color', colors(Ki,:) )
        
        hax(Ki).XLim =  [ts_Stim_all_sel(1), ts_Stim_all_sel(end)];
        hax(Ki).XAxis.Visible = 'off';
        hax(Ki).Box = 'off';
        
        %     ylims(Ki,:) = hax(Ki).YLim;
        drawnow
    end
    
    % ylims = [floor(min(ylims(:,1))/5)*5, ceil(max(ylims(:,2))/5)*5];
    
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
    export_fig(fullfile(folder2figure4, sprintf('resp_songs_%s.eps',num2str(sel_stN))))
%     
%     figure
%     spectrogram(stimTraces(i_St),10000,500,10000,4e4,'yaxis');    % used to be 20000 - 500 - 20000. Trying to distort less the temporal dimension.
%     ax = gca;
%     ax.TickDir = 'out';
%     ax.YLim = [0,0.6];
%     ax.YTick = [0, 0.3, 0.6];
%     ax.YTickLabel = 0:300:600;
%     ax.YLabel.String = 'Frequency (Hz)';
%     % ax.XTick = [];
%     % ax.XLabel.String = '';
%     ax.Box = 'off';
%     np = 128;
%     expn = 1.75;
%     cmap = brewermap(np, 'Oranges');
%     x = 1:np;
%     y = x.^expn;
%     y = ceil(y./max(y)*np);
%     cmap2 = interp1(x, cmap, y);
%     colormap(cmap2)
%     ax.XLim =  [ts_Stim_all_sel(1), ts_Stim_all_sel(end)];
%     export_fig(fullfile(folder2figure4, sprintf('songs_%s_spectrogram.eps',num2str(sel_stN))))
%     colorbar off
%     ax.Box = 'on';
%     ax.XTick = [];
%     ax.YTick = [];
%     ax.XLabel.Visible = 'off';
%     ax.YLabel.Visible = 'off';
%     ax.XAxis.Color = [1,1,1];
%     ax.YAxis.Color = [1,1,1];
%     export_fig(fullfile(folder2figure4, sprintf('songs_%s_spectrogram.jpg',num2str(sel_stN))), '-m4')
%     
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
% [R_steps, ~] = makeSaveR_mf_sortedStimuli( aZ, sets, pixel2keep, sortedStimNs, sortedClassNs, cropbaseline );
% save(saveName, 'R_steps', '-append') %will be updated
load(saveName, 'R_steps')

%% STEPS - plot
% timestamps!
[stimTraces, ts_Stim_all, stN, clN, zeros_st_idxs] = makeSortedStimulusSet_protot(T, sortedStimNs, sortedClassNs, cropbaseline);
stimTraces = stimTraces.*3; %um
YLIM = 15;

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
[hfig, hax] = figureTracesI_PP_pixels( nK, 200, 30 );
hax = flipud(hax); % 1 is bottom, as indendrogram
% plot stimulus
% axes(hax(end)), hold on
% plot(ts_Stim_all, stimTraces, '-k')
% hax(end).YLim = [-YLIM YLIM]; %um max range, fixed across stimuli
% hax(end).XLim =  [ts_Stim_all(1), ts_Stim_all(end)];
% sts = unique(stN, 'stable');
% st_onoff_ts = [];
% xticks = [];
% hax(end).XTick = [];
% for n = 1:length(sts)
%     xticks = cat(2, xticks, ts_Stim_all( find(stN==sortedStimNs(n), 1) ) ); %trial divisor
%     st_on = ts_Stim_all(zeros_st_idxs(n)); 
%     st_off =ts_Stim_all(zeros_st_idxs(n))+ T(sortedClassNs).stimDur;
%     plot([st_on, st_off, st_off, st_on, st_on], [-YLIM -YLIM YLIM YLIM -YLIM], 'Color', [0.7 0.7 0.7])
% end
% hax(end).XTick = unique(xticks);
% hax(end).XMinorTick='off';
% hax(end).XGrid = 'on';
% hax(end).XAxis.Visible = 'off';
% hax(end).YAxis.Visible = 'off';
% 

for Ki = 1:nK
    axes(hax(Ki))
    hold off
    K2use = find(id_colors==Ki);
    allR = R_steps(Tparent==K2use,:);
    avg = mean(allR,1);
%     cmap(K2use,:)
    plot([ts_dec_all(1); ts_dec_all; ts_dec_all(end); ts_dec_all(1)], [0, avg, 0, 0], 'Color',  cmap(K2use,:)), hold on, axis tight
    
    hax(Ki).XLim =  [ts_Stim_all(1), ts_Stim_all(end)];
    hax(Ki).XAxis.Visible = 'off';
    hax(Ki).Box = 'off';
    
%     ylims(Ki,:) = hax(Ki).YLim;
    drawnow
end

% ylims = [floor(min(ylims(:,1))/5)*5, ceil(max(ylims(:,2))/5)*5];
% ylims = [-50, 400];
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
export_fig(fullfile(folder2figure4, 'resp_steps.eps'))


% 
% figure
% spectrogram(stimTraces,10000,500,10000,4e4,'yaxis');    % used to be 20000 - 500 - 20000. Trying to distort less the temporal dimension.
% ax = gca;
% ax.TickDir = 'out';
% ax.YLim = [0,0.6];
% ax.YTick = [0, 0.3, 0.6];
% ax.YTickLabel = 0:300:600;
% ax.YLabel.String = 'Frequency (Hz)';
% ax.XTick = [];
% ax.XLabel.String = '';
% ax.Box = 'off';
% np = 128;
% expn = 1.75;
% cmap = brewermap(np, 'Oranges');
% x = 1:np;
% y = x.^expn;
% y = ceil(y./max(y)*np);
% cmap2 = interp1(x, cmap, y);
% colormap(cmap2)
% ax.XLim =  [ts_Stim_all(1), ts_Stim_all(end)];
% export_fig(fullfile(folder2figure4, 'steps_spectrogram.eps'))
% colorbar off
% ax.Box = 'on';
% ax.XTick = [];
% ax.YTick = [];
% ax.XLabel.Visible = 'off';
% ax.YLabel.Visible = 'off';
% ax.XAxis.Color = [1,1,1];
% ax.YAxis.Color = [1,1,1];
% export_fig(fullfile(folder2figure4, 'steps_spectrogram.jpg'), '-m4')
% 





%% chunk n: frequency-response curves at 3 amplitudes (pip data)
sizeMap = [60, 86];
clust = remapIndexedPointsToRun_RoisMaps(Tparent, pixel2keep, sizeMap);
ROIs = clust.ROIs; %these are still unsorted, so do not remove the last one!
CN = 'pips';
chunks.(CN)
% [timestamps, out, freqTunCurves] = trace_ChunkSubstims(classnum, trialTypes, aZ, ROIs) 
[chunks.(CN).timestamps, chunks.(CN).amplitude_1_out, chunks.(CN).amplitude_1_ft] = ...
    trace_ChunkSubstims(chunks.(CN).classnum, chunks.(CN).amplitude_1, aZ, ROIs);
%
[~, chunks.(CN).amplitude_2_out, chunks.(CN).amplitude_2_ft] = ...
    trace_ChunkSubstims(chunks.(CN).classnum, chunks.(CN).amplitude_2, aZ, ROIs);

[~, chunks.(CN).amplitude_3_out, chunks.(CN).amplitude_3_ft] = ...
    trace_ChunkSubstims(chunks.(CN).classnum, chunks.(CN).amplitude_3, aZ, ROIs);

R_pips = chunks.pips;
save(saveName, 'R_pips', '-append')
%% add stimulus (pip)
sortedStimNs = {42};
sortedClassNs = zeros(1, length(sortedStimNs));
for i = 1: length(sortedStimNs)
    sortedClassNs(i) = aZ{1}.stimuli2Class(sortedStimNs{i}(1), 1);
end
[stimTraces, ts_Stim_all, stN, clN, zeros_st_idxs] = makeSortedStimulusSet_protot(T, sortedStimNs, sortedClassNs, cropbaseline);
stimTraces = stimTraces.*3; %um
YLIM = 15;

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
[hfig, hax] = figureTracesI_PP_pixels( nK+1, 200, 30 );
hax = flipud(hax); % 1 is bottom, as indendrogram
% plot stimulus
axes(hax(end)), hold on
plot(ts_Stim_all, stimTraces, '-k')
hax(end).YLim = [-YLIM YLIM]; %um max range, fixed across stimuli
hax(end).XLim =  [ts_Stim_all(1), ts_Stim_all(end)];


resp_on = ts_Stim_all(zeros_st_idxs) + ts_dec(find(ts_dec>=0.25,1));
resp_off =ts_Stim_all(zeros_st_idxs) + ts_dec(find(ts_dec>=T(sortedClassNs).stimDur+0.25,1));
plot([resp_on, resp_off], [-10 -10], 'Color', [1 0 0])
bsl_on = 0;
bsl_off =ts_Stim_all(zeros_st_idxs)+ ts_dec(find(ts_dec>=0.12,1));
plot([bsl_on, bsl_off], [-10 -10], 'Color', [0 0 0])
% hax(end).XTick = unique(xticks);
% hax(end).XMinorTick='off';
% hax(end).XGrid = 'on';
hax(end).XAxis.Visible = 'off';
hax(end).YAxis.Visible = 'off';


XaxReference = min([floor(hax(end).XLim(end)/2)*2, 10]);
plot([0,XaxReference], [-YLIM -YLIM], '-k')
hax(end).XTick = [];
hax(end).XAxis.Label.String = sprintf('%d / %.2f',XaxReference, hax(end).XLim(end));
hax(end).XAxis.Label.Visible = 'on';
hax(end).XAxis.Label.FontSize = 9;

%% plot pips freq tuning

ylimFT = [];
for Ki = 1:nK
    axes(hax(Ki))
    hold off
    K2use = find(id_colors==Ki);
    xvar = chunks.(CN).carrierFreqLevels;
    yvar = squeeze(chunks.(CN).amplitude_1_ft.zsBasExt(K2use,:,:));
    plot(xvar, nanmean(yvar,2),'Color', cmap(K2use,:)*0.3), hold on
    
    xvar = chunks.(CN).carrierFreqLevels;
    yvar = squeeze(chunks.(CN).amplitude_2_ft.zsBasExt(K2use,:,:));
    plot(xvar, nanmean(yvar,2),'Color', cmap(K2use,:)*0.65), hold on
    
    xvar = chunks.(CN).carrierFreqLevels;
    yvar = squeeze(chunks.(CN).amplitude_3_ft.zsBasExt(K2use,:,:));
    plot(xvar, nanmean(yvar,2),'Color', cmap(K2use,:)), hold on
    
    plot([0,xvar(end)],[0 0], '-k')
    axis tight
    
    hax(Ki).Box = 'off';
    
    ylimsFT(Ki,:) = hax(Ki).YLim;
    set(gca, 'XTick', xvar)
    set(gca, 'XGrid', 'on')
    hax(Ki).TickDir = 'out';
    hax(Ki).XLim = [0,xvar(end)];
    if Ki == 1
        hax(Ki).XTickLabel = xvar;
        hax(Ki).XTickLabelRotation = 90;
        hax(Ki).XAxis.FontSize = 4;
    else
        hax(Ki).XTickLabel = [];
        hax(Ki).XAxis.Visible = 'off';
    end
    drawnow
end

ND = 1;
ylimsFT = [floor(min(ylimsFT(:,1))/ND)*ND, ceil(max(ylimsFT(:,2))/ND)*ND];
for Ki = 1:nK
    axes(hax(Ki))
    hax(Ki).YLim = ylimsFT;
    hax(Ki).YTick = [];
    hax(Ki).YAxis.Visible = 'off';
end
YaxReference = floor(diff(ylimsFT)/10)*10;
plot([hax(Ki).XLim(1),hax(Ki).XLim(1)], [ylimsFT(1), ylimsFT(1) + YaxReference], '-k')
hax(Ki).YAxis.Label.String = YaxReference;
hax(Ki).YAxis.Label.Visible = 'on';
hax(Ki).YAxis.Label.FontSize = 9;

% export and save!
hax(1).XAxis.FontSize = 4;
export_fig(fullfile(folder2figure4, 'resp_pips_ft.eps'))

%% spectrogram
sortedStimNs = {37:54};
sortedClassNs = zeros(1, length(sortedStimNs));
for i = 1: length(sortedStimNs)
    sortedClassNs(i) = aZ{1}.stimuli2Class(sortedStimNs{i}(1), 1);
end
[stimTraces, ts_Stim_all, stN, clN, zeros_st_idxs] = makeSortedStimulusSet_protot(T, sortedStimNs, sortedClassNs, cropbaseline);

figure
spectrogram(stimTraces,10000,500,10000,4e4,'yaxis');    % used to be 20000 - 500 - 20000. Trying to distort less the temporal dimension.
ax = gca;
ax.TickDir = 'out';
ax.YLim = [0,0.6];
ax.YTick = [0, 0.3, 0.6];
ax.YTickLabel = 0:300:600;
ax.YLabel.String = 'Frequency (Hz)';
ax.XTick = [];
ax.XLabel.String = '';
ax.Box = 'off';
np = 128;
expn = 1.75;
cmap = brewermap(np, 'Oranges');
x = 1:np;
y = x.^expn;
y = ceil(y./max(y)*np);
cmap2 = interp1(x, cmap, y);
colormap(cmap2)
export_fig(fullfile(folder2figure4, 'pips_spectrogram.eps'))
colorbar off
ax.Box = 'on';
ax.XTick = [];
ax.YTick = [];
ax.XLabel.Visible = 'off';
ax.YLabel.Visible = 'off';
ax.XAxis.Color = [1,1,1];
ax.YAxis.Color = [1,1,1];
export_fig(fullfile(folder2figure4, 'pips_spectrogram.jpg'), '-m4')



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
clust = remapIndexedPointsToRun_RoisMaps(Tparent, pixel2keep, sizeMap);
ROIs = clust.ROIs; %these are still unsorted, so do not remove the last one!

trialTypes = {chunks.pips.amplitude_1, chunks.pips.amplitude_2, chunks.pips.amplitude_3, ...
              chunks.tones.amplitude_1_carrPhase_1, chunks.tones.amplitude_2_carrPhase_1 }; %keep separate subsets
saveNameFTData = fullfile(clusterFolder, 'FTdata.mat');
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
export_fig(fullfile(folder2figure4, 'tuningCurves_01_06sec_4HzAndPips_MaxAmplit_100%.eps'))

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
export_fig(fullfile(folder2figure4, 'tuningAmplitude_01_06sec_100%.eps'))






