%% load
%
warning('this has only in small part adapted for re-sorted the 10 clusters (bar plot)')
%

load_panSpec_downsampled;
load('/Users/galileo/Dropbox (HMS)/Data/chunks_stimfinal.mat', 'chunks')
load('/Users/galileo/Dropbox (HMS)/Data_Raw_Metadata_BackedUp/BackedUp_Metadata/fly171_run01_metadata_2016-12-20_163728.mat', 'stimuli')

% % this is currently unsynched. Do I really need the full aZ??
% %load aZ
% clear aZ
% for z = 1:NaZ
%     aZ{z} = matfile(datalist{z});
% end


folder2stimuli = '/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/stimuli';
folder2figure1 = '/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/figure1';

%% GENERAL cmap _ 20170913
load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/Downsampled/rawDataLinkage/ward_ZSCpix/klust_maxClust19.mat', 'cmap')
figure; imagesc(1:size(cmap,1)); colormap(cmap)
ax = gca; ax.YTick = []; ax.XTick = 1:19;
cmap = makeRainbowCmap;

%% sample dataset?
iRun = 14;
if ~exist('aZ', 'var')
    aZ = matfile(datalist{iRun});
    T = aZ.T;
    fastStimMeta = aZ.fastStimMeta;
else
    T = aZ{iRun}.T;
    fastStimMeta = aZ{iRun}.fastStimMeta;
end

%% new fuller dendrogram _ 20170913
load('R_matrix_Downsampled_smallWindow.mat', 'R', 'iKeep', 'pixel2keep', 'pxKeep') %pxKeep have already been selected in the pixel2keep variable
R = R(pxKeep,:);
clusterFolder = '/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix';
treeFileName = 'output_ward_movAVG_ZSCpix_linkage.mat';
cd(clusterFolder)
load('klust_maxClust11.mat')
load(fullfile(clusterFolder, treeFileName), 'Zc');
maxClust = 11;
maxLeavesDend = maxClust;
order11 = [7 8 9 6 5 4 3 10 11 1 2]; %this is actually flipped horizontally
outpermLeaves11 = [5    11     1    10     4     2     3     6     7     8     9];
figure; [~,Tdend,outpermLeaves] = dendrogram(Zc, maxLeavesDend, 'Orientation', 'left');
for i = 1:maxClust
    count_peroutpermLeaves(i) = sum(Tdend==outpermLeaves(i));
end
for i = 1:maxClust
    count_perTParent(i) = sum(Tparent==i);
end  
% the saved dendrogram figure has been resorted accordingly manually, 
% and these would be the corresponding outpermLeaves:
outpermLeaves = fliplr(outpermLeaves11(order11)); % I am not changing the cluster identity here, just moving branches around (physically in the dendrogram)
% so far, this was old stuff. So, there is still correspondence between
% outpermLeaves and Tdend

% NOTE that outpermLeaves11(fliplr(order11)) == fliplr(outpermLeaves11(order11))

% NOTE NOTE: Tparent had been resorted according to dendrogram.


% make a dendrogram with far more leaves than the clusters, though not complete
% showing ~5% of the tree leaves. With 10% the structure does not change at all.
figure; [~,TdendBig,outpermLeavesBig] = dendrogram(Zc, 5134, 'Orientation', 'left'); %, 'ColorThreshold', maxClust); % colorthreshold does not work.
ax = gca;
ax.YTick = [];
export_fig(fullfile(folder2figure1, 'dendrogram_5percentLeaves.eps'))

% map onto existing clusters?
for Ki = 1:length(unique(Tparent)), allElements(Ki, 1:length(klust(Ki).k)) = klust(Ki).k; end
% try again, intuitively -- both ways yield the same result
for countLeaves = 1:length(outpermLeavesBig)
    leaf = outpermLeavesBig(countLeaves);
    table_T01 = crosstab(TdendBig==leaf, Tparent);
    clustNames_allLeavesBig(countLeaves) = find(table_T01(2,:));
end
clustNames = unique(clustNames_allLeavesBig, 'stable');
% 
% countLeaves = 1;
% for Ki = 1:length(unique(Tparent))   %while ~isempty(leavesorder)
%     lengthLeave = 0;
%     leavesorder = find(TdendBig == outpermLeavesBig(countLeaves));
%     countLeaves = countLeaves+1;
%     lengthLeave = lengthLeave + length(leavesorder);
%     firstLeaf(Ki) = leavesorder(1);
%     [dendrOrder(Ki), ~] = find(allElements==firstLeaf(Ki));
%     while sum(Tparent == dendrOrder(Ki)) > lengthLeave
%         leavesorder = find(TdendBig == outpermLeavesBig(countLeaves));
%         countLeaves = countLeaves+1;
%         lengthLeave = lengthLeave + length(leavesorder);
%     end
% end
% % clear leavesorder countLeaves leavesorder

%% pick pixels from clusters 6 and 9 (3-4 each), and make mock dendrogram
Rsubset = load(fullfile(folder2figure1,'R_subsetStim_figure1.mat'), 'R', 'sortedClassNs', 'sortedStimNs', 'cropbaseline');
sortedClassNs = Rsubset.sortedClassNs;
sortedStimNs = Rsubset.sortedStimNs;
cropbaseline = Rsubset.cropbaseline;
ts_dec_all_maker;

for c = 1:length(sortedClassNs)
    majorTicks(c) = length(sortedStimNs{c});
end
majorTicks = cumsum(majorTicks) +1;
majorTicks(1) = 1;
majorTicks(end) = [];
%%
Ki = 5;
R_i = find(Tparent==Ki);
responses = Rsubset.R(R_i,:);
maxResp = max(responses,[],2);
selectResp_i = find(maxResp>mean(maxResp));
selectResp_i = selectResp_i(randperm(length(selectResp_i)));
Nuse = 10;
[hfig, hax] = figureTracesI_PP( Nuse, 1 );
for i = 1:Nuse
    axes(hax(i,1));
    plot(responses(selectResp_i(i),:));
    ylabel(sprintf('%d\n(%d)', R_i(selectResp_i(i)), selectResp_i(i) ))
    axis tight
end
hfig.Position = [0.36    0.0037    0.5510    0.9083];
export_fig(fullfile(folder2figure1, sprintf('example_traces_JONcluster%d_03.eps',Ki)))

%
Ki = 11;
R_i = find(Tparent==Ki);
responses = Rsubset.R(R_i,:);
maxResp = max(responses,[],2);
selectResp_i = find(maxResp>mean(maxResp));
selectResp_i = selectResp_i(randperm(length(selectResp_i)));
Nuse = 10;
[hfig, hax] = figureTracesI_PP( Nuse, 1 );
for i = 1:Nuse
    axes(hax(i,1));
    plot(responses(selectResp_i(i),:));
    ylabel(sprintf('%d\n(%d)', R_i(selectResp_i(i)), selectResp_i(i) ))
    axis tight
end
hfig.Position = [0.36    0.0037    0.5510    0.9083];
export_fig(fullfile(folder2figure1, sprintf('example_traces_JONcluster%d_03.eps',Ki)))

%% select pixels
selPixels = [65205, 72708, 58475, 60529, 20662, 35862, 6110, 72358, 47531, 51387];
selPixels_supervisedLabel = [5*ones(1,5), 11*ones(1,5)];
selTraces = Rsubset.R(selPixels,:);
selTraces_zs = zscore(selTraces, 0, 2);
% compute linkage
Zc = linkage(selTraces_zs, 'ward', 'euclidean');
figure; [~,Tdend,outpermLeaves] = dendrogram(Zc, 'Orientation', 'left');
export_fig(fullfile(folder2figure1, 'selectedPixels_clusts05_11_dendrogram.eps'))

%% plot corresponding traces
[hfig, hax] = figureTracesI_PP( length(selPixels), 1 );
hax = flipud(hax);
for i = 1:5
    axes(hax(i,1));
    pixelUse = outpermLeaves(i);
    plot([ts_dec_all(1); ts_dec_all; ts_dec_all(end); ts_dec_all(1)], [0, selTraces(pixelUse,:), 0, 0], 'Color', colors(11,:)), hold on, axis tight    
    hold on
    plot([0 0], [0 300], '-k')
    hax(i,1).YAxis.Visible = 'off';
    hax(i,1).XMinorGrid = 'on';
    hax(i,1).XGrid = 'on';
    hax(i,1).XTick = ts_dec_all(zero_tr_idxs(majorTicks));
    hax(i,1).XAxis.MinorTickValues = ts_dec_all(zero_tr_idxs);
    hax(i,1).XAxis.TickLength = [0 00];
    ylabel(sprintf('%d\n%d', pixelUse, selPixels(pixelUse) ))
    axis tight
end
for i = 6:10
    axes(hax(i,1));
    pixelUse = outpermLeaves(i);
    plot([ts_dec_all(1); ts_dec_all; ts_dec_all(end); ts_dec_all(1)], [0, selTraces(pixelUse,:), 0, 0], 'Color', colors(5,:)), hold on, axis tight    
    hold on
    plot([0 0], [0 300], '-k')
    hax(i,1).YAxis.Visible = 'off';
    hax(i,1).XMinorGrid = 'on';
    hax(i,1).XGrid = 'on';
    hax(i,1).XTick = ts_dec_all(zero_tr_idxs(majorTicks));
    hax(i,1).XAxis.MinorTickValues = ts_dec_all(zero_tr_idxs);
    hax(i,1).XAxis.TickLength = [0 00];
    ylabel(sprintf('%d\n%d', pixelUse, selPixels(pixelUse) ))
    axis tight
end
hfig.Position = [0.36    0.0037    0.5510    0.9083];


export_fig(fullfile(folder2figure1, 'selectedPixels_clusts05_11_traces.eps'))

%% plot corresponding traces - NO TIMESTAMPS
[hfig, hax] = figureTracesI_PP( length(selPixels), 1 );
hax = flipud(hax);
for i = 1:5
    axes(hax(i,1));
    pixelUse = outpermLeaves(i);
    plot(selTraces(pixelUse,:), 'Color', colors(11,:)), hold on, axis tight
    plot([1, size(selTraces,2)], [0 0], '-k')
%     plot([1, 1:size(selTraces,2), size(selTraces,2), 1], [0, selTraces(pixelUse,:), 0, 0], 'Color', colors(11,:)), hold on, axis tight    
    plot([0 0], [0 300], '-k')
    hax(i,1).YAxis.Visible = 'off';
    hax(i,1).XAxis.Visible = 'off';
%     hax(i,1).XMinorGrid = 'on';
%     hax(i,1).XGrid = 'on';
%     hax(i,1).XTick = ts_dec_all(zero_tr_idxs(majorTicks));
%     hax(i,1).XAxis.MinorTickValues = ts_dec_all(zero_tr_idxs);
%     hax(i,1).XAxis.TickLength = [0 00];
%     ylabel(sprintf('%d\n%d', pixelUse, selPixels(pixelUse) ))
%     axis tight
end
for i = 6:10
    axes(hax(i,1));
    pixelUse = outpermLeaves(i);
%     plot([1, 1:size(selTraces,2), size(selTraces,2), 1], [0, selTraces(pixelUse,:), 0, 0], 'Color', colors(5,:)), hold on, axis tight    
    plot(selTraces(pixelUse,:), 'Color', colors(5,:)), hold on, axis tight 
    plot([1, size(selTraces,2)], [0 0], '-k')
    plot([0 0], [0 300], '-k')
    hax(i,1).YAxis.Visible = 'off';
    hax(i,1).XAxis.Visible = 'off';
%     hax(i,1).XMinorGrid = 'on';
%     hax(i,1).XGrid = 'on';
%     hax(i,1).XTick = ts_dec_all(zero_tr_idxs(majorTicks));
%     hax(i,1).XAxis.MinorTickValues = ts_dec_all(zero_tr_idxs);
%     hax(i,1).XAxis.TickLength = [0 00];
%     ylabel(sprintf('%d\n%d', pixelUse, selPixels(pixelUse) ))
%     axis tight
end
hfig.Position = [0.36    0.0037    0.5510    0.9083];


export_fig(fullfile(folder2figure1, 'selectedPixels_clusts05_11_traces_notimestamps.eps'))


%% stim parameters (sample set)
%Chirp_down, Chirp_up, Tones_ampl2 (excluding phase), pips ampls 1 - 2 - 3,
%songA, songB, steps]
cropbaseline = -1;
sortedStimNs = {62, 61, [86 90 94 76 78 80 82 84], 1:18, 19:36, 37:54, 74, 73, 55:60 };
sortedClassNs = zeros(1, length(sortedStimNs));
for i = 1: length(sortedStimNs)
%     sortedClassNs(i) = aZ{1}.stimuli2Class(sortedStimNs{i}(1), 1);
    sortedClassNs(i) = aZ.stimuli2Class(sortedStimNs{i}(1), 1);
end

%% concatenated stimulus set [may skip]
cd(folder2stimuli)

% make prototype stim set
[stimTraces, ts_Stim_all, stN, clN, zeros_st_idxs] = makeSortedStimulusSet_protot(T, sortedStimNs, sortedClassNs, cropbaseline);
figure;
plot(ts_Stim_all, stimTraces.*3, '-k') %um
ax = gca;
ax.XLim = [ts_Stim_all(1), ts_Stim_all(end)];
ax.XTick = ts_Stim_all(zeros_st_idxs);
ax.XTickLabel =[];
ax.XGrid = 'on';
ax.TickDir = 'out';
ax.Box = 'off';
hold on;
plot([100,110], [10,10], '-k') %10 seconds bar
export_fig('stimFinal_stimulusSet_soretd_cleanedUp.eps')
% plot spectrogram
figure
spectrogram(stimTraces,10000,500,10000,4e4,'yaxis');
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
export_fig(sprintf('stimFinal_stimulusSet_soretd_cleanedUp_spectrogram.eps'))
clear x y cmap cmap2 expn i clN stN

%% make R for single pixels traces - selected stimulus set
cd(Folder2Save)
if ~exist('pixel2keep', 'var')
    load('R_matrix_Downsampled_smallWindow.mat', 'pixel2keep', 'iKeep', 'pxKeep', 'oldpixel2keep')
    if sum(pixel2keep(:)) > sum(pxKeep)
        % una tantum - DONE
        oldpixel2keep = pixel2keep;
        sizePXK = size(pixel2keep);
        pixel2keep = pixel2keep(:);
        pixel2keep(pixel2keep) = pxKeep;
        pixel2keep = reshape(pixel2keep, sizePXK);
        save('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/R_matrix_Downsampled_smallWindow.mat', 'pixel2keep', 'oldpixel2keep', '-append');
    end
end
%remake R
R = makeSaveR_mf_sortedStimuli( aZ, sets, pixel2keep, sortedStimNs, sortedClassNs, cropbaseline );
save(fullfile(folder2figure1,'R_subsetStim_figure1.mat'), 'R', 'sortedClassNs', 'sortedStimNs', 'cropbaseline')

%% sample run image -> registered
iRun = 14;






%% choose pixels 2 show (version 1)
x = [75 64 43 65 69 58 18 12 26 17];
y = [25 45 32 15 16 18 17 20 20 28 ];
clear px px_R traces
close all
pxk_iRun = pixel2keep(:,iRun);
pxk_iRun_withNoise = oldpixel2keep(:,iRun);
for i = 1:length(x)
    px(i) = (x(i)-1)*sizeMap(1) + y(i);
    px_R(i) = sum(pxk_iRun(1:px(i)));
    traces(i,:) = R_iRun(px_R(i),:);
    figure;
    plot(traces(i,:), '-k')
    title(sprintf('%d %d',x(i), y(i)))
    axis tight
end
traces = zscore(traces, 0, 2);
keep = [1,[1,2,4,5,7]+1];
colors = [0 0 0; 0.5 0.5 0.5; 1 1 1; ...
          brewermap(5, 'Set1')];
map = uint8(pxk_iRun_withNoise);
map(pxk_iRun) = 2;
for i = 1:length(keep)
    map(px(keep(i))) = i+2;
end
map = reshape(map, sizeMap);
map = uint8(map);
figure; image(map); colormap(colors); axis image; axis off; colorbar off
export_fig(fullfile(folder2figure1, sprintf('selected_pixels_%s.tif', basenames{iRun})), '-m3')
colorbar

%\% replot and save corresponding traces
clear ax
i = 1; %'noise pixel'
figure;
plot(traces(keep(i),:), 'Color', colors(2,:))
axis tight
ax(i) = gca;

for i = 2:length(keep)
    figure;
    plot(traces(keep(i),:), 'Color', colors(i+2,:))
    axis tight
    ax(i) = gca;
    ax(i).Box = 'off';
    hold on
end
linkaxes([ax(2:length(keep)), ax(1)], 'xy')
for i = 1:length(keep)
    export_fig(ax(i), fullfile(folder2figure1, sprintf('zsc_trace_%d.eps',i)))
end







%% load cluster data and replot histograms
Rzs = zscore(R, 0, 2);
clusterFolder = '/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix';
treeFileName = 'output_ward_movAVG_ZSCpix_linkage.mat';
if exist(clusterFolder) ~= 7
    mkdir(clusterFolder)
end
cd(clusterFolder)
load('klust_maxClust11.mat')
load('alignedMatrix_Images_maxClust11.mat')
colors = Col(2:end,:);
assert(size(R,1) == sum(pxKeep))


%% variante 1: all zs per fly averaged
% first, get clust-indexed aligned Maps - sizeMap will be bigger...
clust30 = remapAlignedMapsToRun_RoisMaps(AllClusterMaps);
% clust30 = remapIndexedPointsToRun_RoisMaps(Tparent, pixel2keep, sizeMap, pxKeep);
sortFliesGenoType = 1;
[~, flyClusterCountDivZs30] = barcluster_panSpec(Tparent, datalist, flyfoldersUnique, clust30.ROIs, sortFliesGenoType);

nK = length(klust);
sortedflyMeta = {'pan' 'pan' 'pan' 'pan' 'pan' 'pan' 'A22' 'A23' 'A23' 'A26' 'AB15' 'AB15' 'AB15' 'BA28' 'B2' 'B2' 'ACE4' 'ACE4' 'ACE4' 'CE32' 'CE32' 'AD29'  'AD29' 'AD29'};
[hfig, hax] = figureTracesI_PP_pixels( nK, [300], 30 );
hax = flipud(hax) % 1 is bottom, as indendrogram
for k = 1:nK
    axes(hax(k))
    bar(1:NaF, flyClusterCountDivZs30(:,k), 'FaceColor', colors(k,:))
    hax(k).XTick = 1:NaF;
    axis tight
    if k == 1
        hax(k).FontSize = 7;
        hax(k).XTickLabel = flyNumUnique(idx_sortingFlies);
        hax(k).XTickLabel = sortedflyMeta;
        hax(k).XTickLabelRotation = 90;
    else
        hax(k).XTickLabel = [];
    end
end
export_fig(fullfile(folder2figure1, 'clusterBars_AlignedMaps_CountOverAllRunsAveraged.eps'))

%% variante 2 - single max Z chosen per cluster, per fly
% first, get clust-indexed aligned Maps - sizeMap will be bigger...
clust30 = remapAlignedMapsToRun_RoisMaps(AllClusterMaps);
% clust30 = remapIndexedPointsToRun_RoisMaps(Tparent, pixel2keep, sizeMap, pxKeep);
sortFliesGenoType = 1;
flyClusterCountDivZs30_maxZout = barcluster_panSpec_maxZout(Tparent, datalist, flyfoldersUnique, clust30.ROIs, sortFliesGenoType);

nK = length(klust);
sortedflyMeta = {'pan' 'pan' 'pan' 'pan' 'pan' 'pan' 'A22' 'A23' 'A23' 'A26' 'AB15' 'AB15' 'AB15' 'BA28' 'B2' 'B2' 'ACE4' 'ACE4' 'ACE4' 'CE32' 'CE32' 'AD29'  'AD29' 'AD29'};
[hfig, hax] = figureTracesI_PP_pixels( nK, [300], 30 );
hax = flipud(hax); % 1 is bottom, as indendrogram
for k = 1:nK
    axes(hax(k))
    bar(1:NaF, flyClusterCountDivZs30_maxZout(:,k), 'FaceColor', colors(k,:));
    hax(k).XTick = 1:NaF;
    axis tight
    if k == 1
        hax(k).FontSize = 7;
        hax(k).XTickLabel = flyNumUnique(idx_sortingFlies);
        hax(k).XTickLabel = sortedflyMeta;
        hax(k).XTickLabelRotation = 90;
    else
        hax(k).XTickLabel = [];
    end
end
export_fig(fullfile(folder2figure1, 'clusterBars_AlignedMaps_CountOverMaxZ.eps'))

%% variante 3 - single max Z chosen per cluster, per fly, over cropped area of aligned matrix
% first, get clust-indexed aligned Maps - sizeMap will be bigger...
figure;
[~, rC]  = imcrop(mean(AllClusterMaps>0, 3));
rC = round(rC);
AllClusterMaps_crop = AllClusterMaps(rC(2):rC(2)+rC(4), rC(1):rC(1)+rC(3), : );

hold on
rectangle('Position', rC, 'EdgeColor', 'r');
export_fig(fullfile(folder2figure1, 'crop_Allmaps.eps'), '-m4')

clust30 = remapAlignedMapsToRun_RoisMaps(AllClusterMaps_crop);
% clust30 = remapIndexedPointsToRun_RoisMaps(Tparent, pixel2keep, sizeMap, pxKeep);
sortFliesGenoType = 1;
flyClusterCountDivZs30_maxZout = barcluster_panSpec_maxZout(Tparent, datalist, flyfoldersUnique, clust30.ROIs, sortFliesGenoType);

nK = length(klust);
sortedflyMeta = {'pan' 'pan' 'pan' 'pan' 'pan' 'pan' 'A22' 'A23' 'A23' 'A26' 'AB15' 'AB15' 'AB15' 'BA28' 'B2' 'B2' 'ACE4' 'ACE4' 'ACE4' 'CE32' 'CE32' 'AD29'  'AD29' 'AD29'};
[hfig, hax] = figureTracesI_PP_pixels( nK, [300], 30 );
hax = flipud(hax); % 1 is bottom, as indendrogram
for k = 1:nK
    axes(hax(k))
    K2use = clusterOrder(Ki);
    bar(1:NaF, flyClusterCountDivZs30_maxZout(:,k), 'FaceColor', cmap10(k,:));
    hax(k).XTick = 1:NaF;
    axis tight
    if k == 1
        hax(k).FontSize = 7;
        hax(k).XTickLabel = flyNumUnique(idx_sortingFlies);
        hax(k).XTickLabel = sortedflyMeta;
        hax(k).XTickLabelRotation = 90;
    else
        hax(k).XTickLabel = [];
    end
end
export_fig(fullfile(folder2figure1, '10K_clusterBars_AlignedMaps_CountOverMaxZ_croppedMap.eps'))


