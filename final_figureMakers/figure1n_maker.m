%% load
%
warning('this has only in small part adapted for re-sorted the 10 clusters (bar plot)')
%

load_panSpec_downsampled;
load('/Users/galileo/Dropbox (HMS)/Data/chunks_stimfinal.mat', 'chunks')
load('/Users/galileo/Dropbox (HMS)/Data_Raw_Metadata_BackedUp/BackedUp_Metadata/fly171_run01_metadata_2016-12-20_163728.mat', 'stimuli')
folder2stimuli = '/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/stimuli';
folder2figure1 = '/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/figure1';

%%
load('/Users/galileo/Dropbox (HMS)/Data/cmap19_spectral.mat', 'cmap10')

iRun = 14;
aZ = matfile(datalist_allJON{iRun}); %174_run01
T = aZ.T;
fastStimMeta = aZ.fastStimMeta;
clear aZ
for z = 1:2% NaZ
    aZ{z} = matfile(datalist_allJON{z});
end

load(fullfile(Folder2Save, 'R_matrix_Downsampled_smallWindow.mat'), 'R', 'pixel2keep', 'iKeep', 'pxKeep', 'oldpixel2keep') %pxKeep have already been selected in the pixel2keep variable
R = R(pxKeep,:);

clusterFolder = '/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix';
treeFileName = 'output_ward_movAVG_ZSCpix_linkage.mat';
cd(clusterFolder)
load('klust_maxClust10.mat')
nK = length(klust);


folder2stimuli = '/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/figure2n';
% folder2figure1 = '/Users/galileo/Dropbox (HMS)/figures/figure1';
folder2figure2 = '/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/figure2n';
% folder2figure3 = '/Users/galileo/Dropbox (HMS)/figures/figure3';
cropbaseline = - 1.7; %first point give errors
saveName = fullfile(Folder2Save, 'R_chunks.mat');
ylims = [-50, 400];

clusterOrder = [4 1 2 3 7 5 6 8 9 10]; % this is K2use!! bottom-up
cOrder2 = [1 2 3 4 8 5 7 6 9 10];
[~, cmapOrder] = sort(cOrder2);
cmap10 = cmap10(cmapOrder, :);      %color resorted based on cdf
cmap10 = cmap10([4:end, 3 2 1], :);     %index with Ki
figure; imagesc(1:size(cmap10,1)); colormap(cmap10); title('resorted cmap10')


%% pick pixels 
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
Ki = 2;
K2use = clusterOrder(Ki);
R_i = find(Tparent==K2use);
responses = R(R_i,:);
maxResp = max(responses,[],2);
selectResp_i = find(maxResp>(mean(maxResp)+3*std(maxResp)));
selectResp_i = selectResp_i(randperm(length(selectResp_i)));

Nuse = 20;
[hfig, hax] = figureTracesI_PP( Nuse, 1 );
for i = 1:Nuse
    axes(hax(i,1));
    plot(responses(selectResp_i(i),:));
    ylabel(sprintf('%d\n(%d)', R_i(selectResp_i(i)), selectResp_i(i) ))
    axis tight
end
hfig.Position = [0.36    0.0037    0.5510    0.9083];
%%
[hfig, hax] = figureTracesI_PP( 8, 1 );

Ki = 4;
K2use = clusterOrder(Ki);
selPixels = [72449 41937 42025 30800];
selTraces = Rsubset.R(selPixels,:);
for i = 1:4
    axes(hax(i,1));
    plot(ts_dec_all, selTraces(i,:), 'Color', cmap10(Ki,:));
    hold on
    plot([ts_dec_all(1), ts_dec_all(end)], [0 0], 'k')
    plot([0 0], [0 300], '-k')
    hax(i,1).YAxis.Visible = 'off';
    hax(i,1).XMinorGrid = 'off';
    hax(i,1).XGrid = 'on';
    hax(i,1).XTick = ts_dec_all(zero_tr_idxs(majorTicks));
    hax(i,1).XAxis.TickLength = [0 00];
    axis tight
end

Ki = 1;
K2use = clusterOrder(Ki);
selPixels = [65205, 72708, 60529, 20662];
selTraces = Rsubset.R(selPixels,:);
for i = 1:4
    axes(hax(i+4,1));
    plot(ts_dec_all, selTraces(i,:), 'Color', cmap10(Ki,:));
    hold on
    plot([ts_dec_all(1), ts_dec_all(end)], [0 0], 'k')
    plot([0 0], [0 300], '-k')
    hax(i+4,1).YAxis.Visible = 'off';
    hax(i+4,1).XMinorGrid = 'off';
    hax(i+4,1).XGrid = 'on';
    hax(i+4,1).XTick = ts_dec_all(zero_tr_idxs(majorTicks));
    hax(i+4,1).XAxis.TickLength = [0 00];
    axis tight
end
export_fig(fullfile(folder2figure1, sprintf('example_NEWtraces_JONcl4_1.eps')))


%% there was a pproblem with R_songs' ts


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
% cd(folder2stimuli)

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
export_fig('stimFinal_stimulusSet_sorted_cleanedUp.eps')
%% plot spectrogram
f = figure;
spectrogram(stimTraces,10000,500,10000,4e4,'yaxis');
ax = gca;
ax.TickDir = 'out';
ax.YLim = [0,0.61];
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

%%
np = 128;
expn = 3;
killpt = 15;
cmap = brewermap(np, 'Greys');
x = 1:np;
y = x.^expn;
y = ceil(y./max(y)*np);
cmap2 = interp1(x, cmap, y);
cmap2(1:killpt,:) = [];
cmap2= cat(1, cmap2, zeros(killpt,3));
% figure; imagesc(1:size(cmap2,1)); colormap(cmap2);
figure(f)
colormap(cmap2)
axis on

%%

export_fig(sprintf('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/figure1/stimFinal_stimulusSet_sorted_cleanedUp_spectrogram.tiff'))
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


