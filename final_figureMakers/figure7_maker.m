Folder2Save = '/Users/galileo/Dropbox (HMS)/Data/WINDS_100/Downsampled/SubsetFlies';
load(fullfile(Folder2Save, 'datalist_selectWindFlies.mat'));
load(fullfile(Folder2Save, 'R_matrix_downsampled_unfiltered_winds_100&86.mat'));    % unfiltered
load(fullfile(Folder2Save,'R_matrix_downsampled_winds_100&86_withBaseline.mat'), 'ts_dec');

R = Rflipped;
clear Rflipped
clear runfolders
for z = 1:length(datalist)
    [runfolders{z}, ~] = fileparts(datalist{z});
end

clusterFolder = fullfile(Folder2Save, '/rawDataLinkage'); % linkage was already calculated on orchestra
load(fullfile(clusterFolder, 'output_ward_Rzsp_linkage_winds_FLIPPED.mat'), 'Zc');
cd(clusterFolder)
pxKeep = true(size(R, 1),1);   % no (easily-detected at least) uncorrelated pixels within this dataset

load('/Users/galileo/Dropbox (HMS)/Data/cmap19_spectral.mat', 'cmapwinds')

%% select a K and plot subclusters (traces)
K = 7;
nK = K;
load(fullfile(clusterFolder, sprintf('klust_maxClust%d.mat', K)));   %'klust', 'dendrOrder', 'Tparent'); % these are sorted as in dedndrogram
orDT = load('/Users/galileo/Dropbox (HMS)/Data/WINDS_100/datalist_allWINDS.mat', 'dataTable');
orDT = orDT.dataTable;
clust = remapIndexedPointsToRun_RoisMaps_oldNew(Tparent, R_100, R_86);      %sukee
% make maps (divided old/new because different dimensions)
klust = makeplot_cMaps_superKlusts_singleIteration_oldNew(clust, klust, dendrOrder, cmap);
regGeneralData = matfile('/Users/galileo/Dropbox (HMS)/Data/WINDS_100/anatomyImages/anatomy_alignment_metadata.mat');
allStacks = regGeneralData.allStacks; %allStacks
allStacks_flyNums = cat(1,allStacks.flyNum);

dWA = matfile('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/datalist_WEDAMMC_piezo.mat');
alignedAMMCmasks = dWA.alignedAMMCmasks;

%% make indices for single flies
% flyOrdinals = zeros(size(pxKeep));
flyOrdinals = ones(sum(R_100.pixel2keep(:)), 1);
fly86Nums = [157, 179, 180];
for i = 1:3
    f = fly86Nums(i);
    zs = dataTable.fly == f;
    zs(1:13) = []; %remove f100
    flyOrdinals = cat(1, flyOrdinals, (i+1)*ones(sum(sum(R_86.pixel2keep(:,zs))),1) );
end

%%
stim_boxcar = ts_dec>=0 & ts_dec<=4;
stim_ipsi = cat(1, stim_boxcar, stim_boxcar, -stim_boxcar, -stim_boxcar);
stim_contra = cat(1, -stim_boxcar, stim_boxcar, stim_boxcar, -stim_boxcar);
%
[hfig, hax] = figureTracesI_PP_pixels( 6, 250, 60 );
hax = flipud(hax);
ks = [3 5 1 7];
for iKi = 1:4
    Ki = ks(iKi);
    % parent cluster
    axes(hax(iKi,1)), hold on
    allR = R(Tparent==Ki,:);
    avg = mean(allR,1);

    plot(avg)
    avg = medfilt1(avg,11);
    plot(avg)
    plot([1,340], [0 0], '-k')
    hax(iKi,1).XAxis.Visible = 'off';
end
axes(hax(iKi+2,1)), hold on
plot(stim_ipsi, '-k')
axes(hax(iKi+1,1)), hold on
plot(stim_contra, '-k')
hax(iKi+2,1).XAxis.Visible = 'off';
hax(iKi+1,1).XAxis.Visible = 'off';

hax(1,1).XLim = [1,340];
hax(1,1).YLim = [-20,60];
linkaxes(hax(1:iKi), 'xy')
export_fig(sprintf('clust%2d_responses_baselines_medfilt11_concatenated_sorted.eps',K))
% don't look so well because of cropped dynamics

%% bar plots with average dff in sustained phase only - no single fly
% stim sorted as before, I'll redesign it in AI
assert(size(R, 2) == 340, 'the following indices re meant for R with baseline')
sustainedIdxs = 24+[19:41];

[hfig, hax] = figureTracesI_PP_pixels( 4, 250, 60 );
hax = flipud(hax);
ks = [3 5 1 7];
for iKi = 1:4
    Ki = ks(iKi);
    allR = R(Tparent==Ki,:);
    avg = mean(allR,1);
    avg = reshape(avg(:), [],4);
    avg = mean(avg(sustainedIdxs,:), 1);
    bar(hax(iKi,1), avg, 'FaceColor', 'none');
    hax(iKi,1).XAxis.Visible = 'off';
end

hax(1,1).YLim = [-20,40];
linkaxes(hax(1:iKi), 'xy')
export_fig(sprintf('clust%2d_responses_bars_concatenated_sorted.eps',K))

%% bar plots with average dff in sustained phase only - single fly separation - histogram
% stim sorted as before, I'll redesign it in AI
assert(size(R, 2) == 340, 'the following indices re meant for R with baseline')
sustainedIdxs = 24+[19:41];

[hfig, hax] = figureTracesI_PP_pixels( 4, 250, 60 );
hax = flipud(hax);
ks = [3 5 1 7];
for iKi = 1:4
    Ki = ks(iKi);
    % split flies
    clear avgflies
    for f = 1:4
        allR = R(Tparent==Ki & flyOrdinals==f,:);
        avgflies(1,:,f) = mean(allR,1);
    end
    % reshape 4 stimuli
    avgflies = reshape(avgflies, [], 4, 4); %time,stim,fly
    avgflies = mean(avgflies(sustainedIdxs,:,:), 1);
    avgflies = squeeze(avgflies); %stim,fly
    
    
    allR = R(Tparent==Ki,:);
    avg = mean(allR,1);
    avg = reshape(avg(:), [],4);
    avg = mean(avg(sustainedIdxs,:), 1);
    bar(hax(iKi,1), avg, 'FaceColor', 'none');
    axes(hax(iKi,1))
    hold on
    for f = 1:4
        plot(1:4, avgflies(:,f))
    end
    
%     bar(hax(iKi,1), avg, 'FaceColor', 'none');
%     hax(iKi,1).XAxis.Visible = 'off';
end

hax(1,1).YLim = [-25,50];
linkaxes(hax(1:iKi), 'xy')
export_fig(sprintf('clust%2d_responses_bars_concatenated_sorted_singleFliesLINES.eps',K))



%% single bar (no fly division), plus lines for flies
assert(size(R, 2) == 340, 'the following indices re meant for R with baseline')
sustainedIdxs = 24+[19:41];

[hfig, hax] = figureTracesI_PP_pixels( 4, 250, 60 );
hax = flipud(hax);
ks = [3 5 1 7];
for iKi = 1:4
    Ki = ks(iKi);
    % split flies
    clear avg
    for f = 1:4
        allR = R(Tparent==Ki & flyOrdinals==f,:);
        avg(1,:,f) = mean(allR,1);
    end
    % reshape 4 stimuli
    avg = reshape(avg, [], 4, 4); %time,stim,fly
    avg = mean(avg(sustainedIdxs,:,:), 1);
    avg = squeeze(avg); %stim,fly
    
    bar(hax(iKi,1), avg, 'FaceColor', 'none');
    hax(iKi,1).XAxis.Visible = 'off';
end

hax(1,1).YLim = [-25,50];
linkaxes(hax(1:iKi), 'xy')
export_fig(sprintf('clust%2d_responses_bars_concatenated_sorted_singleFlies.eps',K))









%% plot unconcatenated responses (unfiltered traces)
[hfig, hax] = figureTracesI_PP_pixels( 4, 120*ones(1,4), 60 );
hax = flipud(hax);
ks = [3 5 1 7];
for iKi = 1:4
    Ki = ks(iKi);
    allR = R(Tparent==Ki,:);
    avg = mean(allR,1);
    avg = reshape(avg(:), [],4);
    for i = 1:4
        axes(hax(iKi,i)), hold on
        plot([0 0], [-20, 60], '-k')
        plot([0 4], [0, 0], '-k')
        plot(ts_dec, avg(:,i))
        axis off
    end
end

hax(1,1).XLim = [ts_dec(1), 7.4920]; %5 times the sustained phase duration
hax(1,1).YLim = [-20, 60];
linkaxes(hax, 'xy')

axes(hax(1,1))
plot(ts_dec([24+19, 24+41]), [20, 20], '-k')
plot([ts_dec(1), 7.4920],[-20 -20], '-k')
export_fig(sprintf('clust%2d_responses_unconcatenated_sorted.eps',K))

%%
close