% load
% based on: clustering_on_rawData_selectWINDS_FLIPPED
% ---> in folder: analysis/final_WINDS
Folder2Save = '/Users/galileo/Dropbox (HMS)/Data/WINDS_100/Downsampled/SubsetFlies';
load(fullfile(Folder2Save, 'datalist_selectWindFlies.mat'));
load(fullfile(Folder2Save, 'R_matrix_downsampled_winds_100&86_withBaseline.mat'));
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
K = 8;
nK = K;
load(fullfile(clusterFolder, sprintf('klust_maxClust%d.mat', K)));   %'klust', 'dendrOrder', 'Tparent'); % these are sorted as in dedndrogram

cmap = cmapwinds(1:8,:);

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


colors = [205, 149 43; ...
    148, 111, 176; ...
    237, 28, 48; ...
    55, 160, 72]./255;


%% traces with standard error not viable, nor great. SKIP
[hfig, hax] = figureTracesI_PP_pixels( K, 120, 60 );
hax = flipud(hax);
for Ki = 1:K
    % parent cluster
    axes(hax(Ki,1)), hold on
    allR = R(Tparent==Ki,:);
    avg = mean(allR,1);
    avg = reshape(avg(:), [],4);
    for i = 1 : 4
        plot(ts_dec, avg(:,i), 'Color', colors(i,:)); hold on
    end
    plot([ts_dec(1), ts_dec(end)], [0 0], '-k')
    ylabel(sprintf('%d\n%d',Ki, sum(Tparent==Ki)))
    hax(Ki,1).YLabel.FontSize = 8.5;
    hax(Ki,1).XTick = [0,4];
end
linkaxes(hax, 'y')
hax(1,1).XTickLabel=[0,4];

export_fig(sprintf('clust%2d_responses_baselines.pdf',K))
% close

%% for each cluster, plot  zs(mean sustained phase) for all 4 stimuli: overlay all single pixels
load('Rzs_winds_downsampled_13FLIPPED.mat')
Rzs_sc = reshape(permute(Rzs,[2,3,1]),[], 4, size(Rzs,1) );
Rzs_sc = squeeze(mean(Rzs_sc));
for Ki = 1:K
    figure, hold on
    hax(Ki) = gca;
    allRzs_sc = Rzs_sc(:,Tparent==Ki);
    plot(allRzs_sc, '-k')
    xlim([0.5,4.5])
    title(Ki)
    hax(Ki).XTick = 1:4;
end
linkaxes(hax,'y')
for Ki = 1:K
    axes(hax(Ki))
    export_fig(sprintf('zs_sust_means_k%02d.eps',Ki))
end

%% boxplot
for Ki = 1:2 %that were split
    figure, hold on
    allRzs_sc = Rzs_sc(:,Tparent==Ki);
    boxplot(allRzs_sc', 'PlotStyle','compact')
    hax(Ki) = gca;
    title(Ki)
    hax(Ki).YLim = [-2,2];
    box off
    export_fig(sprintf('boxplot_sust_means_k%02d.eps',Ki))
end

%% boxplot all
[hfig, hax] = figureTracesI_PP_pixels( 1, repmat(80,1,K), 200 );
for Ki = 1:K
    axes(hax(Ki))
    hold on
    allRzs_sc = Rzs_sc(:,Tparent==Ki);
    boxplot(allRzs_sc', 'Colors', colors, 'Widths',0.6)
    plot(1:4, mean(allRzs_sc'), 'xk')
    hax(Ki) = gca;
    title(Ki)
    hax(Ki).YLim = [-2,2];
    hax(Ki).XGrid = 'off';
    hax(Ki).XAxis.Visible= 'off';
    hax(Ki).YTick = -2:1:2;
    if Ki>1
        hax(Ki).YAxis.Visible= 'off';
    end
    box off
end
export_fig(sprintf('boxplot_sust_means_%dALL.eps',K))


