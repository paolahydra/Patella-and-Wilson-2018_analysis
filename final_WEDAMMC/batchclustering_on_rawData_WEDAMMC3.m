clusterFolder = '/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/Downsampled/rawDataLinkage/ward_ZSCpix';
treeFileName = 'output_ward_Rzsp_linkage_wedammc.mat';
load(fullfile(clusterFolder, treeFileName), 'Zc');
sizeMap = [60, 86];
map = brewermap(256, 'RdBu');
map(1:128,:) = [];

for K = 15:25
    disp(K)
%     clearvars -except R Rzs K pixel2keep iKeep pxKeep total_within Zc clusterFolder sizeMap datalist flyfoldersUnique flyNumUnique
    NaF = length(flyfoldersUnique);
    maxClust = K;
    T1 = load(fullfile(clusterFolder, sprintf('klust_maxClust%d.mat', maxClust)));
    
    T2 = load(fullfile(clusterFolder, sprintf('klust_maxClust%d.mat', maxClust+1)));
    
    
    %%
    [table_T12] = crosstab(T1.Tparent, T2.Tparent);
    figure; imagesc(table_T12); colormap(map)
    axis square
    ax = gca;
    ax.YTick = 1:maxClust;
    ax.XTick = 1:maxClust+1;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    export_fig(sprintf('table_pair_%d_%d.png',maxClust, maxClust+1))  
    
end
