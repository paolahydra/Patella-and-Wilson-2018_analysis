%% load
load_panSpec_downsampled;
load('/Users/galileo/Dropbox (HMS)/Data/chunks_stimfinal.mat', 'chunks')
% load('/Users/galileo/Dropbox (HMS)/Data_Raw_Metadata_BackedUp/BackedUp_Metadata/fly171_run01_metadata_2016-12-20_163728.mat', 'stimuli')
load('/Users/galileo/Dropbox (HMS)/Data/cmap19_spectral.mat', 'cmap10')

iRun = 14;
aZ = matfile(datalist_allJON{iRun}); %174_run01
T = aZ.T;
fastStimMeta = aZ.fastStimMeta;
clear aZ
for z = 1:2% NaZ
    aZ{z} = matfile(datalist_allJON{z});
end

load(fullfile(Folder2Save, 'R_matrix_Downsampled_smallWindow.mat'), 'pixel2keep', 'iKeep', 'pxKeep', 'oldpixel2keep')
    
clusterFolder = '/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix';
treeFileName = 'output_ward_movAVG_ZSCpix_linkage.mat';
cd(clusterFolder)
load('klust_maxClust10.mat')
nK = length(klust);
% load('alignedMatrix_Images_maxClust15.mat')
% cmap10 = Col(2:end,:);
% assert(size(R,1) == sum(pxKeep))


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

%% reference final clustering
cmap10H = flipud(cmap10);
figure; imagesc(1:size(cmap10H,1)); colormap(cmap10H); title('resorted cmap10 - horizontal')

K1 = 10;
T1 = load(fullfile(clusterFolder, sprintf('klust_maxClust%d.mat', K1)));
load(sprintf('singleKmaps_maxClust%d.mat', K1), 'singleKmaps')
T1.singleKmaps = singleKmaps;
clear singleKmaps
% sort these by color/finalAXIS
cOrderHoriz = fliplr(clusterOrder); %black left, A right
TR.Tparent = T1.Tparent;

for Ki = 1:K1
    K2use = cOrderHoriz(Ki);
    TR.Tparent(T1.Tparent==K2use) = Ki;
    TR.klust(Ki) = T1.klust(K2use);
end
TR.singleKmaps = T1.singleKmaps(:,:,:,cOrderHoriz);
% for Ki = 1:10, figure; imagesc(TR.singleKmaps(:,:,13,Ki)), axis image, title(Ki), end

% OK!!
clear T1

%% code for sorting and match clusters across different instantiations
% % specific cluster instantiation: I had previously run it up to 50 clusters, and also made singleKmaps in:
% edit JON_consistency_f_of_k.m
% 
% % sorting up from:
% edit batchclustering_on_rawData_PANSPEC3.m
folder2figure1 = '/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix/subclusters_SX1';
mkdir(folder2figure1)
sortFliesGenoType = 1;
load('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/anatomyImages/JON_cropMask.mat')


for k = 12:25
    T2 = load(fullfile(clusterFolder, sprintf('klust_maxClust%d.mat', k)));
    load(sprintf('singleKmaps_maxClust%d.mat', k), 'singleKmaps')
    T2.singleKmaps = singleKmaps;
    clear singleKmaps
    
    [table_T12] = crosstab(TR.Tparent, T2.Tparent);
    [TR_i, ~] = find(table_T12);
    [TR_i2R, T2_i2R] = sort(TR_i);
    
    % % check usage : YES
    % for i = 1:k
    %     cv (i) = length(T2.klust(T2_i2R(i)).k);
    % end
      
    %% code for bar plot
    % edit figure1n_maker
    AllClusterMaps_crop = T2.singleKmaps(   find(sum(cropMask,2),1) : find(sum(cropMask,2),1,'last'),...
        find(sum(cropMask,1),1) : find(sum(cropMask,1),1,'last'),...
        :,:);
    clust.ROIs = zeros(size(AllClusterMaps_crop,1)*size(AllClusterMaps_crop,2), k, size(AllClusterMaps_crop,3));           %2D matrix: rows is your tiff's (x*y), columns is # of z planes considered
    for z = 1:size(AllClusterMaps_crop,3)     %loop through z planes
        for ki = 1:k        %go through clusters
            clust.ROIs(:,ki,z) = reshape(AllClusterMaps_crop(:,:,z,ki), [], 1);
        end
    end
    flyClusterCountDivZs30_maxZout = barcluster_panSpec_maxZout(T2.Tparent, datalist, flyfoldersUnique, clust.ROIs, sortFliesGenoType);
    
    [hfig, hax] = figureTracesI_PP_pixels( k, [200], 30 ); hold on
    for i = 1:k
        axes(hax(i)), hold on
        K2use_T2 = T2_i2R(i);
        K2use_TR = TR_i2R(i);
        bar(1:NaF, flyClusterCountDivZs30_maxZout(:,K2use_T2), 'FaceColor', cmap10H(K2use_TR,:));
        axis off
    end
    
    export_fig(fullfile(folder2figure1, sprintf('clusterBars2K10_AlMaps_CountOverMaxZ_cropMap_%02d.eps',k)))
end


for k = 40
    T2 = load(fullfile(clusterFolder, sprintf('klust_maxClust%d.mat', k)));
    load(sprintf('singleKmaps_maxClust%d.mat', k), 'singleKmaps')
    T2.singleKmaps = singleKmaps;
    clear singleKmaps
    
    [table_T12] = crosstab(TR.Tparent, T2.Tparent);
    [TR_i, ~] = find(table_T12);
    [TR_i2R, T2_i2R] = sort(TR_i);
    
    % % check usage : YES
    % for i = 1:k
    %     cv (i) = length(T2.klust(T2_i2R(i)).k);
    % end
      
    %% code for bar plot
    % edit figure1n_maker
    AllClusterMaps_crop = T2.singleKmaps(   find(sum(cropMask,2),1) : find(sum(cropMask,2),1,'last'),...
        find(sum(cropMask,1),1) : find(sum(cropMask,1),1,'last'),...
        :,:);
    clust.ROIs = zeros(size(AllClusterMaps_crop,1)*size(AllClusterMaps_crop,2), k, size(AllClusterMaps_crop,3));           %2D matrix: rows is your tiff's (x*y), columns is # of z planes considered
    for z = 1:size(AllClusterMaps_crop,3)     %loop through z planes
        for ki = 1:k        %go through clusters
            clust.ROIs(:,ki,z) = reshape(AllClusterMaps_crop(:,:,z,ki), [], 1);
        end
    end
    flyClusterCountDivZs30_maxZout = barcluster_panSpec_maxZout(T2.Tparent, datalist, flyfoldersUnique, clust.ROIs, sortFliesGenoType);
    
    [hfig, hax] = figureTracesI_PP_pixels( 20, [200], 30 ); hold on
    for i = 1:20
        axes(hax(i)), hold on
        K2use_T2 = T2_i2R(i);
        K2use_TR = TR_i2R(i);
        bar(1:NaF, flyClusterCountDivZs30_maxZout(:,K2use_T2), 'FaceColor', cmap10H(K2use_TR,:));
        axis off
    end
    export_fig(fullfile(folder2figure1, sprintf('clusterBars2K10_AlMaps_CountOverMaxZ_cropMap_%02d_a.eps',k)))
    
    [hfig, hax] = figureTracesI_PP_pixels( 20, [200], 30 ); hold on
    for i = 21:40
        axes(hax(i-20)), hold on
        K2use_T2 = T2_i2R(i);
        K2use_TR = TR_i2R(i);
        bar(1:NaF, flyClusterCountDivZs30_maxZout(:,K2use_T2), 'FaceColor', cmap10H(K2use_TR,:));
        axis off
    end
    export_fig(fullfile(folder2figure1, sprintf('clusterBars2K10_AlMaps_CountOverMaxZ_cropMap_%02d_b.eps',k)))
end


