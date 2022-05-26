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

% do I want to separate between AMMC/WED? no
% do I care about differentially normalizing ipsi and contra? maybe, but as
% a first step no (for prob map)
% do I want to make a probability map, or a proportion/prevalence map? the
% latter.
% what do I need for that? singlepixel_counts, custCount (registered)

%%
alClustCount = zeros(length(datalist), nK, 3); % third dimension has two levels: 1 takes count of AMMC pixels, second of ipsilateral notAMMC pixels, third of contralateral notAMMC pixels
for z = 11: size(dataTable,1)
    disp(z)
    orZ = find( orDT.fly == dataTable.fly(z) & orDT.run == dataTable.run(z) ); %used to reference regGeneralData
    flyNum = dataTable.fly(z);
    i_allSt = find(allStacks_flyNums == flyNum); 
    flyfolder = flyfolders{z};
    % load specific run's transformation
    regData = matfile(fullfile(runfolders{z}, ['stackRegistration_' basenames{z} '.mat']));
    if ~isfield(regData, 'Roriginal')
        regData.Properties.Writable = true;
        regData.Roriginal = imref2d(size(regData.target));
    end
    
    al2fly119_mapAll = [];
    for Ki = 1:length(unique(Tparent))
        
        if Ki == 1
            if flyNum == 175
                mapClean = klust(Ki).mapsZold(:,:,z);
                mapClean(mapClean>=1) = 1;  %binary - all px2keep
                mapClean = cat(1, mapClean, zeros(1,size(mapClean, 2)));    % 70 x 86
                mapClean = cat(2, mapClean, zeros(size(mapClean, 1), 4));   % 70 x 90
                mapClean = imresize(mapClean, [100 128]);                   % 100 x 128
                %re-crop
                mapClean(end,:) = [];
                mapClean(:,end-5:end) = [];                                 % 99 x 122
                mapClean(mapClean<0.3) = 0;
                mapClean(mapClean>=0.1) = 1;                        %note 1 and not K2use!
            else
                z_86 = z - NaZ_100;         % now this is relative to the R_86 subset
                mapClean = klust(Ki).mapsZnew(:,:,z_86);
                mapClean(mapClean>=1) = 1;  %binary - all px2keep
                mapClean = cat(1, mapClean, zeros(1,size(mapClean, 2))); %61 x 86
                mapClean = cat(2, mapClean, zeros(size(mapClean, 1), 4)); %61 x 90
                mapClean = imresize(mapClean, [86 128]); %should become 86 x 128. The quality sucks, but so be it
                %re-crop
                mapClean(end,:) = [];
                mapClean(:,end-5:end) = [];
                mapClean(mapClean<0.3) = 0;
                mapClean(mapClean>=0.1) = 1;                        %note 1 and not Ki!
            end
            mapClean = imrotate(mapClean, regGeneralData.totAng(1,orZ)); %preliminary transformation
            recovered_mapClean = imwarp(mapClean,regData.tform,'OutputView',regData.Roriginal,'FillValues',0); %indexed image
            al2fly119_mapClean = imwarp(recovered_mapClean,allStacks(i_allSt).tform,'OutputView',allStacks(i_allSt).Roriginal,'FillValues',0);
            al2fly119_mapClean(al2fly119_mapClean<=0.3) = 0;
            al2fly119_mapClean(al2fly119_mapClean>=0.1) = 1; %dendrogram is colored in order from 1st to nK-esim  in colormap. Clusters are colored accordingy.
            singleKmaps_px2keep(:,:,z) = al2fly119_mapClean;
        end
        
        if flyNum == 175
            mapClean = klust(Ki).mapsZold(:,:,z);
            mapClean(mapClean == nK+1) = 0;
            mapClean(mapClean>=1) = 1; %binary so far
            %convert from downsampled to fullsize
            mapClean = cat(1, mapClean, zeros(1,size(mapClean, 2)));    % 70 x 86
            mapClean = cat(2, mapClean, zeros(size(mapClean, 1), 4));   % 70 x 90
            mapClean = imresize(mapClean, [100 128]);                   % 100 x 128
            %re-crop
            mapClean(end,:) = [];
            mapClean(:,end-5:end) = [];                                 % 99 x 122
            mapClean(mapClean<0.3) = 0;
            mapClean(mapClean>=0.1) = 1;                        %note 1 and not K2use!
        else
            z_86 = z - NaZ_100;         % now this is relative to the R_86 subset
            mapClean = klust(Ki).mapsZnew(:,:,z_86);
            mapClean(mapClean == nK+1) = 0;
            mapClean(mapClean>=1) = 1; %binary so far
            %convert from downsampled to fullsize
            mapClean = cat(1, mapClean, zeros(1,size(mapClean, 2))); %61 x 86
            mapClean = cat(2, mapClean, zeros(size(mapClean, 1), 4)); %61 x 90
            mapClean = imresize(mapClean, [86 128]); %should become 86 x 128. The quality sucks, but so be it
            %re-crop
            mapClean(end,:) = [];
            mapClean(:,end-5:end) = [];
            mapClean(mapClean<0.3) = 0;
            mapClean(mapClean>=0.1) = 1;                        %note 1 and not Ki!
        end
        mapClean = imrotate(mapClean, regGeneralData.totAng(1,orZ)); %preliminary transformation
        recovered_mapClean = imwarp(mapClean,regData.tform,'OutputView',regData.Roriginal,'FillValues',0); %indexed image
        al2fly119_mapClean = imwarp(recovered_mapClean,allStacks(i_allSt).tform,'OutputView',allStacks(i_allSt).Roriginal,'FillValues',0);
        al2fly119_mapClean(al2fly119_mapClean<=0.3) = 0;
        al2fly119_mapClean(al2fly119_mapClean>=0.1) = 1; %dendrogram is colored in order from 1st to nK-esim  in colormap. Clusters are colored accordingy.
        
        singleKmaps(:,:,z,Ki) = al2fly119_mapClean;
%         % apply corresponding mask and count
%         thisMZ = dataTable.matrixsZ(z) + 11;
%         
%         if dataTable.ipsilateral(z)
%             alClustCount(z, Ki, 1) = sum(sum( al2fly119_mapClean(alignedAMMCmasks(:,:,thisMZ)) )); %AMMC count ipsi (1)
%             alClustCount(z, Ki, 3) = sum(sum( al2fly119_mapClean(~alignedAMMCmasks(:,:,thisMZ)) )); %ipsi notAMMC count
%         else
%             alClustCount(z, Ki, 2) = sum(sum( al2fly119_mapClean(alignedAMMCmasks(:,:,thisMZ)) )); %AMMC count contra (2)
%             alClustCount(z, Ki, 4) = sum(sum( al2fly119_mapClean(~alignedAMMCmasks(:,:,thisMZ)) )); %contra notAMMC count
%         end
    end
end
save(sprintf('/Users/galileo/Dropbox (HMS)/Data/WINDS_100/Downsampled/SubsetFlies/rawDataLinkage/%dclusters/clusterHeatMaps/alignedData.mat',K), ...
    'alClustCount', 'singleKmaps', 'singleKmaps_px2keep')

%%  heatmaps 
% load('/Users/galileo/Dropbox (HMS)/Data/WINDS_100/Downsampled/SubsetFlies/rawDataLinkage/7clusters/clusterHeatMaps/alignedData.mat');



fall = sprintf('/Users/galileo/Dropbox (HMS)/Data/WINDS_100/Downsampled/SubsetFlies/rawDataLinkage/%dclusters/clusterHeatMaps/conditionalP',K);
mkdir(fall)
cd(fall)
prob_px2keep = mean(singleKmaps_px2keep,3);
for Ki = 1:nK
    map = mean(singleKmaps(:,:,:,Ki), 3) .* prob_px2keep;
    map = map.^(1/4);
    map = 1 - mat2gray(map);
    
    savename = sprintf('condProb_heatmap_cluster%02d.tif',Ki);
    imwrite(map, savename)
    savename = sprintf('norm_heatmap_cluster%02d.tif',Ki);
    imwrite(map, savename)
end


fall = sprintf('/Users/galileo/Dropbox (HMS)/Data/WINDS_100/Downsampled/SubsetFlies/rawDataLinkage/%dclusters/clusterHeatMaps/regularProbabilityHeatmap',K);
mkdir(fall)
cd(fall)
for Ki = 1:nK
    map = mean(singleKmaps(:,:,:,Ki), 3);
    map = map.^(1/4);
    map = 1 - mat2gray(map);
    
    savename = sprintf('regProb_heatmap_cluster%02d.tif',Ki);
    imwrite(map, savename)
    savename = sprintf('norm_heatmap_cluster%02d.tif',Ki);
    imwrite(map, savename)
end


fall = sprintf('/Users/galileo/Dropbox (HMS)/Data/WINDS_100/Downsampled/SubsetFlies/rawDataLinkage/%dclusters/clusterHeatMaps/proportion',K);
mkdir(fall)
cd(fall)
alMaps_pixelMax = sum(singleKmaps_px2keep,3);
for Ki = 1:nK
    map = sum(singleKmaps(:,:,:,Ki), 3) ./ alMaps_pixelMax;
    map(~isfinite(map)) = 0;
    map = map.^(1/4);
    map = 1 - mat2gray(map);
    
    savename = sprintf('norm_heatmap_cluster%02d.tif',Ki);
    imwrite(map, savename)
end




