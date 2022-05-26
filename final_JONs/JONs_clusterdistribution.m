% load everything from figure2n_mmaker

%% fix dataTable with z alignment information - UNATANTUM
% % manually split into 4 levels: very ventral, ventral, central, dorsal
% dataTable.zLevel
% save('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/datalist_allJON.mat', 'dataTable', '-append')

%% re-make mask to cropped area of aligned matrix as in fly-specific bar plots - UNATANTUM
% figure;
% [~, rC]  = imcrop(mean(singleKmaps_px2keep, 3));
% rC = round(rC);
% cropMask = false(size(mean(singleKmaps_px2keep, 3)));
% cropMask(rC(2):rC(2)+rC(4), rC(1):rC(1)+rC(3)) = 1;
% save('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/anatomyImages/JON_cropMask.mat', 'cropMask', 'rC')

%%
for z = 1:length(datalist)
    [runfolders{z}, ~] = fileparts(datalist{z});
end
clust = remapIndexedPointsToRun_RoisMaps(Tparent, pixel2keep, sizeMap);

% make maps (divided old/new because different dimensions)
klust = makeplot_cMaps_superKlusts_singleIteration(clust, klust, dendrOrder);

regGeneralData = matfile('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/anatomyImages/anatomy_alignment_metadata.mat');

%% one cluster at the time, for each run count how many aligned pixels fit within the aligned AMMC mask, and those that don't.
% ALSO< STORE THOSE SINGLE_CLUSTER ALIGNED MAPS
% simplified to just follow datalist/dataTable order by run

load('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/anatomyImages/JON_cropMask.mat')

alClustCount = zeros(length(datalist), nK, length(unique(dataTable.zLevel))); 

for z = 1: size(dataTable,1)
    disp(z)
    flyNum = dataTable.fly(z);
    flyfolder = flyfolders{z};
    
    % load fly's stack2fly150 transformation
    stack2fly150_info = matfile(fullfile(flyfolder, ['stack2fly150sStack_' basenames{z}(1:6) '.mat']));
    if ~isfield(stack2fly150_info, 'Roriginal') && ~strcmp(basenames{z}(1:6), 'fly150')
        stack2fly150_info.Properties.Writable = true;
        stack2fly150_info.Roriginal = imref2d(size(stack2fly150_info.recovered));
    end
    
    
    % load specific run's transformation
    regData = matfile(fullfile(runfolders{z}, ['stackRegistration_' basenames{z} '.mat']));
    if ~isfield(regData, 'Roriginal')
        regData.Properties.Writable = true;
        regData.Roriginal = imref2d(size(regData.target));
    end
    
    al2fly150_mapAll = [];
    for Ki = 1:length(unique(Tparent))  %stick to Ki while building. Use K2use when retrieving 
        
        if Ki == 1
            mapClean = klust(Ki).mapsZold(:,:,z);
            mapClean(mapClean>=1) = 1; %binary - all px2keep
            mapClean = cat(1, mapClean, zeros(1,size(mapClean, 2))); %61 x 86
            mapClean = cat(2, mapClean, zeros(size(mapClean, 1), 4)); %61 x 90
            mapClean = imresize(mapClean, [86 128]); %should become 86 x 128. The quality sucks, but so be it
            %re-crop
            mapClean(end,:) = [];
            mapClean(:,end-5:end) = [];
            mapClean(mapClean<0.3) = 0;
            mapClean(mapClean>=0.1) = 1;                        %note 1 and not Ki!
            
            mapClean = imrotate(mapClean, regGeneralData.totAng(1,z)); %preliminary transformation
            recovered_mapClean = imwarp(mapClean,regData.tform,'OutputView',regData.Roriginal,'FillValues',0); %indexed image
            if strcmp(basenames{z}(1:6), 'fly150')
                al2fly150_mapClean = recovered_mapClean;
            else
                al2fly150_mapClean = imwarp(recovered_mapClean,stack2fly150_info.tform,'OutputView',stack2fly150_info.Roriginal,'FillValues',0);
            end
            al2fly150_mapClean(al2fly150_mapClean<=0.3) = 0;
            al2fly150_mapClean(al2fly150_mapClean>=0.1) = 1; %dendrogram is colored in order from 1st to nK-esim  in colormap. Clusters are colored accordingy.
            singleKmaps_px2keep(:,:,z) = al2fly150_mapClean;
        end
        
        
        mapClean = klust(Ki).mapsZold(:,:,z);    
        mapClean(mapClean == nK+1) = 0;
        mapClean(mapClean>=1) = 1; %binary so far
        
%         figure; imshow(mapClean, [])
        
        %convert from downsampled to fullsize
        mapClean = cat(1, mapClean, zeros(1,size(mapClean, 2))); %61 x 86
        mapClean = cat(2, mapClean, zeros(size(mapClean, 1), 4)); %61 x 90
        mapClean = imresize(mapClean, [86 128]); %should become 86 x 128. The quality sucks, but so be it
        %re-crop
        mapClean(end,:) = [];
        mapClean(:,end-5:end) = [];
        mapClean(mapClean<0.3) = 0;
        mapClean(mapClean>=0.1) = 1;                        %note 1 and not Ki!
        
        mapClean = imrotate(mapClean, regGeneralData.totAng(1,z)); %preliminary transformation
        recovered_mapClean = imwarp(mapClean,regData.tform,'OutputView',regData.Roriginal,'FillValues',0); %indexed image
        if strcmp(basenames{z}(1:6), 'fly150')
            al2fly150_mapClean = recovered_mapClean;
        else
            al2fly150_mapClean = imwarp(recovered_mapClean,stack2fly150_info.tform,'OutputView',stack2fly150_info.Roriginal,'FillValues',0);
        end
            
        al2fly150_mapClean(al2fly150_mapClean<=0.3) = 0;
        al2fly150_mapClean(al2fly150_mapClean>=0.1) = 1; %dendrogram is colored in order from 1st to nK-esim  in colormap. Clusters are colored accordingy.
        
        singleKmaps(:,:,z,Ki) = al2fly150_mapClean;
        
        % split up different z levels: split more - you can lump later 
        % let's only consider a specific registered window (as already done
        % before) -- where?: bar plots distributions
        alClustCount(z, Ki, dataTable.zLevel(z)) = sum(sum( al2fly150_mapClean(cropMask) ));
        
    end
end

save('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix/10clusters/clusterHeatMaps/alignedData.mat', ...
    'alClustCount', 'singleKmaps', 'singleKmaps_px2keep')

%% make new distribution plots - NOW START USING K2USE!!!!
% only count in a cropped region to exclude outliers/ nerve recordings
for z = 1: size(dataTable,1)
    for Ki = 1:nK
        al2fly150_mapClean = singleKmaps(:,:,z,Ki);
        alClustCount(z, Ki, dataTable.zLevel(z)) = sum(sum( al2fly150_mapClean(cropMask) ));
    end
end


alClustCount_splits = squeeze(sum(alClustCount, 1)); 
alClustCount_splits = alClustCount_splits(clusterOrder,:); %resorted as in axis from bottom to up


% alClustCount_splits = alClustCount_splits(:,[1 3 2 4]);  % IPSI AMMC, ipsi wed, CONTRA AMMC, contra wed
% alClustCount_splits = bsxfun(@rdivide, alClustCount_splits, reshape(repmat(max(reshape(alClustCount_splits, [], 2)), 2,1), 1,4)  );

figure; 
b = bar(alClustCount_splits);
h = gca;
h.XTick = 1:nK;
h.Box = 'off';
% for i = 1:4
%     b(i).BarWidth = 1;
% end
legend('very ventral', 'ventral', 'central', 'dorsal')
legend boxoff
ylabel('registered pixel count over cropped area')
xlabel('cluster number as in axis, bottom-up')
export_fig(fullfile(folder2figure2, 'clusterDistrib_pixelCount.eps'))

% USE --- normalize by the maximum count within each level
alClustCount_splits_n = bsxfun(@rdivide, alClustCount_splits, max(alClustCount_splits,[],1) );

figure; 
b = bar(alClustCount_splits_n);
h = gca;
h.XTick = 1:nK;
h.Box = 'off';
% for i = 1:4
%     b(i).BarWidth = 1;
% end
legend('very ventral', 'ventral', 'central', 'dorsal')
legend boxoff
ylabel(sprintf('registered pixel count over cropped area,\nnormalized to maximum within each level'))
xlabel('cluster number as in axis, bottom-up')
export_fig(fullfile(folder2figure2, 'clusterDistrib_ClusterDominanceWithinEachLevel.eps'))


%  normalize both within levels and within a cluster
alClustCount_splits_n = bsxfun(@rdivide, alClustCount_splits, max(alClustCount_splits,[],1) );
alClustCount_splits_nC = bsxfun(@rdivide, alClustCount_splits_n, max(alClustCount_splits_n,[],2) );
alClustCount_splits_nC = fliplr(alClustCount_splits_nC); %dorsal to ventral

figure; 
b = bar(alClustCount_splits_nC);
h = gca;
h.XTick = 1:nK;
h.Box = 'off';
% for i = 1:4
%     b(i).BarWidth = 1;
% end
% legend boxoff
export_fig(fullfile(folder2figure2, 'clusterDistrib_ClustersDorsal2Ventral.eps'))



% 
% %
% stackData = squeeze(sum(alClustCount, 1)); % AMMC IPSI CONTRA, wed ipsi contra
% stackData_splits = stackData(:,[1 3 2 4]);  % IPSI AMMC, ipsi wed, CONTRA AMMC, contra wed %gola: normalize within all-ipsi and all-contra
% stackData_splits = max(reshape(stackData_splits, [], 2));
% divis_stackData = repmat( stackData_splits, 1, 2);  %ipsi contra ipsi contra
% stackData = bsxfun(@rdivide, stackData, divis_stackData );
% stackData = permute(reshape(stackData,[], 2, 2), [1,3,2] );
% groupLabels = num2cell(1:19);
% hb = plotBarStackGroups(stackData, groupLabels);
% hb(1,1).FaceColor = [1 0.4 0.4];   %AMMC ipsi
% hb(2,1).FaceColor = [0.4 0.4 1];   %wed ipsi
% hb(1,2).FaceColor = [0.6  0  0];   %AMMC contra
% hb(2,2).FaceColor = [0  0  0.6];   %wed contra
% export_fig(fullfile(folder2figure, 'clusterDistrib_AMMC_WED_ipsicontraStacked.eps'))


%% NEW-2 DISTRIB PLOTS
load('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix/10clusters/clusterHeatMaps/alignedData.mat')
load('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/anatomyImages/JON_cropMask.mat')
alClustCount = zeros(length(datalist), nK, length(unique(dataTable.zLevel)));
for z = 1: size(dataTable,1)
    for Ki = 1:nK
        al2fly150_mapClean = singleKmaps(:,:,z,Ki);
        alClustCount(z, Ki, dataTable.zLevel(z)) = sum(sum( al2fly150_mapClean(cropMask) ));
    end
end
h = histogram(dataTable.zLevel);
alClustCount_splits = squeeze(sum(alClustCount, 1));
alClustCount_splits = bsxfun(@rdivide, alClustCount_splits, h.Values); %not a pure count, but an average over the planes that actually contributed
alClustCount_splits = alClustCount_splits(clusterOrder,:); %resorted as in axis from bottom to up
close;

% un-normalized mean representation (averaged within Z level)
figure; 
b = bar(fliplr(alClustCount_splits));
h = gca;
h.XTick = 1:nK;
h.Box = 'off';
% for i = 1:4
%     b(i).BarWidth = 1;
% end
legend(fliplr({'very ventral', 'ventral', 'central', 'dorsal'}))
legend boxoff
ylabel('registered pixel mean over cropped area')
xlabel('cluster number as in axis, bottom-up')
export_fig(fullfile(folder2figure2, 'clusterDistrib_pixelMean_unNorm.eps'))





%% now make and save heatmaps -K2USE!!!!
load('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix/10clusters/clusterHeatMaps/alignedData.mat');

fall = '/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix/10clusters/clusterHeatMaps';
cd(fall)

for Z = 1:4
    levMaps_pixelMax = sum(singleKmaps_px2keep(:,:,dataTable.zLevel==Z),3);
    for Ki = 1:nK
        K2use = clusterOrder(Ki);
        map = sum(singleKmaps(:,:,dataTable.zLevel==Z,K2use), 3) ./ levMaps_pixelMax;
        map(~isfinite(map)) = 0;
        map = map.^(1/4);
        map = 1 - mat2gray(map);
        savename = sprintf('norm_heatmap_cluster%02d_lev%d.tif',Ki,Z);
        imwrite(map, savename)
    end
end


alMaps_pixelMax = sum(singleKmaps_px2keep,3);
for Ki = 1:nK  
    K2use = clusterOrder(Ki);
    map = sum(singleKmaps(:,:,:,K2use), 3) ./ alMaps_pixelMax;
    map(~isfinite(map)) = 0;
    map = map.^(1/4);
    map = 1 - mat2gray(map);
    savename = sprintf('norm_heatmap_cluster%02d_lev9.tif',Ki);
    imwrite(map, savename)
end


%% save a pmaxproj_px2keepMap for level 3
Z = 3;
levMaps_pixelMax = sum(singleKmaps_px2keep(:,:,dataTable.zLevel==Z),3);
levMaps_pixelMax(levMaps_pixelMax>=1) = 1;
levMaps_pixelMax = logical(levMaps_pixelMax);
levMaps_pixelMax = ~levMaps_pixelMax;
savename = sprintf('px2keep_maxProj_lev%d.tif',Z);
imwrite(levMaps_pixelMax, savename)

%% tentative colormap
x = 0:0.001:1;
x = x.^(2/5);
x = 1-x;
x = repmat(x, 100, 1);
figure; imshow(x)
export_fig(fullfile(folder2figure, 'tentativeColormap_normHeatMaps.tif'))




