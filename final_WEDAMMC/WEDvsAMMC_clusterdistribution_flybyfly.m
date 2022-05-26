load_wedAmmc_downsampled;
load('/Users/galileo/Dropbox (HMS)/Data/chunks_stimfinal.mat', 'chunks')
% load('/Users/galileo/Dropbox (HMS)/Data/cmap19_spectral.mat', 'cmap')
% save('klust_maxClust19.mat', 'cmap', '-append')
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

%% separate data between AMMC and WED
% Any way wuold be approximative, but the best way I (already) have, is to
% use the AMMC masks I created.

%% 1. find where I saved the most recent masks, and load them
% -- SAVED INTO the DATALIST file as alignedAMMCmasks - this is the final version
% in z it follows the same numbering as matrixsZ+11.

%% align to specific runs and split labeled pixels into ipsi-wed, contra-wed, ipsiAMMC
%% align to fly 119
for z = 1:length(datalist)
    [runfolders{z}, ~] = fileparts(datalist{z});
end
clust = remapIndexedPointsToRun_RoisMaps(Tparent, pixel2keep, sizeMap);

% make maps (divided old/new because different dimensions)
klust = makeplot_cMaps_superKlusts_singleIteration(clust, klust, dendrOrder);
% folder = 'alignedMaps';
% mkdir(folder)
% cd(folder)

Col = cat(1,[1,1,1],cmap);

regGeneralData = matfile('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/anatomyImages/anatomy_alignment_metadata.mat');
load('/Users/galileo/Dropbox (HMS)/Data/TDTstacks/anatomy_alignment_metadata.mat') %allStacks
allStacks_flyNums = cat(1,allStacks.flyNum);

% %% first, rerun cluster maps with new sortedcluster and no brown
% cd('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/Downsampled/rawDataLinkage/ward_ZSCpix/19clusters/alignedMaps')
% 
% for f = 1 : NaF
%     flyNum = flyNumUnique(f);
%     flyfolder = flyfoldersUnique{f};
%     ZZ = strfind(datalist, flyfolder);
%     zs_fly = zeros(size(ZZ));
%     for i = 1:length(zs_fly)
%         if ~isempty(ZZ{i})
%             zs_fly(i) = 1; % indices into datalist relative to fly f
%         end
%     end
%     zs_fly = logical(zs_fly);
%     
%     % find corresponding allStacks row
%     if flyNum == 118
%         i_allSt = find(allStacks_flyNums == 126);
%     else
%         i_allSt = find(allStacks_flyNums == flyNum);
%     end
%     
%     %% sort runs from ventral to dorsal
%     zetas = find(zs_fly);
%     sliceNums = regGeneralData.alignedSliceNumbers(1,zetas);
%     [sliceNums,b] = sort(sliceNums);
%     zetas = zetas(b);
% 
%     %% load specific runs data and transformations and compute, correcting for downsampling
%     for iz = 1:sum(zs_fly)
%         z = zetas(iz);
%         sliceNum = sliceNums(iz);
%         % load specific run's transformation
%         regData = matfile(fullfile(runfolders{z}, ['stackRegistration_' basenames{z} '.mat']));
%         if ~isfield(regData, 'Roriginal')
%             regData.Properties.Writable = true;
%             regData.Roriginal = imref2d(size(regData.target));
%         end
% 
%         al2fly119_mapAll = [];
%         for Ki = 1:length(unique(Tparent))
%             mapClean = klust(Ki).mapsZold(:,:,z);
%             mapClean(mapClean == nK+1) = 0;
%             mapClean(mapClean>=1) = 1; %binary so far
% 
%             %convert from downsampled to fullsize
%             mapClean = cat(1, mapClean, zeros(1,size(mapClean, 2))); %61 x 86
%             mapClean = cat(2, mapClean, zeros(size(mapClean, 1), 4)); %61 x 90
%             mapClean = imresize(mapClean, [86 128]); %should become 86 x 128. The quality sucks, but so be it
%             %re-crop
%             mapClean(end,:) = [];
%             mapClean(:,end-5:end) = [];
%             mapClean(mapClean<0.3) = 0;
%             mapClean(mapClean>=0.1) = 1;                        %note 1 and not Ki!
% 
%             mapClean = imrotate(mapClean, regGeneralData.totAng(1,z)); %preliminary transformation
%             recovered_mapClean = imwarp(mapClean,regData.tform,'OutputView',regData.Roriginal,'FillValues',0); %indexed image
%             al2fly119_mapClean = imwarp(recovered_mapClean,allStacks(i_allSt).tform,'OutputView',allStacks(i_allSt).Roriginal,'FillValues',0);
%             al2fly119_mapClean(al2fly119_mapClean<=0.5) = 0;
%             al2fly119_mapClean(al2fly119_mapClean>=0.1) = Ki; %dendrogram is colored in order from 1st to nK-esim  in colormap. Clusters are colored accordingy.
%             if isempty(al2fly119_mapAll)
%                 al2fly119_mapAll = al2fly119_mapClean;
%             else
%                 al2fly119_mapAll(al2fly119_mapClean == Ki) = Ki; %dendrogram is colored in order from 1st to nK-esim  in colormap. Clusters are colored accordingy.
%             end
%         end
% 
%         AllClusterMaps(:,:,z) = al2fly119_mapAll;
% 
%         AlphaMatrix = zeros(size(al2fly119_mapAll(:,:,1)));
%         AlphaMatrix(al2fly119_mapAll(:,:,1) > 0) = 1;
% 
%         al2fly119_mapAll = ind2rgb(uint8(al2fly119_mapAll), Col);
%         savename = sprintf('alignALLklustsMap_fly%02d_ordZ%02d_%s.png', f, z, basenames{z});
%         imwrite(al2fly119_mapAll, savename, 'Alpha', double(AlphaMatrix))
%     end
% end

%% one cluster at the time, for each run count how many aligned pixels fit within the aligned AMMC mask, and those that don't AND that FIT THE WED MASK (new)
% ALSO< STORE THOSE SINGLE_CLUSTER ALIGNED MAPS
% simplified to just follow datalist/dataTable order by run

alClustCount = zeros(length(datalist), nK, 3); % third dimension has two levels: 1 takes count of AMMC pixels, second of ipsilateral notAMMC pixels, third of contralateral notAMMC pixels

for z = 1: size(dataTable,1)
    disp(z)
    flyNum = dataTable.fly(z);
    flyfolder = flyfolders{z};
    % find corresponding allStacks row
    if flyNum == 118
        i_allSt = find(allStacks_flyNums == 126);
    else
        i_allSt = find(allStacks_flyNums == flyNum);
    end
    
    % load specific run's transformation
    regData = matfile(fullfile(runfolders{z}, ['stackRegistration_' basenames{z} '.mat']));
    if ~isfield(regData, 'Roriginal')
        regData.Properties.Writable = true;
        regData.Roriginal = imref2d(size(regData.target));
    end
    
    al2fly119_mapAll = [];
    for Ki = 1:length(unique(Tparent))
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
        al2fly119_mapClean = imwarp(recovered_mapClean,allStacks(i_allSt).tform,'OutputView',allStacks(i_allSt).Roriginal,'FillValues',0);
        al2fly119_mapClean(al2fly119_mapClean<=0.3) = 0;
        al2fly119_mapClean(al2fly119_mapClean>=0.1) = 1; %dendrogram is colored in order from 1st to nK-esim  in colormap. Clusters are colored accordingy.
        
        singleKmaps(:,:,z,Ki) = al2fly119_mapClean;
        % apply corresponding mask and count
        thisMZ = dataTable.matrixsZ(z) + 11;
        
        if dataTable.ipsilateral(z)
            alClustCount(z, Ki, 1) = sum(sum( al2fly119_mapClean(alignedAMMCmasks(:,:,thisMZ)) )); %AMMC count ipsi (1)
            alClustCount(z, Ki, 3) = sum(sum( al2fly119_mapClean(alignedWEDmask & ~alignedAMMCmasks(:,:,thisMZ)) )); %ipsi notAMMC count
        else
            alClustCount(z, Ki, 2) = sum(sum( al2fly119_mapClean(alignedAMMCmasks(:,:,thisMZ)) )); %AMMC count contra (2)
            alClustCount(z, Ki, 4) = sum(sum( al2fly119_mapClean(alignedWEDmask & ~alignedAMMCmasks(:,:,thisMZ)) )); %contra notAMMC count
        end
    end
   
end

%% actual cluster maps split into AMMC WED
for z = 1: size(dataTable,1)
    singleKmapsAMMC(:,:,z,:) = bsxfun( @and, singleKmaps(:,:,z,:),  alignedAMMCmasks(:,:,dataTable.matrixsZ(z)+11) );
    singleKmapsWED(:,:,z,:) = bsxfun( @and, singleKmaps(:,:,z,:),  (alignedWEDmask & ~ alignedAMMCmasks(:,:,dataTable.matrixsZ(z)+11)) );
end

% max project each cluster and each region
%% AMMC
maxAMMC = sum( sum(singleKmapsAMMC,4), 3);
maxAMMC(maxAMMC>1) = 1;
figure; imshow(maxAMMC,[]); title('max projection AMMC')
areaAMMC = sum(maxAMMC(:));

%% WED
maxWED = sum( sum(singleKmapsWED,4), 3);
maxWED(maxWED>1) = 1;
figure; imshow(maxWED,[]); title('max projection WED')
areaWED = sum(maxWED(:));

%% consider also each cluster max project within each region
for Ki = 1:length(unique(Tparent)) 
    % AMMC
    maxAMMC_K(:,:,Ki) = max(singleKmapsAMMC(:,:,:,Ki), [], 3);
    areaAMMC_K(Ki) = sum(sum(maxAMMC_K(:,:,Ki)));
    % WED
    maxWED_K(:,:,Ki) = max(singleKmapsWED(:,:,:,Ki), [], 3);
    areaWED_K(Ki) = sum(sum(maxWED_K(:,:,Ki)));
end

areaAMMC_K = areaAMMC_K./sum(areaAMMC_K);
areaWED_K = areaWED_K./sum(areaWED_K);

figure; bar( cat(2, areaAMMC_K', areaWED_K') )

save('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix/relativeWeightClusters_AMMC_WED.mat', 'areaAMMC_K', 'areaWED_K')
%% make new distribution plots
alClustCount_splits = squeeze(sum(alClustCount, 1));
alClustCount_splits = alClustCount_splits(:,[1 3 2 4]);  % IPSI AMMC, ipsi wed, CONTRA AMMC, contra wed
 
alClustCount_splits = bsxfun(@rdivide, alClustCount_splits, reshape(repmat(max(reshape(alClustCount_splits, [], 2)), 2,1), 1,4)  );

figure; 
b = bar(alClustCount_splits);
h = gca;
h.XTick = 1:nK;
h.Box = 'off';
b(1).FaceColor = [1 0.4 0.4];   %AMMC ipsi
b(2).FaceColor = [0.4 0.4 1];   %wed ipsi
b(3).FaceColor = [0.6  0  0];   %AMMC contra
b(4).FaceColor = [0  0  0.6];   %wed contra
for i = 1:4
    b(i).BarWidth = 1;
end
legend('ipsi AMMC', 'ipsi WED', 'contra AMMC', 'contra WED')
legend boxoff
export_fig(fullfile(folder2figure, 'clusterDistrib_AMMC_WED_ipsi_contra.eps'))


%
stackData = squeeze(sum(alClustCount, 1)); % AMMC IPSI CONTRA, wed ipsi contra
stackData_splits = stackData(:,[1 3 2 4]);  % IPSI AMMC, ipsi wed, CONTRA AMMC, contra wed %gola: normalize within all-ipsi and all-contra
stackData_splits = max(reshape(stackData_splits, [], 2));
divis_stackData = repmat( stackData_splits, 1, 2);  %ipsi contra ipsi contra
stackData = bsxfun(@rdivide, stackData, divis_stackData );
stackData = permute(reshape(stackData,[], 2, 2), [1,3,2] );
groupLabels = num2cell(1:19);
hb = plotBarStackGroups(stackData, groupLabels);
hb(1,1).FaceColor = [1 0.4 0.4];   %AMMC ipsi
hb(2,1).FaceColor = [0.4 0.4 1];   %wed ipsi
hb(1,2).FaceColor = [0.6  0  0];   %AMMC contra
hb(2,2).FaceColor = [0  0  0.6];   %wed contra
export_fig(fullfile(folder2figure, 'clusterDistrib_AMMC_WED_ipsicontraStacked.eps'))


%% make aligned px2keep maps for pixel-wise normalization
for z = 1: size(dataTable,1)
    disp(z)
    flyNum = dataTable.fly(z);
    flyfolder = flyfolders{z};
    % find corresponding allStacks row
    if flyNum == 118
        i_allSt = find(allStacks_flyNums == 126);
    else
        i_allSt = find(allStacks_flyNums == flyNum);
    end
    
    % load specific run's transformation
    regData = matfile(fullfile(runfolders{z}, ['stackRegistration_' basenames{z} '.mat']));
    if ~isfield(regData, 'Roriginal')
        regData.Properties.Writable = true;
        regData.Roriginal = imref2d(size(regData.target));
    end
    
    al2fly119_mapAll = [];
    Ki = 1;
    mapClean = klust(Ki).mapsZold(:,:,z);
    mapClean(mapClean>=1) = 1; %binary - all px2keep
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
    al2fly119_mapClean = imwarp(recovered_mapClean,allStacks(i_allSt).tform,'OutputView',allStacks(i_allSt).Roriginal,'FillValues',0);
    al2fly119_mapClean(al2fly119_mapClean<=0.3) = 0;
    al2fly119_mapClean(al2fly119_mapClean>=0.1) = 1; %dendrogram is colored in order from 1st to nK-esim  in colormap. Clusters are colored accordingy.
    
    singleKmaps_px2keep(:,:,z) = al2fly119_mapClean;
end
save('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/Downsampled/rawDataLinkage/ward_ZSCpix/19clusters/clusterHeatMaps/alignedData.mat', ...
    'alClustCount', 'singleKmaps', 'singleKmaps_px2keep')

%% now make and save heatmaps
load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/Downsampled/rawDataLinkage/ward_ZSCpix/19clusters/clusterHeatMaps/alignedData.mat');



% % ventral portion
% M = (dataTable.matrixsZ+11) <=21;
% ventralMaps_pixelMax = sum(singleKmaps_px2keep(:,:,M),3);
% for Ki = 1:nK
%     ventralMaps(:,:,Ki) = squeeze(sum(singleKmaps(:,:,M,Ki), 3) );
%     ventralMaps_norm(:,:,Ki) = ventralMaps(:,:,Ki) ./ ventralMaps_pixelMax;
% end
% 
% 
% M = (dataTable.matrixsZ+11) > 21;
% dorsalMaps_pixelMax = sum(singleKmaps_px2keep(:,:,M),3);
% for Ki = 1:nK
%     dorsalMaps(:,:,Ki) = squeeze(sum(singleKmaps(:,:,M,Ki), 3) );
%     dorsalMaps_norm(:,:,Ki) = dorsalMaps(:,:,Ki) ./ dorsalMaps_pixelMax;
% end
% 
% 
% mkdir('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/Downsampled/rawDataLinkage/ward_ZSCpix/19clusters/clusterHeatMaps')
% fventral = '/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/Downsampled/rawDataLinkage/ward_ZSCpix/19clusters/clusterHeatMaps/ventralPlane';
% mkdir(fventral)
% cd(fventral)
% for Ki = 1:nK
%     map = ventralMaps(:,:,Ki);
%     map = 1 - mat2gray(map);
%     savename = sprintf('heatmap_cluster%02d.tif',Ki);
%     imwrite(map, savename)
%     
%     map = ventralMaps_norm(:,:,Ki);
%     map(~isfinite(map)) = 0;
%     map = 1 - mat2gray(map);
%     savename = sprintf('norm_heatmap_cluster%02d.tif',Ki);
%     imwrite(map, savename)
% end
% 
% fdorsal = '/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/Downsampled/rawDataLinkage/ward_ZSCpix/19clusters/clusterHeatMaps/dorsalPlane';
% mkdir(fdorsal)
% cd(fdorsal)
% for Ki = 1:nK
%     map = dorsalMaps(:,:,Ki);
%     map = 1 - mat2gray(map);
%     savename = sprintf('heatmap_cluster%02d.tif',Ki);
%     imwrite(map, savename)
%     
%     map = dorsalMaps_norm(:,:,Ki);
%     map(~isfinite(map)) = 0;
%     map = 1 - mat2gray(map);
%     savename = sprintf('norm_heatmap_cluster%02d.tif',Ki);
%     imwrite(map, savename)
% end


fall = '/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/Downsampled/rawDataLinkage/ward_ZSCpix/19clusters/clusterHeatMaps/allPlanes';
mkdir(fall)
cd(fall)

alMaps_pixelMax = sum(singleKmaps_px2keep,3);
for Ki = 1:nK
    map = sum(singleKmaps(:,:,:,Ki), 3);
    map = 1 - mat2gray(map);
    savename = sprintf('heatmap_cluster%02d.tif',Ki);
    imwrite(map, savename)
    
    map = sum(singleKmaps(:,:,:,Ki), 3) ./ alMaps_pixelMax;
    map(~isfinite(map)) = 0;
    map = 1 - mat2gray(map);
    savename = sprintf('norm_heatmap_cluster%02d.tif',Ki);
    imwrite(map, savename)
end


% make more saturated images
for Ki = 1:nK  
    map = sum(singleKmaps(:,:,:,Ki), 3) ./ alMaps_pixelMax;
    map(~isfinite(map)) = 0;
    map = map.^(1/4);
    map = 1 - mat2gray(map);
    savename = sprintf('norm_heatmap_cluster%02d.tif',Ki);
    imwrite(map, savename)
end


%% tentative colormap
x = 0:0.001:1;
x = x.^(2/5);
x = 1-x;
x = repmat(x, 100, 1);
figure; imshow(x)
export_fig(fullfile(folder2figure, 'tentativeColormap_normHeatMaps.tif'))

%% plot the two distributions (AMMC vs WED)
alClustCount_splits = squeeze(sum(alClustCount, 1));
alClustCount_splits = alClustCount_splits(:,[1 3 2 4]);  % IPSI AMMC, ipsi wed, CONTRA AMMC, contra wed
alClustCount_splits = bsxfun(@rdivide, alClustCount_splits, reshape(repmat(max(reshape(alClustCount_splits, [], 2)), 2,1), 1,4)  );

% % solution #1 - pool ipsi and contra only for WED, keep only ipsi AMMC
% distrAMMC = alClustCount_splits(:,1);
% distrWED = sum(alClustCount_splits(:,[2,4]),2);
% 
% figure; hold on, plot(0:1,0:1, '-k'); axis image
% scatter(distrAMMC, distrWED, 120, cmap, 'filled')
% h = gca;
% h.TickDir = 'out';
% h.XTick = 0:0.25:1;
% h.YTick = 0:0.25:1;
% xlabel('extent into ipsilateral AMMC')
% ylabel('extent into ipsi- or contra- WED')

% solution #2 - pool ipsi and contra only for WED, ipsi and contra AMMC
distrAMMC = sum(alClustCount_splits(:,[1,3]),2);
distrWED = sum(alClustCount_splits(:,[2,4]),2);

figure; hold on, plot(0:2,0:2, '-k'); axis image
scatter(distrAMMC, distrWED, 175, cmap, 'filled')
h = gca;
h.TickDir = 'out';
h.XTick = 0:0.25:1.5;
h.YTick = 0:0.25:1.5;
xlabel('extent into ipsi- or contra- AMMC')
ylabel('extent into ipsi- or contra- WED')
% plot([0,2],[0,1], '--r')
% plot([0,1],[0,2], '--b')
% plot([0,2],[0,0.5], '.-r')
% plot([0,0.5],[0,2], '-.b')
% plot([0,0.25],[0,2], ':b')
h.XLim = [0, 1.1];
h.YLim = [0, 1.1];
export_fig(fullfile(folder2figure, 'scatterAMMCWED.eps'))

% % solution #3 - pool ipsi and contra plus contra AMMC for WED, vs ipsi AMMC only
% distrAMMC = sum(alClustCount_splits(:,[1]),2);
% distrWED = sum(alClustCount_splits(:,[2:4]),2);
% 
% figure; hold on, plot(0:1,0:1, '-k'); axis image
% scatter(distrAMMC, distrWED, 120, cmap, 'filled')
% h = gca;
% h.TickDir = 'out';
% h.XTick = 0:0.25:1;
% h.YTick = 0:0.25:1;
% xlabel('extent into ipsi- or contra- AMMC')
% ylabel('extent into ipsi- or contra- WED')


%% sort tuning curves
clustersAMMC = find(distrAMMC > distrWED);
clustersWED = find(distrAMMC < distrWED);
TCnorm = bsxfun(@rdivide, TCfinal, max(TCfinal,[],2) );
% [xvar, ia, ic] = unique([chunks.pips.carrierFreqLevels; chunks.tones.carrierFreqLevels(1)]);




crop_N_lastFreqs = 0;
[hfig, hax] = figureTracesI_PP_pixels( 2, 400, 100 );

% AMMC clusters
Kgroup = 1;
axes(hax(Kgroup)), hold on
for Ki = clustersAMMC'
    plot(xvar(1:end-crop_N_lastFreqs), TCnorm(Ki,1:end-crop_N_lastFreqs), 'Color', cmap(Ki,:))
end
% plot([xvar(1), xvar(end-crop_N_lastFreqs)],[0 0], '-k'); 

hax(Kgroup).XLim = [xvar(1), xvar(end-crop_N_lastFreqs)];
hax(Kgroup).YLim = [min( [0, min(min(TCnorm(clustersAMMC,:)))] ), 1 ];
hax(Kgroup).XAxis.Visible = 'on';
hax(Kgroup).XTick = xvar(1:end-crop_N_lastFreqs);
hax(Kgroup).YTick = 0:0.5:1;
hax(Kgroup).XAxis.Scale = 'log';
hax(Kgroup).XMinorGrid = 'off';
hax(Kgroup).XMinorTick = 'off';
hax(Kgroup).XGrid = 'off';
hax(Kgroup).XTickLabel  = [];
% hax(Kgroup).XTickLabel = xvar(1:end-crop_N_lastFreqs);


% WED clusters
Kgroup = 2;
axes(hax(Kgroup)), hold on
for Ki = clustersWED'
    plot(xvar(1:end-crop_N_lastFreqs), TCnorm(Ki,1:end-crop_N_lastFreqs), 'Color', cmap(Ki,:))
end
% plot([xvar(1), xvar(end-crop_N_lastFreqs)],[0 0], '-k'); 

hax(Kgroup).XLim = [xvar(1), xvar(end-crop_N_lastFreqs)];
hax(Kgroup).YLim = [min( [0, min(min(TCnorm(clustersWED,:)))] ), 1 ];
hax(Kgroup).XAxis.Visible = 'on';
hax(Kgroup).XTick = xvar(1:end-crop_N_lastFreqs);
hax(Kgroup).YTick = 0:0.5:1;
hax(Kgroup).XAxis.Scale = 'log';
hax(Kgroup).XMinorGrid = 'off';
hax(Kgroup).XMinorTick = 'off';
hax(Kgroup).XGrid = 'off';
hax(Kgroup).XTickLabel = xvar(1:end-crop_N_lastFreqs);
hax(Kgroup).XTickLabelRotation = 90;

export_fig(fullfile(folder2figure, 'normalizedTuningCurves_AMMCvsWED.eps'))


%%
folder2figure = '/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/figure4';
Z = 4;
map = singleKmaps_px2keep(:,:,Z);
map = logical(map);
map2 = imfill(map, 'holes');
se = strel('disk',10);
map3 = imclose(map2, se);
figure; imshow(map3)

se = strel('line',7,0);
map4 = imerode(map3, se);
se = strel('line',7,90);
map4 = imerode(map4, se);
se = strel('disk',3);
map4 = imdilate(map4, se);
se = strel('line',7,45);
map4 = imerode(map4, se);
se = strel('line',7,135);
map4 = imerode(map4, se);
figure; imshow(map4)

savename = fullfile(folder2figure, sprintf('binaryMap_AMMC_ROI_z%02d.tif',Z));
    imwrite(map4, savename)


%% remake heatmaps as real conditional probability maps
load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/Downsampled/rawDataLinkage/ward_ZSCpix/19clusters/clusterHeatMaps/alignedData.mat');
fall = '/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/Downsampled/rawDataLinkage/ward_ZSCpix/19clusters/clusterHeatMaps/allPlanes/conditionalP';
mkdir(fall)
cd(fall)


% alMaps_pixelMax = sum(singleKmaps_px2keep,3);
prob_px2keep = mean(singleKmaps_px2keep,3);
for Ki = 1:nK
    map = mean(singleKmaps(:,:,:,Ki), 3) .* prob_px2keep;
    map = map.^(1/4);
%     figure; imshow(map, []), title('condit prob K==1 map - .^(1/3)', 'Interpreter', 'none'); colorbar
%     colormap(flipud(gray))
    map = 1 - mat2gray(map);
    
    savename = sprintf('condProb_heatmap_cluster%02d.tif',Ki);
    imwrite(map, savename)
    savename = sprintf('norm_heatmap_cluster%02d.tif',Ki);
    imwrite(map, savename)
end


fall = '/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/Downsampled/rawDataLinkage/ward_ZSCpix/19clusters/clusterHeatMaps/allPlanes/regularProbabilityHeatmap';
mkdir(fall)
cd(fall)

for Ki = 1:nK
    map = mean(singleKmaps(:,:,:,Ki), 3);
    map = map.^(1/4);
%     figure; imshow(map, []), title('condit prob K==1 map - .^(1/3)', 'Interpreter', 'none'); colorbar
%     colormap(flipud(gray))
    map = 1 - mat2gray(map);
    
    savename = sprintf('regProb_heatmap_cluster%02d.tif',Ki);
    imwrite(map, savename)
    savename = sprintf('norm_heatmap_cluster%02d.tif',Ki);
    imwrite(map, savename)
end


%% remake normalization and splitting
load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/Downsampled/rawDataLinkage/ward_ZSCpix/19clusters/clusterHeatMaps/alignedData.mat');
% 'alClustCount', 'singleKmaps', 'singleKmaps_px2keep'
alClustCount_shuf = squeeze(sum(alClustCount, 1));
alClustCount_shuf = alClustCount_shuf(:,[1 3 2 4]);  % IPSI AMMC, ipsi wed, CONTRA AMMC, contra wed
alClust_shuf_mean(:,1:2) = alClustCount_shuf(:,1:2) / sum(dataTable.ipsilateral);
alClust_shuf_mean(:,3:4) = alClustCount_shuf(:,3:4) / sum(dataTable.ipsilateral==0);

%% do not norm - just average z planes counted in each case (ipsi vs contra)
% maxProjSigPx =  max(singleKmaps_px2keep,[],3); %across all z only gives you ipsi vs contra
% imwrite(logical(maxProjSigPx), fullfile(folder2figure, 'maxProjection_sigPixels.tif'))
maxProjSigPx = false(size(singleKmaps_px2keep,1), size(singleKmaps_px2keep,2), 4);
for z = 1: size(dataTable,1)
    % apply corresponding mask and count
    thisMZ = dataTable.matrixsZ(z) + 11;
    al2fly119_mapClean = logical(singleKmaps_px2keep(:,:,z));
    if dataTable.ipsilateral(z)
        maxProjSigPx(:,:,1) =  maxProjSigPx(:,:,1) | (al2fly119_mapClean & alignedAMMCmasks(:,:,thisMZ)) ; %AMMC count ipsi (1)
        maxProjSigPx(:,:,3) =  maxProjSigPx(:,:,3) | (al2fly119_mapClean & ~alignedAMMCmasks(:,:,thisMZ)) ; %ipsi notAMMC count
    else
        maxProjSigPx(:,:,2) = maxProjSigPx(:,:,2) | (al2fly119_mapClean & alignedAMMCmasks(:,:,thisMZ)); %AMMC count contra (2)
        maxProjSigPx(:,:,4) = maxProjSigPx(:,:,4) | (al2fly119_mapClean & ~alignedAMMCmasks(:,:,thisMZ)); %contra notAMMC count
    end
end
figure;
for z = 1:4
    subplot(2,2,z)
    imshow(maxProjSigPx(:,:,z))
    maxProjSigPx_counts(z) = sum(sum(maxProjSigPx(:,:,z)));
    title(maxProjSigPx_counts(z))
end
export_fig(fullfile(folder2figure, 'maxProjection_sigPixels_splits.eps'))

figure;
subplot(1,3,1)
bar(maxProjSigPx_counts([1 3 2 4]))
axis tight
set(gca, 'XTickLabel', {'Ai', 'Wi', 'Ac', 'Wc'})
title('max proj area')

subplot(1,3,2)
bar(round(sum(alClustCount_shuf)))
axis tight
set(gca, 'XTickLabel', {'Ai', 'Wi', 'Ac', 'Wc'})
title('total count across zs')

subplot(1,3,3)
bar(round(sum(alClust_shuf_mean)))
axis tight
set(gca, 'XTickLabel', {'Ai', 'Wi', 'Ac', 'Wc'})
title('mean count over z')
export_fig(fullfile(folder2figure, 'maxProjection_sum_mean_counts_splits.eps'))


% scatter plot - area of registered pixels
alClustCount_splits = alClust_shuf_mean;    % IPSI AMMC, ipsi wed, CONTRA AMMC, contra wed
distrAMMC = sum(alClustCount_splits(:,[1,3]),2);
distrWED = sum(alClustCount_splits(:,[2,4]),2);

maxV = max(cat(1, distrAMMC, distrWED));

figure; hold on, plot(0:maxV,0:maxV, '-k'); axis image
scatter(distrAMMC, distrWED, 175, cmap, 'filled')
h = gca;
h.TickDir = 'out';
% h.XTick = 0:0.25:1;
% h.YTick = 0:0.25:1;
xlabel(sprintf('representation in the AMMC\n(average registered pixels)'))
ylabel(sprintf('representation in the WED\n(average registered pixels)'))
% plot([0,2],[0,1], '--r')
% plot([0,1],[0,2], '--b')
% plot([0,2],[0,0.5], '.-r')
% plot([0,0.5],[0,2], '-.b')
% plot([0,0.25],[0,2], ':b')
h.XLim = [0, maxV];
h.YLim = [0, maxV];
export_fig(fullfile(folder2figure, 'scatterAMMCWED_regPixels.eps'))


% remake bar plots as well
% for alClustCount_splits, stackData needs to have:
% size(stackData,1) == nK          % n groups per XAxis
% (stackdata, 2):  [ ipsi, contra] % within a stacked bar
% (stackdata, 3):  [ A ,W]         % bars within a group
stackData = alClustCount_splits; % IPSI AMMC, ipsi wed, CONTRA AMMC, contra wed
stackData = permute(stackData, [1,3,2,4]);
stackData = reshape(stackData, [],2,2 );

groupLabels = num2cell(1:19);
hb = plotBarStackGroups(stackData, groupLabels);
hb(1,1).FaceColor = [1 0.4 0.4];   %AMMC ipsi
hb(2,1).FaceColor = [0.4 0.4 1];   %wed ipsi
hb(1,2).FaceColor = [0.6  0  0];   %AMMC contra
hb(2,2).FaceColor = [0  0  0.6];   %wed contra
export_fig(fullfile(folder2figure, 'clusterDistrib_AMMC_WED_ipsicontraStacked_means.eps'))


%% now try cluster-spec max proj
%% do not norm - just average z planes counted in each case (ipsi vs contra)
maxProjK = false(size(singleKmaps,1), size(singleKmaps,2), 4, nK);
for Ki = 1:nK
    for z = 1: size(dataTable,1)
        % apply corresponding mask and count
        thisMZ = dataTable.matrixsZ(z) + 11;
        al2fly119_mapClean = logical(singleKmaps(:,:,z,Ki));
        if dataTable.ipsilateral(z)
            maxProjK(:,:,1, Ki) =  maxProjK(:,:,1, Ki) | (al2fly119_mapClean & alignedAMMCmasks(:,:,thisMZ)) ; %AMMC count ipsi (1)
            maxProjK(:,:,3, Ki) =  maxProjK(:,:,3, Ki) | (al2fly119_mapClean & ~alignedAMMCmasks(:,:,thisMZ)) ; %ipsi notAMMC count
        else
            maxProjK(:,:,2, Ki) = maxProjK(:,:,2, Ki) | (al2fly119_mapClean & alignedAMMCmasks(:,:,thisMZ)); %AMMC count contra (2)
            maxProjK(:,:,4, Ki) = maxProjK(:,:,4, Ki) | (al2fly119_mapClean & ~alignedAMMCmasks(:,:,thisMZ)); %contra notAMMC count
        end
    end
end


for Ki = 1:nK
%     %maps
    mP_wholeK = logical(sum(maxProjK(:,:,:,Ki),3));
    fname = fullfile('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/figure5n/mp_maps', sprintf('mp_wholeMap_k%02d.tif',Ki));
    imwrite(~mP_wholeK, fname)
% end
%     mP_AMMC_K = logical(sum(maxProjK(:,:,1:2,Ki),3));
%     mP_WED_K = logical(sum(maxProjK(:,:,3:4,Ki),3));
%     %save to disk 
    
    % counts
    alClust_shuf_MP(Ki,1) = sum(sum(maxProjK(:,:,1,Ki)));
    alClust_shuf_MP(Ki,2) = sum(sum(maxProjK(:,:,3,Ki)));
    alClust_shuf_MP(Ki,3) = sum(sum(maxProjK(:,:,2,Ki)));
    alClust_shuf_MP(Ki,4) = sum(sum(maxProjK(:,:,4,Ki)));
end




% scatter plot - area of registered pixels
alClustCount_splits = alClust_shuf_MP;    % IPSI AMMC, ipsi wed, CONTRA AMMC, contra wed
distrAMMC = sum(alClustCount_splits(:,[1,3]),2);
distrWED = sum(alClustCount_splits(:,[2,4]),2);

maxV = max(cat(1, distrAMMC, distrWED));

figure; hold on, plot(0:maxV,0:maxV, '-k'); axis image
scatter(distrAMMC, distrWED, 175, cmap, 'filled')
h = gca;
h.TickDir = 'out';
% h.XTick = 0:0.25:1;
% h.YTick = 0:0.25:1;
xlabel(sprintf('representation in the AMMC\n(max proj pixels)'))
ylabel(sprintf('representation in the WED\n(max proj pixels)'))
% plot([0,2],[0,1], '--r')
% plot([0,1],[0,2], '--b')
% plot([0,2],[0,0.5], '.-r')
% plot([0,0.5],[0,2], '-.b')
% plot([0,0.25],[0,2], ':b')
h.XLim = [0, maxV];
h.YLim = [0, maxV];
export_fig(fullfile(folder2figure, 'scatterAMMCWED_MAXPROJ_Pixels.eps'))


% remake bar plots as well
% for alClustCount_splits, stackData needs to have:
% size(stackData,1) == nK          % n groups per XAxis
% (stackdata, 2):  [ ipsi, contra] % within a stacked bar
% (stackdata, 3):  [ A ,W]         % bars within a group
stackData = alClustCount_splits; % IPSI AMMC, ipsi wed, CONTRA AMMC, contra wed
stackData = permute(stackData, [1,3,2,4]);
stackData = reshape(stackData, [],2,2 );

groupLabels = num2cell(1:19);
hb = plotBarStackGroups(stackData, groupLabels);
hb(1,1).FaceColor = [1 0.4 0.4];   %AMMC ipsi
hb(2,1).FaceColor = [0.4 0.4 1];   %wed ipsi
hb(1,2).FaceColor = [0.6  0  0];   %AMMC contra
hb(2,2).FaceColor = [0  0  0.6];   %wed contra
export_fig(fullfile(folder2figure, 'clusterDistrib_AMMC_WED_ipsicontraStacked_MAXPROJ.eps'))




