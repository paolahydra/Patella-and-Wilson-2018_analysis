clearvars -except aZ R pixel2keep iKeep pxKeep
load_panSpec_downsampled;
if ~exist('R', 'var')
    load('R_matrix_Downsampled_smallWindow.mat', 'R', 'pixel2keep', 'iKeep', 'pxKeep')
    R = R(pxKeep,:);
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
Rzs = zscore(R, 0, 2);
clusterFolder = '/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix';
treeFileName = 'output_ward_movAVG_ZSCpix_linkage.mat';
if exist(clusterFolder) ~= 7
    mkdir(clusterFolder)
end
cd(clusterFolder)
load('klust_maxClust10.mat')
assert(size(R,1) == sum(pxKeep))


load('/Users/galileo/Dropbox (HMS)/Data/cmap19_spectral.mat', 'cmap')
cmap10 = cmap([2 1 19 17 14 12 10 8 6 5],:); %sorted from bottom to up as in dendrogram
save('/Users/galileo/Dropbox (HMS)/Data/cmap19_spectral.mat', 'cmap10', '-append')
colors = cmap10;

clusterOrder = [10 9 8 4 1 2 3 7 5 6 ]; % K2use!!

%% prepare R matrix for refined linkage calculation
% % remove some noise
% tic, [W,Y,latent,~,explained,mu] = pca(R, 'Algorithm', 'eig', 'Centered', true, 'VariableWeights','variance'); toc
% % 
% figure; subplot(2,1,1)
% plot(latent(1:100))
% ylabel('eigenvalues')
% subplot(2,1,2)
% expld = cumsum(explained(1:100));
% plot(expld)
% ylabel('cumulative explained variance')
% pca_dims = 50;      %less would be enough
% for ch_i = 16 : 6 : pca_dims
%     figure; hold on
%     for i = 1:6
%         subplot(6, 1, i)
%         plot(W(:, ch_i + i - 1))
%         ax = gca;
%         ax.XGrid = 'on';
%         ax.XLim = [0, size(R,2)];
%         ylabel(sprintf('PC %d', ch_i + i - 1))
%     end
% end
% pca_dims = 50;
% recR = Y(:,1:pca_dims)*W(:,1:pca_dims)';
% recR = bsxfun(@plus, recR, mu);
% recR = zscore(recR, 0, 2);
% save('R_matrix_Downsampled_noNoisyPxs_PCA_centVar_50pcs_cleaned_zscoredByPixel.mat', 'recR')
% 
% % alternatively, try using a moving mean to further smooth the signal
% recR2 = (movmean(R', 9))';
% recR2 = zscore(recR2, 0, 2);
% save('R_matrix_Downsampled_noNoisyPxs_movAvg9_zscoredByPixel.mat', 'recR2')


%% linkage was calculated on orchestra
load(fullfile(clusterFolder, treeFileName), 'Zc');
maxClust = 11;
maxClustComponents = 150;
maxLeavesDend = maxClust;
distanceMetrics = 'euclidean';
linkageMethod = 'ward'; 
% temporary dendr
order11 = [7 8 9 6 5 4 3 10 11 1 2];
outpermLeaves11 = [5    11     1    10     4     2     3     6     7     8     9];
figure; [~,Tdend,outpermLeaves] = dendrogram(Zc, maxLeavesDend, 'Orientation', 'left');
% export_fig(sprintf('dendrogram_maxLeaves%d_originalLabels.pdf',maxLeavesDend))
% will resort the dendrogram manually, and this will be the corresponding
% outpermLeaves
outpermLeaves = fliplr(outpermLeaves11(order11)); % I am not changing the cluster identity here, just moving branches around (physically in the dendrogram)

% axDend = gca;
% axDend.YTickLabel = 1:maxLeavesDend;

%% cluster (main), sort klusts based on dendrogram -- AND THEN RE-SORT again [DO]
Tparent = cluster(Zc,'maxclust',maxClust);

% GOAL: Ki order of superclusts is fixed. Starting from the left-most
% cluster-branch in dendrogram, we need to figure out which klust(Ki) it
% corresponds to.
for Ki = 1:length(unique(Tparent))
    tempKs = find(Tparent==Ki);
    klust(Ki).k = (tempKs(:))';
end
for Ki = 1:length(unique(Tparent)), allElements(Ki, 1:length(klust(Ki).k)) = klust(Ki).k; end

countLeaves = 1;
for Ki = 1:length(unique(Tparent))   %while ~isempty(leavesorder)
    lengthLeave = 0;
    leavesorder = find(Tdend == outpermLeaves(countLeaves));
    countLeaves = countLeaves+1;
    lengthLeave = lengthLeave + length(leavesorder);
    firstLeaf(Ki) = leavesorder(1);
    [dendrOrder(Ki), ~] = find(allElements==firstLeaf(Ki));
    while sum(Tparent == dendrOrder(Ki)) > lengthLeave
        leavesorder = find(Tdend == outpermLeaves(countLeaves));
        countLeaves = countLeaves+1;
        lengthLeave = lengthLeave + length(leavesorder);
    end
end
clear leavesorder countLeaves leavesorder

% Also rename T once and for all according to dendrogram
T = zeros(size(Tparent));
for Ki = 1:length(unique(Tparent))
    T(Tparent==dendrOrder(Ki)) = Ki;
    tempKs = find(T==Ki);
    klustSorted(Ki).k = (tempKs(:))';
end
Tparent = T;
klust = klustSorted;
dendrOrder = 1:length(unique(Tparent));
clear klustSorted T
clear outpermLeaves %to avoid errors
save(sprintf('klust_maxClust%d.mat', maxClust), 'klust', 'dendrOrder', 'Tparent');
nK = length(unique(Tparent));
cmap10 = distinguishable_colors(nK);
cmap10 = flipud(cmap10);
% tColors = distinguishable_colors(14);
% colors(7,:) = tColors(end,:);
save(sprintf('klust_maxClust%d.mat', maxClust), 'colors', '-append')


%% cluster components
Tcomponents = cluster(Zc,'maxclust',maxClustComponents);
[table_T10030] = crosstab(Tcomponents, Tparent);
clust30 = remapIndexedPointsToRun_RoisMaps(Tparent, pixel2keep, sizeMap, pxKeep);
clust100 = remapIndexedPointsToRun_RoisMaps(Tcomponents, pixel2keep, sizeMap, pxKeep);
sortFliesGenoType = 1;
[~, flyClusterCountDivZs100] = barcluster_panSpec(Tcomponents, datalist, flyfoldersUnique, clust100.ROIs, sortFliesGenoType);
[~, flyClusterCountDivZs30] = barcluster_panSpec(Tparent, datalist, flyfoldersUnique, clust30.ROIs, sortFliesGenoType);

figure; imagesc(table_T10030)
axis square
ax = gca;
ax.XTick = 1:nK;
ax.XTickLabel = sum(table_T10030);
ax.XTickLabelRotation = 90;
export_fig('table_main_components.png')
%% display component traces
sortedflyMeta = {'pan' 'pan' 'pan' 'pan' 'pan' 'pan' 'A22' 'A23' 'A23' 'A26' 'AB15' 'AB15' 'AB15' 'BA28' 'B2' 'B2' 'ACE4' 'ACE4' 'ACE4' 'CE32' 'CE32' 'AD29'  'AD29' 'AD29'};
for k = 1:nK
% [hfig, hax] = figureTracesI_PP( sum(table_T10030(:,k)>0)+1, [0.66 0.34] );
[hfig, hax] = figureTracesI_PP_pixels( sum(table_T10030(:,k)>0)+1, [600 300], 30 );
axes(hax(1,1)), hold on
allR = R(Tparent==k,:);
avg = mean(allR,1);
upperbound = prctile(allR, 75, 1);
lowerbound = prctile(allR, 25, 1);
plot(avg, '-m', 'LineWidth', 2), hold on, axis tight
plot(upperbound, '-k')
plot(lowerbound, '-k')
axis tight
ylabel(sprintf('T30 k %d\n(%d)',k, sum(Tparent==k)))

axes(hax(1,2))
bar(1:NaF, flyClusterCountDivZs30(:,k))
hax(1,2).XTick = 1:NaF;
hax(1,2).FontSize = 7;
hax(1,2).XTickLabel = flyNumUnique(idx_sortingFlies);
hax(1,2).XTickLabelRotation = 90;


compNums = find(table_T10030(:,k));
for i = 1 : sum(table_T10030(:,k)>0) %number of component clusters
    axes(hax(i+1,1)), hold on
    allR = R(Tcomponents==compNums(i),:);
    avg = mean(allR,1);
    upperbound = prctile(allR, 75, 1);
    lowerbound = prctile(allR, 25, 1);
    plot(avg, '-r', 'LineWidth', 2), hold on, axis tight
    plot(upperbound, '-k')
    plot(lowerbound, '-k')
    axis tight
    ylabel(sprintf('T60 k %d\n(%d)',compNums(i), sum(Tcomponents==compNums(i))))
    drawnow
    
    axes(hax(i+1,2))
    bar(1:NaF, flyClusterCountDivZs100(:,compNums(i)))
    hax(i+1,2).XTick = 1:NaF;
    hax(i+1,2).FontSize = 7;
    hax(i+1,2).XTickLabel = [];
    hax(i+1,2).XTickLabelRotation = 90;
end
hax(i+1,2).XTickLabel = sortedflyMeta; 

export_fig(sprintf('clust30_%02d_R_ResponseAndDistribution.pdf',k))
close
end

% hhist = figure; h = histogram(categorical(T));
% % compare the two outputs
% [table_TMovavgTpc, chi2_TMovavgTpc, p_TMovavgTpc] = crosstab(TMovAvg, T);
% for k = 1:30
%     disp(k)
%     pw = pdist(Rzs(TMovAvg==k,:));
%     wr_mean(k) = mean(pw);
%     wr_median(k) = median(pw);
%     wr_total(k) = sum(pw);
% end
% for k = 1:30
%     disp(k)
%     pw = pdist(Rzs(T==k,:));
%     wrPC_mean(k) = mean(pw);
%     wrPC_median(k) = median(pw);
%     wrPC_total(k) = sum(pw);
% end
% 
% % compare
% fprintf('mean(movAvg_total): \t %d \n', mean(wrPC_total) > mean(wr_total))
% fprintf('mean(movAvg_mean): \t %d \n', mean(wrPC_mean) > mean(wr_mean))
% fprintf('mean(movAvg_median): \t %d \n', mean(wrPC_median) > mean(wr_median))
% % RESULT: movAvg wins. Is that so surprising? All I wanted to do was avoid
% % to overfit noise... 



%% plot all average traces with errors, sorted from bottom to up as in dendrogram! [may skip if already done]
[hfig, hax] = figureTracesI_PP( length(unique(Tparent)), 1 );
hax = flipud(hax); % index 1 is bottom axis
% hfig.Position   =   [1,5,1800,700];
for Ki = 1:length(unique(Tparent))
    axes(hax(Ki)), hold on
%     hax(Ki).XTick = 0:38:228;
    hax(Ki).XGrid = 'on';
    hax(Ki).FontSize = 3;
    hax(Ki).YLabel.FontSize = 15;
K2use = clusterOrder(Ki);
    allR = R(Tparent==K2use,:);
    avg = mean(allR,1);
%     upperbound = prctile(allR, 75, 1);
%     lowerbound = prctile(allR, 25, 1);
    plot(avg, 'Color', cmap10(Ki,:)), hold on, axis tight
%     plot(upperbound, '-k')
%     plot(lowerbound, '-k')
    axis tight
    ylabel(Ki)
    ylabel(sprintf('clust %d\n(px %d)',Ki, sum(Tparent==Ki)))
    drawnow
end
hfig.Position = [ 0    0.0044    0.3257    0.8900];
export_fig(sprintf('clustMain_d_mtraces_resorted.pdf'))

%% make maps
for z = 1:length(datalist)
    [runfolders{z}, ~] = fileparts(datalist{z});
end
clust = remapIndexedPointsToRun_RoisMaps(Tparent, pixel2keep, sizeMap);

% make maps (divided old/new because different dimensions)
klust = makeplot_cMaps_superKlusts_singleIteration(clust, klust, dendrOrder);
folder = 'alignedMaps';
mkdir(folder)
cd(folder)

nK = length(klust);

% Col = cat(1,[1,1,1],distinguishable_colors(nK));
% cmap10 = cmap11([2,3,1,4:11],:); %new cmap 20170914
Col = cat(1,[1,1,1],cmap10);

regGeneralData = matfile('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/anatomyImages/anatomy_alignment_metadata.mat');
Zclust = 1;
for f = 1 : NaF

    flyfolder = flyfoldersUnique{f};
    ZZ = strfind(datalist, flyfolder);
    zs_fly = zeros(size(ZZ));
    for i = 1:length(zs_fly)
        if ~isempty(ZZ{i})
            zs_fly(i) = 1; % indices into datalist relative to fly f
        end
    end
    zs_fly = logical(zs_fly);

    % load fly's stack2fly150 transformation
    stack2fly150_info = matfile(fullfile(flyfolder, ['stack2fly150sStack_' basenames{find(zs_fly,1)}(1:6) '.mat']));
    if ~isfield(stack2fly150_info, 'Roriginal') && ~strcmp(basenames{find(zs_fly,1)}(1:6), 'fly150')
        stack2fly150_info.Properties.Writable = true;
        stack2fly150_info.Roriginal = imref2d(size(stack2fly150_info.recovered));
    end
    %% sort runs from ventral to dorsal
    zetas = find(zs_fly);
    sliceNums = regGeneralData.alignedSliceNumbers(1,zetas);
    [sliceNums,b] = sort(sliceNums);
    zetas = zetas(b);

    %% load specific runs data and transformations and compute, correcting for downsampling
    for iz = 1:sum(zs_fly)
        z = zetas(iz);
        sliceNum = sliceNums(iz);
        % load specific run's transformation
        regData = matfile(fullfile(runfolders{z}, ['stackRegistration_' basenames{z} '.mat']));
        if ~isfield(regData, 'Roriginal')
            regData.Properties.Writable = true;
            regData.Roriginal = imref2d(size(regData.target));
        end

        al2fly150_mapAll = [];
        for Ki = 1:length(unique(Tparent))
            K2use = clusterOrder(Ki);
            mapClean = klust(K2use).mapsZold(:,:,z);
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
            mapClean(mapClean>=0.1) = 1;                        %note 1 and not K2use!

            mapClean = imrotate(mapClean, regGeneralData.totAng(1,z)); %preliminary transformation
            recovered_mapClean = imwarp(mapClean,regData.tform,'OutputView',regData.Roriginal,'FillValues',0); %indexed image
            if strcmp(basenames{find(zs_fly,1)}(1:6), 'fly150')
                al2fly150_mapClean = recovered_mapClean;
            else
                al2fly150_mapClean = imwarp(recovered_mapClean,stack2fly150_info.tform,'OutputView',stack2fly150_info.Roriginal,'FillValues',0);
            end
            al2fly150_mapClean(al2fly150_mapClean<=0.5) = 0;
            al2fly150_mapClean(al2fly150_mapClean>=0.1) = Ki; %dendrogram is colored in order from 1st to nK-esim  in colormap. Clusters are colored accordingy.
            if isempty(al2fly150_mapAll)
                al2fly150_mapAll = al2fly150_mapClean;
            else
                al2fly150_mapAll(al2fly150_mapClean == Ki) = Ki; %dendrogram is colored in order from 1st to nK-esim  in colormap. Clusters are colored accordingy.
            end
        end

        AllClusterMaps(:,:,z) = al2fly150_mapAll;

        AlphaMatrix = zeros(size(al2fly150_mapAll(:,:,1)));
        AlphaMatrix(al2fly150_mapAll(:,:,1) > 0) = 1;

        al2fly150_mapAll = ind2rgb(uint8(al2fly150_mapAll), Col);
        savename = sprintf('alignALLklustsMap_fly%02d_ordZ%02d_%s.png', f, z, basenames{z});
        imwrite(al2fly150_mapAll, savename, 'Alpha', double(AlphaMatrix))
    end
end
cd ..
save(sprintf('alignedMatrix_Images_maxClust%d.mat',nK), 'AllClusterMaps', 'Col')



%% calculate total within-clusts-distance for this model
clearvars -except Rzs Tparent nK maxClust R pixel2keep iKeep pxKeep
for k = 1:nK
    disp(k)
    pw = pdist(Rzs(Tparent==k,:));
    wr_total(k) = sum(pw);
end
sum(wr_total)
save(sprintf('klust_maxClust%d.mat', maxClust), 'wr_total', '-append');
% save(sprintf('klust_maxClust%d.mat', nK), 'wr_total', '-append');


%% try recomputing linkage and dendrogram based on cluster centroids
for Ki = 1:length(unique(Tparent))
    centroids(Ki, :) = mean(R(Tparent==Ki,:)); 
end
%recompute linkage
Czs = zscore(centroids, 0, 2);
D = pdist(Czs);

map2 = brewermap(30, 'PuBuGn'); %'OrRd');
map2 = flipud(map2);
figure; imagesc(squareform(D)); axis image
colormap(map2)
colorbar
axis xy

CZs_tree = linkage(D, 'complete');
leafOrder = optimalleaforder(CZs_tree,D);
figure
[~, Tdend, outpermL] = dendrogram(CZs_tree,'reorder',leafOrder, 'Orientation', 'left');
title('complete linkage - Optimal Leaf Order')


CZs_tree = linkage(D, 'centroid');
leafOrder = optimalleaforder(CZs_tree,D);
figure
[~, Tdend, outpermL] = dendrogram(CZs_tree,'reorder',leafOrder, 'Orientation', 'left');
title('centroid linkage - Optimal Leaf Order')

CZs_tree = linkage(D, 'single');
leafOrder = optimalleaforder(CZs_tree,D);
figure
[~, Tdend, outpermL] = dendrogram(CZs_tree,'reorder',leafOrder, 'Orientation', 'left');
title('single linkage - Optimal Leaf Order')

%% Part I -- identify a threshold that separates best noisy (to be removed) from good SNR pixels. (iteratively from high to low)
% % pxKeep = true(size(T));   % (first time only)
% % iteratio = 1;  % (first time only)
% distance = 0.85; 
% T = cluster(Zc,'cutoff', distance, 'criterion', 'distance');  
% selClusts = unique(T(pxKeep));
% nK = length(selClusts);
% percPixelsRetained = sum(pxKeep) / size(R,1) * 100;
% 
% hhist = figure; h = histogram(categorical(T(pxKeep)));
% countSelClusters = h.Values;
% ax = gca;
% title(sprintf('Iter %d: Retained %3.1f%% pixels, divided into %d clusters by threshold distance %.3f', iteratio, percPixelsRetained, nK, distance))
% 
% threshcount = 20;
% selClusts = selClusts(h.Values>threshcount);
% nK = length(selClusts)
% countSelClusters = h.Values(h.Values>threshcount);
% 
% fprintf('Iter %d: Retained %3.1f%% pixels, divided into %d clusters,\nof which %d with more than %d pixels, by threshold distance %.3f \n', iteratio, percPixelsRetained, length(unique(T(pxKeep))), nK, threshcount, distance)
% 
% %% Part 1: plot
% [hfig, hax] = figureTracesI_PP( length(selClusts), 1 );
% hax = flipud(hax);
% clear centroids
% for ik = 1:length(selClusts)
%     k = selClusts(ik);
%     centroids(ik,:) = mean(R(T==k,:));
%     axes(hax(ik)), hold on
%     upperbound = prctile(R(T==k,:), 75, 1);
%     lowerbound = prctile(R(T==k,:), 25, 1);
%     plot(upperbound, '-k'), hold on
%     plot(lowerbound, '-k')
%     plot(centroids(ik,:), '-r')
%     axis tight
%     hax(ik).YLabel.String = sprintf('k%d\n#%d', k, countSelClusters(ik));
%     hax(ik).YLabel.FontSize = 6;
%     drawnow
% end
% title(sprintf('Iter %d: Retained %3.1f%% pixels, divided into %d clusters,\nof which %d with more than %d pixels, by threshold distance %.3f \n', iteratio, percPixelsRetained, length(unique(T(pxKeep))), nK, threshcount, distance))
% hfig.Position = [ 0    0.0044    0.3257    0.8900];
% export_fig(hfig, sprintf('CutNoise_%d_Thr%.3f.pdf', iteratio, distance))
% export_fig(hhist, sprintf('CutNoise_%d_Thr%.3f_hist.png', iteratio, distance))
% 
% %% Part 1: make a decision and reiterate (run this just once)
% save(sprintf('retainedPixels_AtIteratio%d.mat', iteratio), 'pxKeep', 'iteratio')
% 
% keepClusts = selClusts;      %change here
% 
% newPxKeep = false(size(T));
% for i = 1:length(keepClusts)
%     newPxKeep = newPxKeep | T==keepClusts(i);
% end
% pxKeep = newPxKeep;
% iteratio = iteratio+1;
% close(hfig)
% 
% %% show which pixels have been excluded
% T(pxKeep) = 1;
% T(~pxKeep) = 2;
% 
% clear klust allElements firstLeaf dendrOrder
% for Ki = 1:length(unique(T))
%     tempKs = find(T==Ki);
%     klust(Ki).k = (tempKs(:))';
% end
% nK = length(unique(T));
% colors = distinguishable_colors(nK);
% 
% for z = 1:length(datalist)
%     [runfolders{z}, ~] = fileparts(datalist{z});
% end
% [flyfolders, ~] = extractRunfoldersBasenames(datalist);
% flyfoldersUnique = unique(flyfolders, 'stable');
% NaF = length(flyfoldersUnique);
% 
% sizeMap = [60, 86];
% clust = remapIndexedPointsToRun_RoisMaps(T, pixel2keep, sizeMap);
% 
% % make maps (divided old/new because different dimensions)
% klust = makeplot_cMaps_superKlusts_singleIteration(clust, klust, 1:2);
% folder = 'noisyPixels_alignedMaps';
% mkdir(folder)
% cd(folder)
% 
% snK = length(unique(T));
% % Col = cat(1,[1,1,1],brewermap(snK, 'Set1')); % only gets up to 9 colors
% Col = cat(1,[1,1,1],distinguishable_colors(snK));
% 
% dendrOrder = 1:2;
% regGeneralData = matfile('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/anatomyImages/anatomy_alignment_metadata.mat');
% Zclust = 1;
% for f = 1 : NaF 
%     
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
%     % load fly's stack2fly150 transformation
%     stack2fly150_info = matfile(fullfile(flyfolder, ['stack2fly150sStack_' basenames{find(zs_fly,1)}(1:6) '.mat']));
%     if ~isfield(stack2fly150_info, 'Roriginal') && ~strcmp(basenames{find(zs_fly,1)}(1:6), 'fly150')
%         stack2fly150_info.Properties.Writable = true;
%         stack2fly150_info.Roriginal = imref2d(size(stack2fly150_info.recovered));
%     end
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
%         al2fly150_mapAll = [];
%         for Ki = 1:length(unique(T))
%             K2use = dendrOrder(Ki);
%             mapClean = klust(K2use).mapsZold(:,:,z);
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
%             mapClean(mapClean>=0.1) = 1;                        %note 1 and not K2use!
%             
%             mapClean = imrotate(mapClean, regGeneralData.totAng(1,z)); %preliminary transformation
%             recovered_mapClean = imwarp(mapClean,regData.tform,'OutputView',regData.Roriginal,'FillValues',0); %indexed image
%             if strcmp(basenames{find(zs_fly,1)}(1:6), 'fly150')
%                 al2fly150_mapClean = recovered_mapClean;
%             else
%                 al2fly150_mapClean = imwarp(recovered_mapClean,stack2fly150_info.tform,'OutputView',stack2fly150_info.Roriginal,'FillValues',0);
%             end
%             al2fly150_mapClean(al2fly150_mapClean<=0.5) = 0;
%             al2fly150_mapClean(al2fly150_mapClean>=0.1) = Ki; %dendrogram is colored in order from 1st to nK-esim  in colormap. Clusters are colored accordingy.
%             if isempty(al2fly150_mapAll)
%                 al2fly150_mapAll = al2fly150_mapClean;
%             else
%                 al2fly150_mapAll(al2fly150_mapClean == Ki) = Ki; %dendrogram is colored in order from 1st to nK-esim  in colormap. Clusters are colored accordingy.
%             end
%         end
%         
%         AllClusterMaps(:,:,z) = al2fly150_mapAll;
%         
%         AlphaMatrix = zeros(size(al2fly150_mapAll(:,:,1)));
%         AlphaMatrix(al2fly150_mapAll(:,:,1) > 0) = 1;
%         
%         al2fly150_mapAll = ind2rgb(uint8(al2fly150_mapAll), Col);
%         savename = sprintf('alignALLklustsMap_fly%02d_ordZ%02d_%s.png', f, z, basenames{z});
%         imwrite(al2fly150_mapAll, savename, 'Alpha', double(AlphaMatrix))
%     end
% end





% %% Part II -- iteratively lower the threshold and isolate clusters to separate further. Here clusters will NOT be sorted according to dendrogram
% level = 1;
% distance = 0.7;
% T = cluster(Zc,'cutoff', distance, 'criterion', 'distance');  % 24 clusters
% nK = length(unique(T))
% hhist = figure; h = histogram(categorical(T))
% ax = gca;
% ax.XTick = [];
% title(distance)
% pxcount = sort(h.Values);
% hcount = figure; plot(pxcount(pxcount>100))
% title(distance)%
% threshcount = 100;
% selClusts = find(h.Values>threshcount);
% length(selClusts)
% percPixelsRetained = sum(pxcount(pxcount>100)) / size(R,1) * 100;
% %%
% [hfig, hax] = figureTracesI_PP( length(selClusts), 1 );
% hax = flipud(hax);
% clear centroids
% for ik = 1:length(selClusts)
%     k = selClusts(ik);
%     centroids(ik,:) = mean(R(T==k,:));
%     axes(hax(ik)), hold on
%     upperbound = prctile(R(T==k,:), 75, 1);
%     lowerbound = prctile(R(T==k,:), 25, 1);
%     plot(upperbound, '-k'), hold on
%     plot(lowerbound, '-k')
%     plot(centroids(ik,:), '-r')
%     axis tight
%     hax(ik).YLabel.String = sprintf('k%d\n#%d', k, h.Values(ik));
%     drawnow
% end
% title(sprintf('threshold distance %.3f, clusters with at least %d pixels - %3.1f%% of pixels included', distance, threshcount, percPixelsRetained))
% hfig.Position = [ 0    0.0044    0.3257    0.8900];
% export_fig(hfig, sprintf('Thr%.3f_thershCount%d.pdf', distance, threshcount))
% export_fig(hhist, sprintf('Thr%.3f_thershCount%d_hist.png', distance, threshcount))
% %%
% 
% 
% 
% 
% 
% 
% %%
% export_fig(sprintf('Lev%d_hist.png', level))
% 
% [hfig, hax] = figureTracesI_PP( length(unique(T)), 1 );
% hax = flipud(hax);
% clear centroids
% for k = 1:nK
%     centroids(k,:) = mean(R(T==k,:));
%     axes(hax(k)), hold on
%     upperbound = prctile(R(T==k,:), 75, 1);
%     lowerbound = prctile(R(T==k,:), 25, 1);
%     plot(upperbound, '-k'), hold on
%     plot(lowerbound, '-k')
%     plot(centroids(k,:), '-r')
%     axis tight
%     ylabel(k)
%     drawnow
% end
% hfig.Position = [ 0    0.0044    0.3257    0.8900];
% export_fig(sprintf('Lev%d_mtraces.png', level))
% 
% 
% pxs = T==18;
% %% iteratively lower the threshold and isolate clusters to separate further. Here clusters will NOT be sorted according to dendrogram
% level = 2;
% distance = 0.98;
% newT = cluster(Zc,'cutoff', distance, 'criterion', 'distance'); 
% 
% T = newT( PixHier(1).pxs );
% nK = length(unique(T))
% hhist = figure; h = histogram(categorical(T))
% h.Values
% export_fig(sprintf('Lev%d_hist.png', level))
% 
% [hfig, hax] = figureTracesI_PP( length(unique(T)), 1 );
% hax = flipud(hax);
% clear centroids
% subclusters = unique(T);
% for ik = 1:nK
%     k = subclusters(ik);
%     centroids(ik,:) = mean(R(T==k,:));
%     axes(hax(ik)), hold on
%     upperbound = prctile(R(T==k,:), 75, 1);
%     lowerbound = prctile(R(T==k,:), 25, 1);
%     plot(upperbound, '-k'), hold on
%     plot(lowerbound, '-k')
%     plot(centroids(ik,:), '-r')
%     axis tight
%     ylabel(k)
%     drawnow
% end
% hfig.Position = [ 0    0.0044    0.3257    0.8900];
% export_fig(sprintf('Lev%d_mtraces.png', level))
% 
% PixHier(1).sub(1).pxs = newT==71;
% PixHier(1).sub(1).distance = distance;
% close(hfig) 
% close(hhist)
% 
% % try to automatically zoom in dendrogram
% hdend = figure;
% [~,Tdend,outpermLeaves] = dendrogram(Zc, 300, 'ColorThreshold', distance, 'Orientation', 'left');
% ax = gca;
% %goal is to restrict the ylim to previous cluster 18 (higher hierarchical
% %level)
% for l = 1: length(outpermLeaves)
%     pixelsThisLeaf = Tdend == outpermLeaves(l);
%     if sum(pixelsThisLeaf .* PixHier(1).pxs)
%        firstLeafThisCluster = l %true
%        lastLeafThisCluster = l; %check if true
%        %check next until you find the last one
%        l = l+1;
%        pixelsThisLeaf = Tdend == outpermLeaves(l);
%        while sum(pixelsThisLeaf .* PixHier(1).pxs)
%            lastLeafThisCluster = l;
%            l = l+1;
%            pixelsThisLeaf = Tdend == outpermLeaves(l);
%        end
%        lastLeafThisCluster
%        break
%     end
% end
% ax.YLim = [firstLeafThisCluster-0.5, lastLeafThisCluster+0.5];
% ax.XLim(2) = PixHier(1).distance;
% export_fig(sprintf('Lev%d_dendrogram.png', level))
% 

