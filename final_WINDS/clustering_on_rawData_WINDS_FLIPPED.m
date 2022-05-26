load_winds_downsampled;
% for the purpose of this script, what is called R is going to contain
% mirror data for the two hemispheres. (in left runs, stimuli 1 and 3 had
% been flipped):
R = Rflipped;

load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/Downsampled/rawDataLinkage/ward_ZSCpix/klust_maxClust19.mat', 'cmap')
% R_100.pixel2keep = load('/Users/galileo/Dropbox (HMS)/Data/WINDS_100/Downsampled/R_matrix_Downsampled_100.mat', 'pixel2keep');
% R_86.pixel2keep = load('/Users/galileo/Dropbox (HMS)/Data/WINDS_100/Downsampled/R_matrix_Downsampled_86.mat', 'pixel2keep');
% R_100.pixel2keep = R_100.pixel2keep.pixel2keep;
% R_86.pixel2keep = R_86.pixel2keep.pixel2keep;
runfolders = sets.analysislist;

treeFileName = 'output_ward_Rzsp_linkage_winds_FLIPPED.mat';
clusterFolder = fullfile(Folder2Save, '/rawDataLinkage'); % linkage was already calculated on orchestra
load(fullfile(clusterFolder, treeFileName), 'Zc');
cd(clusterFolder)
pxKeep = true(size(R, 1),1);   % no (easily-detected at least) uncorrelated pixels within this dataset

colors = [208,28,139; ...   %fucsia
          31,120,180; ...   %blu
          230,97,1; ...     % orange
          51,160,44]./255;  % green
      
% if exist(clusterFolder) ~= 7
%     mkdir(clusterFolder)
% end
% cd(clusterFolder)
% load('klust_maxClust19.mat'), nK = length(unique(Tparent));
% assert(size(R,1) == sum(pxKeep))
% clusterFolder = fullfile(clusterFolder, sprintf('%dclusters', nK));
% if exist(clusterFolder) ~= 7
%     mkdir(clusterFolder);
% end
% cd(clusterFolder)


%% batch-cluster for different Ks
map = brewermap(256, 'RdBu');
map(1:128,:) = [];
K1 = 7;
Kend = 16;
for K = K1:Kend
    disp(K)
    % maxClustComponents = 150;
    maxLeavesDend = K;
    figure; [~,Tdend,outpermLeaves] = dendrogram(Zc, maxLeavesDend, 'Orientation', 'left');
    axDend = gca;
    axDend.YTickLabel = 1:maxLeavesDend;
    export_fig(sprintf('dendrogram_maxLeaves%d.pdf',maxLeavesDend))
    
    %% cluster (main), sort klusts based on dendrogram -- AND THEN RE-SORT again [DO]
    Tparent = cluster(Zc,'maxclust',K);
    
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
    save(sprintf('klust_maxClust%d.mat', K), 'klust', 'dendrOrder', 'Tparent');
    nK = length(unique(Tparent));
    
    %% plot unconcatenated version of traces (rough, from R)
    [hfig, hax] = figureTracesI_PP_pixels( length(unique(Tparent)), 100, 60 );  %figureTracesI_PP( length(unique(Tparent)), 1 );
    hax = flipud(hax); % index 1 is bottom axis
    for Ki = 1:nK
        axes(hax(Ki)), hold on
        allR = R(Tparent==Ki,:);
        avg = mean(allR,1);
        avg = reshape(avg(:), [],4);
        for i = 1 : 4 
            plot(avg(:,i), 'Color', colors(i,:))
        end
        plot([0, size(avg,1)], [0 0], '-k')
        ylabel(sprintf('%d\n%d',Ki, sum(Tparent==Ki)))
        hax(Ki).YLabel.FontSize = 7.5;
    end
    export_fig(sprintf('clustMain_%d_mtraces.pdf', nK))
    close
    %%
    if K > K1
        T1 = load(fullfile(clusterFolder, sprintf('klust_maxClust%d.mat', K-1)));
        [table_T12] = crosstab(T1.Tparent, Tparent);
        figure; imagesc(table_T12); colormap(map)
        axis square
        ax = gca;
        ax.YTick = 1:K-1;
        ax.XTick = 1:K;
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        export_fig(sprintf('table_pair_%d_%d.png',K-1, K))
        close
    end
end


%% only remake tables
% K1 = 7;
% Kend = 16;
% for K = K1:Kend
%     load(fullfile(clusterFolder, sprintf('klust_maxClust%d.mat', K)))
%     disp(K)
%     if K > K1
%         load(fullfile(clusterFolder, sprintf('klust_maxClust%d.mat', K)))
%         T1 = load(fullfile(clusterFolder, sprintf('klust_maxClust%d.mat', K-1)));
%         [table_T12] = crosstab(T1.Tparent, Tparent);
%         figure; imagesc(table_T12); colormap(map)
%         axis square
%         ax = gca;
%         ax.YTick = 1:K-1;
%         ax.XTick = 1:K;
%         ax.XGrid = 'on';
%         ax.YGrid = 'on';
%         export_fig(sprintf('table_pair_%d_%d.png',K-1, K))
%         close
%     end
% end


%% plot a higher N dendrogram for prectical reference
maxLeavesDend = 300;
figure; [~,Tdend,outpermLeaves] = dendrogram(Zc, maxLeavesDend, 'Orientation', 'left');
axDend = gca;
axDend.YTickLabel = [];
export_fig(sprintf('dendrogram_maxLeaves%d.pdf',maxLeavesDend))


%% select a K and plot subclusters (traces)
K = 12;
load(fullfile(clusterFolder, sprintf('klust_maxClust%d.mat', K)));   %'klust', 'dendrOrder', 'Tparent'); % these are sorted as in dedndrogram

maxClustComponents = 60;
folder = fullfile(clusterFolder, sprintf('%dclusters', K));
if exist(folder) ~= 7
    mkdir(folder);
end
cd(folder)
    
Tcomponents = cluster(Zc,'maxclust',maxClustComponents);
[table_T10030] = crosstab(Tcomponents, Tparent);
figure; imagesc(table_T10030)
axis square
ax = gca;
ax.XTick = 1:K;
ax.XTickLabel = sum(table_T10030);
ax.XTickLabelRotation = 90;
export_fig('table_main_components.png')

nrows = K;
ncols = max(sum(table_T10030>0)) + 1;
[hfig, hax] = figureTracesI_PP_pixels( nrows, 100*ones(1,ncols), 60 );
hax = flipud(hax);
for Ki = 1:K
    % parent cluster
    axes(hax(Ki,1)), hold on
    allR = R(Tparent==Ki,:);
    avg = mean(allR,1);
    avg = reshape(avg(:), [],4);
    for i = 1 : 4
        plot(avg(:,i), 'Color', colors(i,:))
    end
    plot([0, size(avg,1)], [0 0], '-k')
    ylabel(sprintf('%d\n%d',Ki, sum(Tparent==Ki)))
    hax(Ki,1).YLabel.FontSize = 8.5;
    
    % children clusters
    
    compNums = find(table_T10030(:,Ki));
    for i = 1 : sum(table_T10030(:,Ki)>0) %number of component clusters
        axes(hax(Ki,i+1)), hold on
        allR = R(Tcomponents==compNums(i),:);
        avg = mean(allR,1);
        avg = reshape(avg(:), [],4);
        for icol = 1 : 4
            plot(avg(:,icol), 'Color', colors(icol,:))
        end
        plot([0, size(avg,1)], [0 0], '-k')
        ylabel(sprintf('%d',sum(Tcomponents==compNums(i))))
        hax(Ki,i+1).YLabel.FontSize = 6.5;
        drawnow
    end
    linkaxes(hax(Ki,:), 'y')
end

export_fig(sprintf('clust%2d_subclusters_responses.pdf',K))
close


%% cluster colors (12)
nK = 12;
id_colors = 1:nK;
% cmap = distinguishable_colors(nK); % will change it later
cmap = cmap(id_colors, :);

% save('klust_maxClust12.mat', 'id_colors', 'cmap', '-append');


%% set for cluster maps
K = 12;
load(fullfile(clusterFolder, sprintf('klust_maxClust%d.mat', K)));   %'klust', 'dendrOrder', 'Tparent'); % these are sorted as in dedndrogram

folderCluK = fullfile(clusterFolder, sprintf('%dclusters', K));
if exist(folderCluK, 'dir') ~= 7
    mkdir(folderCluK);
end
cd(folderCluK)


clust = remapIndexedPointsToRun_RoisMaps_oldNew(Tparent, R_100, R_86);      %sukee

% make maps (divided old/new because different dimensions)
klust = makeplot_cMaps_superKlusts_singleIteration_oldNew(clust, klust, dendrOrder);
folder = fullfile(folderCluK, 'alignedMaps');
if exist(folder, 'dir') ~= 7
    mkdir(folder)
end
cd(folder)

Col = cat(1,[1,1,1],cmap);

regGeneralData = matfile('/Users/galileo/Dropbox (HMS)/Data/WINDS_100/anatomyImages/anatomy_alignment_metadata.mat');
allStacks = regGeneralData.allStacks; %allStacks
allStacks_flyNums = cat(1,allStacks.flyNum);


figure; imagesc(1:size(cmap,1)); colormap(cmap)
ax = gca;
ax.YTick = [];
ax.XTick = 1:nK;
ax.Box = 'off';
ax.XAxis.FontSize = 25;
export_fig('cluster_colorbar.tif')


%% make cluster maps - flies 100 only
theseFlyNums = [166   169   175 ]; %  157   1701  1700   179   180]; 

for f = 1:length(theseFlyNums)
    flyNum = theseFlyNums(f);
    i_allSt = find(allStacks_flyNums == flyNum); %ef fix for fly 118
    if flyNum > 1000
        hemisph = mod(flyNum,10);
        flyNum = floor(flyNum/10);
        zs_fly = dataTable.fly==flyNum & dataTable.ipsilateral == hemisph; %ordinal numberwithin datalist
        runNames = dataTable.run(dataTable.fly==flyNum & dataTable.ipsilateral == hemisph); %actual run number, to reconstruct flyname
    else
        zs_fly = dataTable.fly==flyNum; %ordinal numberwithin datalist
        runNames = dataTable.run(dataTable.fly==flyNum);
    end  
    flyNum = theseFlyNums(f);
    disp('--- fly: -------------------------------')
    disp(flyNum)
    disp('-')
    %% sort runs from ventral to dorsal
    zetas = find(zs_fly);
    sliceNums = regGeneralData.alignedSliceNumbers(1,zetas);
    [sliceNums,b] = sort(sliceNums);
    zetas = zetas(b);
        
    %% load transformation data and make maps
    for iz = 1:sum(zs_fly)
        z = zetas(iz);          % change for flies 86
        sliceNum = sliceNums(iz);
        % load specific run's transformation
        regData = matfile(fullfile(runfolders{z}, ['stackRegistration_' basenames{z} '.mat']));
        
        al2fly119_mapAll = [];
        for Ki = 1:length(unique(Tparent))
            K2use = dendrOrder(Ki);
            mapClean = klust(K2use).mapsZold(:,:,z);    % change for flies 86
            mapClean(mapClean == nK+1) = 0;
            mapClean(mapClean>=1) = 1; %binary so far                   % 69 x 86

            %convert from downsampled to fullsize
            mapClean = cat(1, mapClean, zeros(1,size(mapClean, 2)));    % 70 x 86
            mapClean = cat(2, mapClean, zeros(size(mapClean, 1), 4));   % 70 x 90
            mapClean = imresize(mapClean, [100 128]);                   % 100 x 128
            %re-crop
            mapClean(end,:) = [];
            mapClean(:,end-5:end) = [];                                 % 99 x 122
            mapClean(mapClean<0.3) = 0;
            mapClean(mapClean>=0.1) = 1;                        %note 1 and not K2use!
            
            
            %transformations here were NOT!!!! calculated over already downsampled data
            mapClean = imrotate(mapClean, regGeneralData.totAng(1,z)); %preliminary transformation
            recovered_mapClean = imwarp(mapClean,regData.tform,'OutputView',regData.Roriginal,'FillValues',0); %indexed image
            al2fly119_mapClean = imwarp(recovered_mapClean,allStacks(i_allSt).tform,'OutputView',allStacks(i_allSt).Roriginal,'FillValues',0); %ok
            al2fly119_mapClean(al2fly119_mapClean<=0.5) = 0;
            al2fly119_mapClean(al2fly119_mapClean>=0.1) = Ki; %dendrogram is colored in order from 1st to nK-esim  in colormap. Clusters are colored accordingy.
            if isempty(al2fly119_mapAll)
                al2fly119_mapAll = al2fly119_mapClean;
            else
                al2fly119_mapAll(al2fly119_mapClean == Ki) = Ki; %dendrogram is colored in order from 1st to nK-esim  in colormap. Clusters are colored accordingy.
            end
        end

        AllClusterMaps(:,:,z) = al2fly119_mapAll;

        AlphaMatrix = zeros(size(al2fly119_mapAll(:,:,1)));
        AlphaMatrix(al2fly119_mapAll(:,:,1) > 0) = 1;

        al2fly119_mapAll = ind2rgb(uint8(al2fly119_mapAll), Col);
        savename = sprintf('alignALLklustsMap_fly%02d_ordZ%02d_%s.png', f, z, basenames{z});
        imwrite(al2fly119_mapAll, savename, 'Alpha', double(AlphaMatrix))
    end
end


% make cluster maps - flies 86 only
theseFlyNums = [ 157   1701  1700   179   180]; 

for f = 1:length(theseFlyNums)
    flyNum = theseFlyNums(f);
    if flyNum > 1000
        hemisph = mod(flyNum,10);
        flyNum = floor(flyNum/10);
        zs_fly = dataTable.fly==flyNum & dataTable.ipsilateral == hemisph;      %logical indexing within datalist
        runNames = dataTable.run(dataTable.fly==flyNum & dataTable.ipsilateral == hemisph);     %actual run number, to reconstruct flyname
    else
        zs_fly = dataTable.fly==flyNum;     %logical indexing within datalist
        runNames = dataTable.run(dataTable.fly==flyNum);
    end 
    i_allSt = find(allStacks_flyNums == flyNum); %ef fix for fly 118
    
    %% sort runs from ventral to dorsal
    zetas = find(zs_fly); 
    sliceNums = dataTable.stacksZ(zs_fly);   
    [sliceNums,b] = sort(sliceNums);        % stack slices paired to any functional run, (sorted ventral to dorsal)
    zetas = zetas(b);                       % ordinal numbers relative to datalist, (sorted ventral to dorsal)
    
    %% load transformation data and make maps
    for iz = 1:sum(zs_fly)          %iz is local counter within fly
        z = zetas(iz);              % relative to datalist
        z_86 = z - NaZ_100;         % now this is relative to the R_86 subset
        sliceNum = sliceNums(iz);   % still relative to its own stack
        % load specific run's transformation
        regData = matfile(fullfile(runfolders{z}, ['stackRegistration_' basenames{z} '.mat']));
        
        
        al2fly119_mapAll = [];
        for Ki = 1:length(unique(Tparent))
            K2use = dendrOrder(Ki);
            mapClean = klust(K2use).mapsZnew(:,:,z_86);    % change for flies 86
            mapClean(mapClean == nK+1) = 0;
            mapClean(mapClean>=1) = 1; %binary so far

            %convert from downsampled to fullsize
            mapClean = cat(1, mapClean, zeros(1,size(mapClean, 2)));    
            mapClean = cat(2, mapClean, zeros(size(mapClean, 1), 4));  
            mapClean = imresize(mapClean, [86 128]);                   % 86 x 128
            %re-crop
            mapClean(end,:) = [];
            mapClean(:,end-5:end) = [];                             
            mapClean(mapClean<0.3) = 0;
            mapClean(mapClean>=0.1) = 1;                        %note 1 and not K2use!
            
            
            %transformations here were calculated over already downsampled data
            mapClean = imrotate(mapClean, regGeneralData.totAng(1,z)); %preliminary transformation
            recovered_mapClean = imwarp(mapClean,regData.tform,'OutputView',regData.Roriginal,'FillValues',0); %indexed image
            al2fly119_mapClean = imwarp(recovered_mapClean,allStacks(i_allSt).tform,'OutputView',allStacks(i_allSt).Roriginal,'FillValues',0);
            al2fly119_mapClean(al2fly119_mapClean<=0.5) = 0;
            al2fly119_mapClean(al2fly119_mapClean>=0.1) = Ki; %dendrogram is colored in order from 1st to nK-esim  in colormap. Clusters are colored accordingy.
            if isempty(al2fly119_mapAll)
                al2fly119_mapAll = al2fly119_mapClean;
            else
                al2fly119_mapAll(al2fly119_mapClean == Ki) = Ki; %dendrogram is colored in order from 1st to nK-esim  in colormap. Clusters are colored accordingy.
            end
        end

        AllClusterMaps(:,:,z) = al2fly119_mapAll;

        AlphaMatrix = zeros(size(al2fly119_mapAll(:,:,1)));
        AlphaMatrix(al2fly119_mapAll(:,:,1) > 0) = 1;

        al2fly119_mapAll = ind2rgb(uint8(al2fly119_mapAll), Col);
        savename = sprintf('alignALLklustsMap_fly%02d_ordZ%02d_%s.png', f+3, z, basenames{z});  % change for flies 86
        imwrite(al2fly119_mapAll, savename, 'Alpha', double(AlphaMatrix))
    end
end



%% barplots - to be adapted yet
% clust30 = remapIndexedPointsToRun_RoisMaps(Tparent, pixel2keep, sizeMap, pxKeep);
% clust100 = remapIndexedPointsToRun_RoisMaps(Tcomponents, pixel2keep, sizeMap, pxKeep);
% sortFliesGenoType = 0;
% [~, flyClusterCountDivZs100] = barcluster_panSpec(Tcomponents, datalist, flyfoldersUnique, clust100.ROIs, sortFliesGenoType);
% [~, flyClusterCountDivZs30] = barcluster_panSpec(Tparent, datalist, flyfoldersUnique, clust30.ROIs, sortFliesGenoType);


% R hemisph
hemisph = 1;
[~, flyClusterCountDivZsR] = barcluster_winds(Tparent, dataTable, hemisph, clust);


% L hemisph
hemisph = 0;
[~, flyClusterCountDivZsL] = barcluster_winds(Tparent, dataTable, hemisph, clust);

%%
[hfig, hax] = figureTracesI_PP_pixels( nK, 300*ones(1,2), 60 );
hax = flipud(hax);
for K = 1:nK
    axes(hax(K,1))
    bar(1:NaF, flyClusterCountDivZsR(:,K))
    ylabel(K)
    axes(hax(K,2))
    bar(1:NaF, flyClusterCountDivZsL(:,K))
    for i = 1:2
        hax(K,i).XTick = 1:NaF;
        hax(K,i).FontSize = 7;
        hax(K,i).XTickLabel = flyNumUnique;
        hax(K,i).XTickLabelRotation = 90;
    end
end
axes(hax(nK,1))
title('Right hemisphere (ipsi)')
axes(hax(nK,2))
title('Left hemisphere (contra)')

export_fig(fullfile(folderCluK, 'barClusterDistrib_RvsL.eps'))


