%% WINDS (flipped). EXPLORE AVERAGE LINKAGE to search for uncorrelated pixels

cd(fullfile(Folder2Save, 'rawDataLinkage'))
load('output_average_Rzsp_linkage_windsFLIPPED.mat', 'Zc');
figure; dendrogram(Zc, 300, 'Orientation', 'left');
distance = 12.2; 
T = cluster(Zc,'cutoff', distance, 'criterion', 'distance'); 
length(unique(T))

pxKeep = true(size(T));   % (first time only)
iteratio = 1;  % (first time only)

%% first and following iterations
distance = 12.2; 
T = cluster(Zc,'cutoff', distance, 'criterion', 'distance');  
selClusts = unique(T(pxKeep));
nK = length(selClusts);
percPixelsRetained = sum(pxKeep) / size(R,1) * 100;

hhist = figure; h = histogram(categorical(T(pxKeep)));
countSelClusters = h.Values;
ax = gca;
title(sprintf('Iter %d: Retained %3.1f%% pixels, divided into %d clusters by threshold distance %.3f', iteratio, percPixelsRetained, nK, distance))

threshcount = 500;
selClusts = selClusts(h.Values>threshcount);
nK = length(selClusts)
countSelClusters = h.Values(h.Values>threshcount);

fprintf('Iter %d: Retained %3.1f%% pixels, divided into %d clusters,\nof which %d with more than %d pixels, by threshold distance %.3f \n', iteratio, percPixelsRetained, length(unique(T(pxKeep))), nK, threshcount, distance)

%% Part 1: plot
[hfig, hax] = figureTracesI_PP( length(selClusts), 1 );
hax = flipud(hax);
clear centroids
for ik = 1:length(selClusts)
    k = selClusts(ik);
    centroids(ik,:) = mean(R(T==k,:));
    axes(hax(ik)), hold on
    upperbound = prctile(R(T==k,:), 75, 1);
    lowerbound = prctile(R(T==k,:), 25, 1);
    plot(upperbound, '-k'), hold on
    plot(lowerbound, '-k')
    plot(centroids(ik,:), '-r')
    axis tight
    hax(ik).YLabel.String = sprintf('k%d\n#%d', k, countSelClusters(ik));
    hax(ik).YLabel.FontSize = 6;
    drawnow
end
title(sprintf('Iter %d: Retained %3.1f%% pixels, divided into %d clusters,\nof which %d with more than %d pixels, by threshold distance %.3f \n', iteratio, percPixelsRetained, length(unique(T(pxKeep))), nK, threshcount, distance))
hfig.Position = [ 0    0.0044    0.3257    0.8900];
export_fig(hfig, sprintf('CutNoise_%d_Thr%.3f.pdf', iteratio, distance))
export_fig(hhist, sprintf('CutNoise_%d_Thr%.3f_hist.png', iteratio, distance))

%% Part 1: make a decision and reiterate (run this just once)
save(sprintf('retainedPixels_AtIteratio%d.mat', iteratio), 'pxKeep', 'iteratio')

keepClusts = selClusts; %[1 293 298 859 864 928 936 949];        %selClusts;      %change here

newPxKeep = false(size(T));
for i = 1:length(keepClusts)
    newPxKeep = newPxKeep | T==keepClusts(i);
end
pxKeep = newPxKeep;
iteratio = iteratio+1;
close(hfig)


%% show which pixels have been excluded
T(pxKeep) = 1;
T(~pxKeep) = 2;


clear klust allElements firstLeaf dendrOrder
for Ki = 1:length(unique(T))
    tempKs = find(T==Ki);
    klust(Ki).k = (tempKs(:))';
end
nK = length(unique(T));
colors = distinguishable_colors(nK);
runfolders = sets.analysislist;

Tparent = T;
dendrOrder = 1:nK;
clust = remapIndexedPointsToRun_RoisMaps_oldNew(Tparent, R_100, R_86);      %sukee

% make maps (divided old/new because different dimensions)
klust = makeplot_cMaps_superKlusts_singleIteration_oldNew(clust, klust, dendrOrder);
folder = '/Users/galileo/Dropbox (HMS)/Data/WINDS_100/Downsampled/alignedMaps';
if exist(folder, 'dir') ~= 7
    mkdir(folder)
end
cd(folder)

Col = cat(1,[1,1,1],colors);

regGeneralData = matfile('/Users/galileo/Dropbox (HMS)/Data/WINDS_100/anatomyImages/anatomy_alignment_metadata.mat');
allStacks = regGeneralData.allStacks; %allStacks
allStacks_flyNums = cat(1,allStacks.flyNum);
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


