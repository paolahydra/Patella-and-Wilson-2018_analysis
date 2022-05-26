% % 2 : Can I retrieve the right data?  -- YES
folder = '/Users/galileo/Dropbox (HMS)/Data/JON_pan_old&new/DownsampledFinal';
cd(folder)
data = '/Users/galileo/Dropbox (HMS)/Data/JON_pan_old&new/DownsampledFinal/R_matrix_old&new_panJONs_steps.mat'; %there is no baseline here: from t=0 to last point available.
load(data); % this load R
% % 1 : Had I downsampled it all spatially?
% % YES
% % template code for clustering:
% edit clustering_on_rawData.m

nK = 2;
load('klust_maxClust2.mat')
% figure; hold on
% for k = 1:nK
%     plot(mean(R(Tparent==k,:)));
% end
% 
% export_fig('clusters_1_2_PullPush.eps')

%% una tantum - clustering - done
% %% zscore and operate on this for now ( I need orchestra, or my new computer...)
% 
% % alternatively, try using a moving mean to further smooth the signal
% recR2 = (movmean(R', 9))';
% recR2 = zscore(recR2, 0, 2);
% save('R_matrix_Downsampled_movAvg9_zscoredByPixel.mat', 'recR2')
% 
% % calculating on orchestra (180107) 
% % bsub -q priority -W 44:00 -R "rusage[mem=60000]" -n 1 'matlab -nojvm -singleCompThread -nodisplay -r orchestra_ward_output_steps_panOLDNEW > output_ward_Rzsp_linkage_steps_panOLDNEW'
% 
% %% retreive linkage calulated on orchestr and decide how many clusters
% treeFileName = 'output_ward_Rzsp_linkage_steps_panOLDNEW.mat';
% load(fullfile(folder, treeFileName), 'Zc');
% 
% 
% maxClust = 2;
% % maxClustComponents = 150;
% maxLeavesDend = maxClust;
% distanceMetrics = 'euclidean';
% linkageMethod = 'ward'; 
% figure; [~,Tdend,outpermLeaves] = dendrogram(Zc, maxLeavesDend, 'Orientation', 'left');
% 
% %
% Tparent = cluster(Zc,'maxclust',maxClust);
% 
% 
% for Ki = 1:length(unique(Tparent))
%     tempKs = find(Tparent==Ki);
%     klust(Ki).k = (tempKs(:))';
% end
% for Ki = 1:length(unique(Tparent)), allElements(Ki, 1:length(klust(Ki).k)) = klust(Ki).k; end
% 
% countLeaves = 1;
% for Ki = 1:length(unique(Tparent))   %while ~isempty(leavesorder)
%     lengthLeave = 0;
%     leavesorder = find(Tdend == outpermLeaves(countLeaves));
%     countLeaves = countLeaves+1;
%     lengthLeave = lengthLeave + length(leavesorder);
%     firstLeaf(Ki) = leavesorder(1);
%     [dendrOrder(Ki), ~] = find(allElements==firstLeaf(Ki));
%     while sum(Tparent == dendrOrder(Ki)) > lengthLeave
%         leavesorder = find(Tdend == outpermLeaves(countLeaves));
%         countLeaves = countLeaves+1;
%         lengthLeave = lengthLeave + length(leavesorder);
%     end
% end
% clear leavesorder countLeaves leavesorder
% 
% % Also rename T once and for all according to dendrogram
% T = zeros(size(Tparent));
% for Ki = 1:length(unique(Tparent))
%     T(Tparent==dendrOrder(Ki)) = Ki;
%     tempKs = find(T==Ki);
%     klustSorted(Ki).k = (tempKs(:))';
% end
% Tparent = T;
% klust = klustSorted;
% dendrOrder = 1:length(unique(Tparent));
% clear klustSorted T
% clear outpermLeaves %to avoid errors
% 
% save(sprintf('klust_maxClust%d.mat', maxClust), 'klust', 'dendrOrder', 'Tparent');
% nK = length(unique(Tparent));
% 
% % display mean dff
% for k = 1:nK
%     figure;
%     plot(mean(R(Tparent==k,:)));
%     title(sprintf('cluster %d, count %d', k, sum(Tparent==k)))
% end



% ORDER IN R and Tparent is old , new

%% recovered material that works like charm:
% of use: 
% edit select_tSNEpoints_andShowMapsAndSignals.m
[datalist, basenames] = loadDatalist_oldNew_stepsOnly;
% retrieve runs/maps info
[clust] = remapIndexedPointsToRun_RoisMaps_oldNew(Tparent, Rold, Rnew);


load('/Users/galileo/Dropbox (HMS)/Data/cmap19_spectral.mat', 'cmap')

% cmap10 = cmap([2 1 19 17 14 12 10 8 6 5],:);
% colors = cmap10([2,1],:);
colors = [0.55, 0.55, 0.55; 0.1200,    0.1200,    0.1200]; %brighter gray

Col = cat(1,[1,1,1],colors);
clear cmap cmap10
clust(1).cMAPs = makeColoredMaps(clust(1).MAPs, colors);
clust(2).cMAPs = makeColoredMaps(clust(2).MAPs, colors);

for z = 1:length(datalist)
    [runfolders{z}, ~] = fileparts(datalist{z});
end

%% needed
[flyfolders, ~] = extractRunfoldersBasenames(datalist);
flyfoldersUnique = unique(flyfolders, 'stable');
flyfoldersUnique_old = flyfoldersUnique(1:5);
flyfoldersUnique_new = flyfoldersUnique(6:end);


%% now just map 'em out!
klust = makeplot_cMaps_superKlusts_singleIteration_oldNew(clust, klust, dendrOrder, colors);
mfolder = fullfile(folder, 'alignedMaps');
mkdir(mfolder)
cd(mfolder)

%% maps old (aligned)
regGeneralData = matfile('/Users/galileo/Dropbox (HMS)/Data/JON_old/anatomyImages/anatomy_alignment_metadata.mat');
flyfoldersUnique = flyfoldersUnique_old;
NaF = length(flyfoldersUnique_old);

Zclust = 1;
AllSigPxlsMaps = [];


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
            K2use = Ki;
            mapClean = klust(K2use).mapsZold(:,:,z);
            mapClean(mapClean == nK+1) = 0;
            mapClean(mapClean>=1) = 1; %binary so far

            %convert from downsampled to fullsize --  NOTE: OLD dataset had been
            %(RE-)calculated over the downsampled maps directly.
            
%             mapClean = cat(1, mapClean, zeros(1,size(mapClean, 2))); %61 x 86
%             mapClean = cat(2, mapClean, zeros(size(mapClean, 1), 4)); %61 x 90
%             mapClean = imresize(mapClean, [80 128]); %should become 86 x 128. The quality sucks, but so be it
%             %re-crop
%             mapClean(end,:) = [];
%             mapClean(:,end-5:end) = [];
%             mapClean(mapClean<0.3) = 0;
%             mapClean(mapClean>=0.1) = 1;                        %note 1 and not K2use!

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
        
        
        %% repeat alignement for all sign pixels outline (if I had it)
        mapClean = reshape(Rold_chirpSteps.pixel2keep(:,z), Rold_chirpSteps.sizeMAP);
        mapClean = imresize(mapClean, [55 86]);
        mapClean = imrotate(mapClean, regGeneralData.totAng(1,z)); %preliminary transformation
        recovered_mapClean = imwarp(mapClean,regData.tform,'OutputView',regData.Roriginal,'FillValues',0); %indexed image
        if strcmp(basenames{find(zs_fly,1)}(1:6), 'fly150')
            al2fly150_mapClean = recovered_mapClean;
        else
            al2fly150_mapClean = imwarp(recovered_mapClean,stack2fly150_info.tform,'OutputView',stack2fly150_info.Roriginal,'FillValues',0);
        end
        al2fly150_mapClean(al2fly150_mapClean<=0.5) = 0;
        al2fly150_mapClean(al2fly150_mapClean>=0.1) = 1; %dendrogram is colored in order from 1st to nK-esim  in colormap. Clusters are colored accordingy.
        
        
        AllSigPxlsMaps(:,:,z) = al2fly150_mapClean;
        
        BW = imfill(al2fly150_mapClean,'holes');
        
        BW = imdilate(BW, strel('disk',4));
        BW = imclose(BW, strel('disk',4) );
        BW = imerode(BW, strel('disk',4));
        BW = imopen(BW, strel('disk',1) );
        
        BW = imclose(BW, strel('disk',24) );
        BW = imdilate(BW, strel('disk',4));
        BW = imfill(BW,'holes');
        
        se = strel('line',4,0);
        BW = imerode(BW, se);
        se = strel('line',4,90);
        BW = imerode(BW, se);
        se = strel('line',4,45);
        BW = imerode(BW, se);
        se = strel('line',4,135);
        BW = imerode(BW, se);
        BW = imdilate(BW, strel('disk',2));
        
        
        figure; imshow(ones(size(BW))); hold on
        [B,L,N] = bwboundaries(BW); 
        for k = 1:length(B)
            boundary = B{k};
            if k <= N %do not include boundaries of any holes....
                plot(boundary(:,2), boundary(:,1), 'r','LineWidth',2);
            end
        end
        savename = sprintf('alignALLklustsMap_outlineCHIRPSSTEPS_fly%02d_ordZ%02d_%s.eps', f, z, basenames{z});
        export_fig(savename)
        close
        
        %% repeat alignement for all sign pixels outline 
        mapClean = reshape(Rold.pixel2keep(:,z), Rold.sizeMAP);
        mapClean = imrotate(mapClean, regGeneralData.totAng(1,z)); %preliminary transformation
        recovered_mapClean = imwarp(mapClean,regData.tform,'OutputView',regData.Roriginal,'FillValues',0); %indexed image
        if strcmp(basenames{find(zs_fly,1)}(1:6), 'fly150')
            al2fly150_mapClean = recovered_mapClean;
        else
            al2fly150_mapClean = imwarp(recovered_mapClean,stack2fly150_info.tform,'OutputView',stack2fly150_info.Roriginal,'FillValues',0);
        end
        al2fly150_mapClean(al2fly150_mapClean<=0.5) = 0;
        al2fly150_mapClean(al2fly150_mapClean>=0.1) = 1; %dendrogram is colored in order from 1st to nK-esim  in colormap. Clusters are colored accordingy.
        
        
        AllSigPxlsMaps(:,:,z) = al2fly150_mapClean;
        
        BW = imfill(al2fly150_mapClean,'holes');
        
        BW = imdilate(BW, strel('disk',3));
        BW = imclose(BW, strel('disk',3) );
        BW = imerode(BW, strel('disk',3));
        BW = imopen(BW, strel('disk',1) );
        
        BW = imclose(BW, strel('disk',20) );
        BW = imdilate(BW, strel('disk',4));
        BW = imfill(BW,'holes');
        
        se = strel('line',4,0);
        BW = imerode(BW, se);
        se = strel('line',4,90);
        BW = imerode(BW, se);
        se = strel('line',4,45);
        BW = imerode(BW, se);
        se = strel('line',4,135);
        BW = imerode(BW, se);
        BW = imdilate(BW, strel('disk',2));
        
        
        figure; imshow(ones(size(BW))); hold on
        [B,L,N] = bwboundaries(BW); 
        for k = 1:length(B)
            boundary = B{k};
            if k <= N %do not include boundaries of any holes....
                plot(boundary(:,2), boundary(:,1), 'r','LineWidth',2);
            end
        end
        savename = sprintf('alignALLklustsMap_outline_fly%02d_ordZ%02d_%s.eps', f, z, basenames{z});
        export_fig(savename)
        close
    end
end


save(sprintf('alignedMatrix_Images_maxClust%d_oldOnly.mat',nK), 'AllClusterMaps', 'Col')



%% maps new (aligned)
regGeneralData = matfile('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/anatomyImages/anatomy_alignment_metadata.mat'); % includes spec
flyfoldersUnique = flyfoldersUnique_new;
NaF = length(flyfoldersUnique_new);
% better reset zetas
datalist(1:22) = [];
basenames(1:22) = [];
runfolders(1:22) = [];

Zclust = 2;

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
            K2use = Ki;
            mapClean = klust(K2use).mapsZnew(:,:,z);
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
        savename = sprintf('alignALLklustsMap_fly%02d_ordZ%02d_%s.png', f+5, z, basenames{z});
        imwrite(al2fly150_mapAll, savename, 'Alpha', double(AlphaMatrix))
    
        
        
        %% ADAPT FROM OLD
        %% repeat alignement for all sign pixels outline
        mapClean = reshape(Rnew.pixel2keep(:,z), Rnew.sizeMAP);
        %convert from downsampled to fullsize
        mapClean = cat(1, mapClean, zeros(1,size(mapClean, 2))); %61 x 86
        mapClean = cat(2, mapClean, zeros(size(mapClean, 1), 4)); %61 x 90
        mapClean = imresize(mapClean, [86 128]); %should become 86 x 128. The quality sucks, but so be it
        %re-crop
        mapClean(end,:) = [];
        mapClean(:,end-5:end) = [];
        mapClean(mapClean<0.3) = 0;
        mapClean(mapClean>=0.1) = 1;
        
        mapClean = imrotate(mapClean, regGeneralData.totAng(1,z)); %preliminary transformation
        recovered_mapClean = imwarp(mapClean,regData.tform,'OutputView',regData.Roriginal,'FillValues',0); %indexed image
        if strcmp(basenames{find(zs_fly,1)}(1:6), 'fly150')
            al2fly150_mapClean = recovered_mapClean;
        else
            al2fly150_mapClean = imwarp(recovered_mapClean,stack2fly150_info.tform,'OutputView',stack2fly150_info.Roriginal,'FillValues',0);
        end
        al2fly150_mapClean(al2fly150_mapClean<=0.5) = 0;
        al2fly150_mapClean(al2fly150_mapClean>=0.1) = 1; %dendrogram is colored in order from 1st to nK-esim  in colormap. Clusters are colored accordingy.
        
        
        AllSigPxlsMaps(:,:,z) = al2fly150_mapClean;
        
        BW = imfill(al2fly150_mapClean,'holes');
        
        BW = imdilate(BW, strel('disk',3));
        BW = imclose(BW, strel('disk',3) );
        BW = imerode(BW, strel('disk',3));
        BW = imopen(BW, strel('disk',1) );
        
        BW = imclose(BW, strel('disk',20) );
        BW = imdilate(BW, strel('disk',4));
        BW = imfill(BW,'holes');
        
        se = strel('line',4,0);
        BW = imerode(BW, se);
        se = strel('line',4,90);
        BW = imerode(BW, se);
        se = strel('line',4,45);
        BW = imerode(BW, se);
        se = strel('line',4,135);
        BW = imerode(BW, se);
        BW = imdilate(BW, strel('disk',2));
        
        
        figure; imshow(ones(size(BW))); hold on
        [B,L,N] = bwboundaries(BW); 
        for k = 1:length(B)
            boundary = B{k};
            if k <= N %do not include boundaries of any holes....
                plot(boundary(:,2), boundary(:,1), 'r','LineWidth',2);
            end
        end
        savename = sprintf('alignALLklustsMap_outline_fly%02d_ordZ%02d_%s.eps', f+5, z, basenames{z});
        export_fig(savename)
        close
        
        %% retrieve all-sigPixels maps for these new-pan flies
        load('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/R_matrix_Downsampled_smallWindow.mat', 'pixel2keep');
        pixel2keep(:,32:end) = [];
        
        mapClean = reshape(pixel2keep(:,z), Rnew.sizeMAP);
        %convert from downsampled to fullsize
        mapClean = cat(1, mapClean, zeros(1,size(mapClean, 2))); %61 x 86
        mapClean = cat(2, mapClean, zeros(size(mapClean, 1), 4)); %61 x 90
        mapClean = imresize(mapClean, [86 128]); %should become 86 x 128. The quality sucks, but so be it
        %re-crop
        mapClean(end,:) = [];
        mapClean(:,end-5:end) = [];
        mapClean(mapClean<0.3) = 0;
        mapClean(mapClean>=0.1) = 1;
        
        mapClean = imrotate(mapClean, regGeneralData.totAng(1,z)); %preliminary transformation
        recovered_mapClean = imwarp(mapClean,regData.tform,'OutputView',regData.Roriginal,'FillValues',0); %indexed image
        if strcmp(basenames{find(zs_fly,1)}(1:6), 'fly150')
            al2fly150_mapClean = recovered_mapClean;
        else
            al2fly150_mapClean = imwarp(recovered_mapClean,stack2fly150_info.tform,'OutputView',stack2fly150_info.Roriginal,'FillValues',0);
        end
        al2fly150_mapClean(al2fly150_mapClean<=0.5) = 0;
        al2fly150_mapClean(al2fly150_mapClean>=0.1) = 1; %dendrogram is colored in order from 1st to nK-esim  in colormap. Clusters are colored accordingy.
        
        
        AllSigPxlsMaps(:,:,z) = al2fly150_mapClean;
        
        BW = imfill(al2fly150_mapClean,'holes');
        
        BW = imdilate(BW, strel('disk',3));
        BW = imclose(BW, strel('disk',3) );
        BW = imerode(BW, strel('disk',3));
        BW = imopen(BW, strel('disk',1) );
        
        BW = imclose(BW, strel('disk',20) );
        BW = imdilate(BW, strel('disk',4));
        BW = imfill(BW,'holes');
        
        se = strel('line',4,0);
        BW = imerode(BW, se);
        se = strel('line',4,90);
        BW = imerode(BW, se);
        se = strel('line',4,45);
        BW = imerode(BW, se);
        se = strel('line',4,135);
        BW = imerode(BW, se);
        BW = imdilate(BW, strel('disk',2));
        
        
        figure; imshow(ones(size(BW))); hold on
        [B,L,N] = bwboundaries(BW); 
        for k = 1:length(B)
            boundary = B{k};
            if k <= N %do not include boundaries of any holes....
                plot(boundary(:,2), boundary(:,1), 'r','LineWidth',2);
            end
        end
        savename = sprintf('alignALLklustsMap_outlineALLST_fly%02d_ordZ%02d_%s.eps', f+5, z, basenames{z});
        export_fig(savename)
        close
    
    end
end

%%
cd ..
save(sprintf('alignedMatrix_Images_maxClust%d_newOnly.mat',nK), 'AllClusterMaps', 'Col')










