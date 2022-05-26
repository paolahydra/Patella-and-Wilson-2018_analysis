% load up
load_panSpec_downsampled;
load('R_matrix_Downsampled_smallWindow.mat', 'pixel2keep', 'iKeep', 'pxKeep')
clusterFolder = '/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix';
treeFileName = 'output_ward_movAVG_ZSCpix_linkage.mat';
load(fullfile(clusterFolder, treeFileName), 'Zc');

sizeMap = [60, 86];

cd(clusterFolder)
regGeneralData = matfile('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/anatomyImages/anatomy_alignment_metadata.mat');
for z = 1:length(datalist)
    [runfolders{z}, ~] = fileparts(datalist{z});
end
load('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/anatomyImages/JON_cropMask.mat')

%%
for K = 3:50
    nK = K;
    clustName = sprintf('klust_maxClust%d.mat', K);
    load(clustName)
    % make/save registered maps
    clust = remapIndexedPointsToRun_RoisMaps(Tparent, pixel2keep, sizeMap);
    % make maps (divided old/new because different dimensions)
    klust = makeplot_cMaps_superKlusts_singleIteration(clust, klust, dendrOrder);
    
    clear singleKmaps
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
            
            if Ki == 1 && K == 2
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
        end
    end
    saveclustName = sprintf('singleKmaps_maxClust%d.mat', K);
    save(saveclustName, 'singleKmaps')
    
    if K == 2
        saveclustName = sprintf('singleKmaps_maxClust1.mat');
        save(saveclustName, 'singleKmaps_px2keep')
    end
end


% try computing b2b_idx and compact_idx

% use cropmask! even if it's a maxproj
cropRows = sum(cropMask,2) > 0;
cropCols = sum(cropMask,1) > 1;
% do also baseline for K == 1
%%
for K = 2:50
    saveclustName = sprintf('singleKmaps_maxClust%d.mat', K);
    load(saveclustName, 'singleKmaps')
    clear maxProjClust
    % by fly
    for f_i = 1:6   % length(flyNumUnique)          %pan only
        f = flyNumUnique(f_i);
        % take all corresponding zetas
        zetas = dataTable.fly == f;
        singleKmaps_fly = singleKmaps(:,:,zetas,:);
        
        % max proj by clust
        for Ki = 1:K
            % do max proj
            maxProjClust(:,:,Ki,f_i) = max(singleKmaps_fly(:,:,:,Ki),[],3);
        end
    end
    maxProjClust = double(maxProjClust);
    for Ki = 1:K
        % do max proj
        b2b_Clust(:,:,Ki) = mean(maxProjClust(:,:,Ki,:), 4);
    end
    b2b_map(:,:,K) = max(b2b_Clust, [], 3);
    % crop and average
    maxProjClust = max(maxProjClust, [], 4);
    maxProjClust = maxProjClust(cropRows, cropCols,:);
    b2b_ClustMask = b2b_Clust(cropRows, cropCols,:);
    allMeans = [];
    for Ki = 1:K
        cmax = logical(maxProjClust(:,:,Ki));
        cmeans = b2b_ClustMask(:,:,Ki);
        allMeans = cat(1, allMeans, cmeans(cmax));
    end
    b2b_idx(K) = mean(allMeans);
end

%
save('allB2B.mat', 'b2b_idx', 'b2b_map')

