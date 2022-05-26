% this is a patch to create outlines of the entire bundle for select old panJON
% flies


%% 

% run first part of:
edit figure_SSteps_JONPAN_oldnew;

flyfoldersUnique_old = flyfoldersUnique_old([1,2,4,5]);
flyfoldersUnique = flyfoldersUnique_old;
d = dir('/Users/galileo/Dropbox (HMS)/Data/JON_old/oldPanAllSIGPXLS');
d(2) = [];
d(1) = [];

%%
for iZ = 1:length(d)
    
    basename = d(iZ).name(10:21);
    
    ZZ = strfind(regGeneralData.runfolders, basename);
    zs_fly = zeros(size(ZZ));
    for i = 1:length(zs_fly)
        if ~isempty(ZZ{i})
            zs_fly(i) = 1; % indices into datalist relative to fly f
        end
    end
    z = find(zs_fly);
    
    flyfolder = regGeneralData.flyfolders(1, z);
    runfolder = regGeneralData.runfolders(1, z);
    
    
    % load fly's stack2fly150 transformation
    name = cell2mat(fullfile(flyfolder, ['stack2fly150sStack_' basename(1:6) '.mat']) );
    stack2fly150_info = matfile(name);
    
    name = cell2mat(fullfile(runfolder, ['stackRegistration_' basename '.mat']));
    regData = matfile(name);
    if ~isfield(regData, 'Roriginal')
        regData.Properties.Writable = true;
        regData.Roriginal = imref2d(size(regData.target));
    end
    
    %% load signPxl map
    name = fullfile('/Users/galileo/Dropbox (HMS)/Data/JON_old/oldPanAllSIGPXLS', d(iZ).name );
    tifSP = readTiffPP(name,1);
    tifSP(end,:) = [];
    tifSP(:,end-5:end) = [];
    mapClean = imresize(tifSP, [55 86]);
    mapClean(mapClean<0.3) = 0;
    mapClean(mapClean>=0.1) = 1;
    
    mapClean = imrotate(mapClean, regGeneralData.totAng(1,z)); %preliminary transformation
    recovered_mapClean = imwarp(mapClean,regData.tform,'OutputView',regData.Roriginal,'FillValues',0); %indexed image
    al2fly150_mapClean = imwarp(recovered_mapClean,stack2fly150_info.tform,'OutputView',stack2fly150_info.Roriginal,'FillValues',0);

    al2fly150_mapClean(al2fly150_mapClean<=0.5) = 0;
    al2fly150_mapClean(al2fly150_mapClean>=0.1) = 1; %dendrogram is colored in order from 1st to nK-esim  in colormap. Clusters are colored accordingy.
        
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
    savename = sprintf('A_alignALLklustsMap_outlineALLSTIM_%s.eps', basename);
    export_fig(savename)
    close
end
