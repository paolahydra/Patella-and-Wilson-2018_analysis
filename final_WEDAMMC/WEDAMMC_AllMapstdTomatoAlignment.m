load_wedAmmc_downsampled;
regGeneralData = matfile('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/anatomyImages/anatomy_alignment_metadata.mat');
for z = 1:length(datalist)
    [runfolders{z}, ~] = fileparts(datalist{z});
end

folderTD = '/Users/galileo/Dropbox (HMS)/Data/TDTstacks/';
load('/Users/galileo/Dropbox (HMS)/Data/TDTstacks/anatomy_alignment_metadata.mat') %allStacks
saveName = '/Users/galileo/Dropbox (HMS)/Data/TDTstacks/anatomy_alignment_metadata.mat'; %allStacks

allStacks_flyNums = cat(1,allStacks.flyNum);

dataTable.zOffsetAlStack(dataTable.fly==118) = 4;   % this takes fly 126's stack
dataTable.zOffsetAlStack(dataTable.fly==120) = 1;   % this goes by 4um and gets its own alignment (runs 3-1-2: mZ:18-20-21)


%% tdtomato average masks
% you have loaded dataTable and Tnew (Tnew: all variables are nan, but variables list may need to be updated.)
%these next changes are not going to be saed (for now)
Tnew.fly = 121;
Tnew.tdTomato = 1;
Tnew.zOffsetAlStack = 5;
dataTable = cat(1, dataTable, Tnew);

Tnew.fly = 138;
Tnew.tdTomato = 1;
Tnew.zOffsetAlStack = 4;
dataTable = cat(1, dataTable, Tnew);

Tnew.fly = 149;
Tnew.tdTomato = 1;
Tnew.zOffsetAlStack = -1;
dataTable = cat(1, dataTable, Tnew);
Tnew.fly = 152;
Tnew.tdTomato = 1;
Tnew.zOffsetAlStack = 1;
dataTable = cat(1, dataTable, Tnew);


tdFlies = [119 122 126 121 138];
for f = 1:length(tdFlies)
    maskst = allStacks(allStacks_flyNums==tdFlies(f)).masked_alignedstackR;
    zof = unique(dataTable.zOffsetAlStack(dataTable.fly == tdFlies(f)));
    maxiTDTMaskStack(:,:, 1+zof : size(maskst,3)+zof, f) = maskst;
end
meanTDTmasks = mean(maxiTDTMaskStack, 4);
meanTDTmasksN = bsxfun(@minus, 1, meanTDTmasks);
% figure; imagesc(meanTDTmasks(:,:,17)); axis image; colormap(gray)

for Z = 1:size(meanTDTmasks,3)
        savename = fullfile('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/Downsampled/linkedIMAGES/alignedTDTMASKS', sprintf('meanTDmasks_MAIN_%02d.png', Z));
        % imwrite(zeros(size(meanTDTmasks(:,:,Z))), savename, 'Alpha', double(meanTDTmasksN(:,:,Z)))      %transparentAMMC_blackREST
        imwrite(zeros(size(meanTDTmasks(:,:,Z))), savename, 'Alpha', double(meanTDTmasks(:,:,Z)))      %blackAMMC_transparentREST
end

%% maxProjection (except for 1s) of maxiTDT
maxProjbut1_TDTmasks = sum(maxiTDTMaskStack, 4);
maxProjbut1_TDTmasks(maxProjbut1_TDTmasks==1) = 0;
maxProjbut1_TDTmasks(maxProjbut1_TDTmasks>=2) = 1;
MIJ.createImage(maxProjbut1_TDTmasks);

strelSizeDilate = 18;
    se = strel('disk',strelSizeDilate,0);
strelSizeErode = 16;
    seErode = strel('disk',strelSizeErode,0);
strelSizeClose = 14;
    seClose = strel('disk',strelSizeClose,0);

% figure; imshow(maxProjbut1_TDTmasks(:,:,20)); title('mpbut1')

for i = 1:size(maxProjbut1_TDTmasks,3)
    tp(:,:,i) = imdilate(maxProjbut1_TDTmasks(:,:,i), se);
    tp(:,:,i) = imclose(tp(:,:,i), seClose);
    tp(:,:,i) = imerode(tp(:,:,i), seErode);
end
figure; imshow(tp(:,:,20)); title(sprintf('dilate close erode %d %d %d',strelSizeDilate, strelSizeErode, strelSizeClose))

% save tp masks info %-added later
alignedAMMCmasks = logical(tp);
save('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/datalist_WEDAMMC_piezo.mat', 'alignedAMMCmasks', '-append')

%% make boundaries from maxProjBut1 and save
saveBorderFd = '/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/Downsampled/linkedIMAGES/alignedTDTMASKS';
for Z = 1:size(maxProjbut1_TDTmasks,3)
    [B,L] = bwboundaries(tp(:,:,Z),'noholes');
    boundaries{Z} = B;
    if ~isempty(B)
        figure; imshow(tp(:,:,Z)); hold on
        for i = 1:length(B)
            boundary = B{i};
            plot(boundary(:,2), boundary(:,1), '-r',  'LineWidth', 1)
        end
        export_fig(fullfile(saveBorderFd, sprintf('ammcSmoothBorder_Z%02d.eps',Z) ))
    end
end


%%
greenfolder = '/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/Downsampled/linkedIMAGES/alignedAnatomicalMaps_GREEN';
redfolder   = '/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/Downsampled/linkedIMAGES/alignedAnatomicalMaps_RED';
masksFolder = '/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/Downsampled/linkedIMAGES/alignedTDTMASKS';
tdtBorderFd = '/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/Downsampled/linkedIMAGES/alignedBorders';

allZs = dataTable.stacksZ + dataTable.zOffsetAlStack;
allZs(2) = allZs(2) - 1; %correct for fly120_4um
dataTable.matrixsZ = allZs;
%%
flyNumSorted = [120 123:125 181:183]; %[119 122 120 123:125 181:183]; %118 a parte,
for f = 1:length(flyNumSorted)
    F = flyNumSorted(f);
    zetas = dataTable.run(dataTable.fly==F);
    for z = 1:sum(dataTable.fly==F)
        tagName  = sprintf('fly%03d_run%02d', F, zetas(z));
        Zaligned = dataTable.matrixsZ(dataTable.fly==F & dataTable.run==zetas(z));
        Zwithin  = dataTable.stacksZ(dataTable.fly==F & dataTable.run==zetas(z));
        
        
        mapGreen = allStacks(allStacks_flyNums==F).alignedstackG(:,:, Zwithin);
        mapName = fullfile(greenfolder, sprintf('alignALL_fly%03d_absZ%02d_%s_mapGreen.png',F, Zaligned, tagName));
        imwrite(mapGreen, mapName);
        dataTable.reg_avgGREEN(dataTable.fly==F & dataTable.run==zetas(z)) = {mapGreen};
        dataTable.reg_pixels2keep(dataTable.fly==F & dataTable.run==zetas(z)) = {[]}; %old is obsolete and new version is not available here
        
        
        if ~isempty(allStacks(allStacks_flyNums==F).alignedstackR)
            mapRed = allStacks(allStacks_flyNums==F).alignedstackR(:,:, Zwithin);
            mapName = fullfile(redfolder, sprintf('alignALL_fly%03d_absZ%02d_%s_mapRed.png',F, Zaligned, tagName));
            imwrite(mapRed, mapName);
            dataTable.reg_avgRED(dataTable.fly==F & dataTable.run==zetas(z)) = {mapRed};
            dataTable.reg_maskRED(dataTable.fly==F & dataTable.run==zetas(z)) = {allStacks(allStacks_flyNums==F).masked_alignedstackR(:,:,Zwithin)};
        end
        
        
        mapWhite = ones(size(mapGreen));
        figure; imshow(mapWhite); hold on
        if ~isempty(allStacks(allStacks_flyNums==F).aligned_Rboundaries)
            bordRed = allStacks(allStacks_flyNums==F).aligned_Rboundaries(Zwithin);
            for ib = 1 : length(bordRed{:})
                boundary = bordRed{1}{ib};
                plot(boundary(:,2), boundary(:,1), '-r')
            end
        end
        
        % add runs borders to the same figure
        iz = find(strcmp(basenames, tagName));
        regData = matfile(fullfile(runfolders{iz}, ['stackRegistration_' basenames{iz} '.mat']));
        
        mapmask = ones([85, 122]); %non-downsampled, cropped image size
        mapmask = imrotate(mapmask, regGeneralData.totAng(1,iz));
        mapmask = imwarp(mapmask,regData.tform,'OutputView',regData.Roriginal,'FillValues',0);
        mapmask = imwarp(mapmask,allStacks(allStacks_flyNums==F).tform,'OutputView',allStacks(allStacks_flyNums==F).Roriginal,'FillValues',0);
        mapmask(mapmask<0.7) = 0;
        mapmask = logical(mapmask);
        [B,L] = bwboundaries(mapmask,'noholes');
        boundary = B{1};
        dataTable.boundaries(dataTable.fly==F & dataTable.run==zetas(z)) = {boundary};
        plot(boundary(:,2), boundary(:,1), '-c')
        
        t = annotation('textbox');
        t.FontSize = 15;
        t.HorizontalAlignment = 'center';
        t.String = tagName;
        t.Interpreter = 'none';
        t.LineStyle = 'none';
        t.Position =  [0.4852    0.1389    0.0857    0.0320];
        t.Color = 'c';
        
        mapName = fullfile(tdtBorderFd, sprintf('alignALL_fly%03d_absZ%02d_%s_OwnRedBound&rectangleRun.eps',F, Zaligned, tagName));
        export_fig(mapName)
        close
    end
    %functionalMapName:
    % alignALLklustsMap_fly01_ordZ01_fly120_run01.png
    %targetName:
    % alignALL + _fly120_run01 + absZ%02d indexing allZs
end


%% 118 only, stack taken from 126
flyNumSorted = 118;
f = 1;
F = flyNumSorted(f);
    zetas = dataTable.run(dataTable.fly==F);
    for z = 1:sum(dataTable.fly==F)
        tagName  = sprintf('fly%03d_run%02d', F, zetas(z));
        Zaligned = dataTable.matrixsZ(dataTable.fly==F & dataTable.run==zetas(z));
        Zwithin  = dataTable.stacksZ(dataTable.fly==F & dataTable.run==zetas(z));
        
        
        mapGreen = allStacks(allStacks_flyNums==126).alignedstackG(:,:, Zwithin);
        mapName = fullfile(greenfolder, sprintf('alignALL_fly%03d_absZ%02d_%s_mapGreen.png',F, Zaligned, tagName));
        imwrite(mapGreen, mapName);
        dataTable.reg_avgGREEN(dataTable.fly==F & dataTable.run==zetas(z)) = {mapGreen};
        dataTable.reg_pixels2keep(dataTable.fly==F & dataTable.run==zetas(z)) = {[]}; %old is obsolete and new version is not available here
        
        
        if ~isempty(allStacks(allStacks_flyNums==126).alignedstackR)
            mapRed = allStacks(allStacks_flyNums==126).alignedstackR(:,:, Zwithin);
            mapName = fullfile(redfolder, sprintf('alignALL_fly%03d_absZ%02d_%s_mapRed.png',F, Zaligned, tagName));
            imwrite(mapRed, mapName);
            dataTable.reg_avgRED(dataTable.fly==F & dataTable.run==zetas(z)) = {mapRed};
            dataTable.reg_maskRED(dataTable.fly==F & dataTable.run==zetas(z)) = {allStacks(allStacks_flyNums==126).masked_alignedstackR(:,:,Zwithin)};
        end
        
        
        mapWhite = ones(size(mapGreen));
        figure; imshow(mapWhite); hold on
        if ~isempty(allStacks(allStacks_flyNums==126).aligned_Rboundaries)
            bordRed = allStacks(allStacks_flyNums==126).aligned_Rboundaries(Zwithin);
            for ib = 1 : length(bordRed{:})
                boundary = bordRed{1}{ib};
                plot(boundary(:,2), boundary(:,1), '-r')
            end
        end
        
        % add runs borders to the same figure
        iz = find(strcmp(basenames, tagName));
        regData = matfile(fullfile(runfolders{iz}, ['stackRegistration_' basenames{iz} '.mat']));
        
        mapmask = ones([85, 122]); %non-downsampled, cropped image size
        mapmask = imrotate(mapmask, regGeneralData.totAng(1,iz));
        mapmask = imwarp(mapmask,regData.tform,'OutputView',regData.Roriginal,'FillValues',0);
        mapmask = imwarp(mapmask,allStacks(allStacks_flyNums==126).tform,'OutputView',allStacks(allStacks_flyNums==126).Roriginal,'FillValues',0);
        mapmask(mapmask<0.7) = 0;
        mapmask = logical(mapmask);
        [B,L] = bwboundaries(mapmask,'noholes');
        boundary = B{1};
        dataTable.boundaries(dataTable.fly==F & dataTable.run==zetas(z)) = {boundary};
        plot(boundary(:,2), boundary(:,1), '-c')
        
        t = annotation('textbox');
        t.FontSize = 15;
        t.HorizontalAlignment = 'center';
        t.String = tagName;
        t.Interpreter = 'none';
        t.LineStyle = 'none';
        t.Position =  [0.4852    0.1389    0.0857    0.0320];
        t.Color = 'c';
        
        mapName = fullfile(tdtBorderFd, sprintf('alignALL_fly%03d_absZ%02d_%s_OwnRedBound&rectangleRun.eps',F, Zaligned, tagName));
        export_fig(mapName)
        close
    end
    %functionalMapName:
    % alignALLklustsMap_fly01_ordZ01_fly120_run01.png
    %targetName:
    % alignALL + _fly120_run01 + absZ%02d indexing allZs




%%
for f = 1: length(allStacks)
    allZGaps_um(f) = -allStacks(f).metaimage.zStepSize/4;
end


table(:,1) = allStacks_flyNums;
table(:,2) = allZGaps_um';

