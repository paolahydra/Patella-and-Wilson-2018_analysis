
%% settings
% start from the 'new' to-be-aligned-stacks. Then add 119 122 and 126 too.
rootname = 'datalist_TDTstacks'; % vai sotto per le altre datalists
rootfolder = 'TDTstacks';

ParentFolder = '/Users/galileo/Dropbox (HMS)/Data';
DataListFolder = fullfile(ParentFolder,rootfolder);
if exist(DataListFolder, 'dir') ~= 7
    mkdir(DataListFolder)
end

% make/get datalist - %provide datalist.mat or otherwise need to DIR##data
if exist(fullfile(DataListFolder, [rootname '.mat']), 'file') ~= 2
    datalist_TDTstacks = uipickfiles('FilterSpec', '/Users/galileo/Dropbox (HMS)/Data/*stack*.mat');
else
    load(fullfile(DataListFolder, [rootname '.mat']));     %DIR#-required list only
end
datalist = eval(rootname); %address of scanimage stacks .tif
clear(rootname)

NaF = length(datalist);

saveName = fullfile(DataListFolder, 'anatomy_alignment_metadata.mat');
if exist(saveName, 'file') == 2 % if part of this datalist has already been processed
    load(saveName)
    flies2process = length(allStacks)+1 : NaF; %default: process all flies which haven't been yet (in the same datalist)
else
    flies2process = 1 : NaF; %defalut: process all the flies in the datalist. New flies are appended at the end of the datalist
end



%% prelimins
for f = flies2process       % this f was indexed as z before.
    disp(f)
    metaimage = scim_openTif(datalist{f});
%     zooms(f) = metaimage.acq.zoomFactor;
%     acqAnglesStacks(f) = metaimage.acq.scanRotation;
    
    [a,b] = fileparts(datalist{f});
    stackR = readTiffPP(fullfile(a, ['preprocessed_R_' b '.tif']), 'double');
    stackG = readTiffPP(fullfile(a, ['preprocessed_G_' b '.tif']), 'double');
    MIJ.createImage(max(stackR,[],3));
    MIJ.createImage(max(stackG,[],3));
    anglesStck_single(f) = input(sprintf('input angle of manual stack rotation to orthogonality: '));
    stackG = imrotate(stackG, -anglesStck_single(f));
    stackR = imrotate(stackR, -anglesStck_single(f));
    clear stackG_ad  stackR_ad
    for z = 1 : size(stackG,3)
        stackG_ad(:,:,z) = adapthisteq(stackG(:,:,z));
        stackG_ad(:,:,z) = imadjust(stackG_ad(:,:,z));
        stackR_ad(:,:,z) = adapthisteq(stackR(:,:,z));
        stackR_ad(:,:,z) = imadjust(stackR_ad(:,:,z));
    end
    MIJ.closeAllWindows();
    
    MIJ.createImage(stackR_ad);
    MIJ.createImage(stackG_ad);
    MPrange = input('slice range to keep for Max Projection ( min:max ): ');
    maxProjectionG = max(stackG(:,:,MPrange),[], 3);
    maxProjectionR = max(stackR(:,:,MPrange),[], 3);
    MIJ.createImage(maxProjectionG)
    MIJ.createImage(maxProjectionR)
    pause()
    MIJ.closeAllWindows();
    
    
    % to be saved
    allStacks(f).flyNum = str2double(datalist{f}(strfind(datalist{f}, 'fly')+3:strfind(datalist{f}, 'fly')+5));
    allStacks(f).stackG = stackG_ad;
    allStacks(f).stackR = stackR_ad;
    allStacks(f).sourceAddress = datalist{f};
    allStacks(f).metaimage = metaimage.acq;
    allStacks(f).manualRotation = anglesStck_single(f);
    allStacks(f).MPrange = MPrange;
    allStacks(f).maxProjectionG = maxProjectionG;
    allStacks(f).maxProjectionR = maxProjectionR;
end

save(saveName, 'allStacks', '-v7.3')

%% add stimFinal main stacks
clear
tdTGal80_old_anatAlign = matfile('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/anatomyImages/anatomy_alignment_metadata.mat');
%tdT+ flies are indexed 6-7-8 (flies 119, 122, 126)

rootname = 'datalist_Gal80stacks'; % exclude 118 for now. Later could show at least part of boundary, or bound ffrom other flies
rootfolder = 'TDTstacks';

ParentFolder = '/Users/galileo/Dropbox (HMS)/Data';
DataListFolder = fullfile(ParentFolder,rootfolder);
if exist(DataListFolder, 'dir') ~= 7
    mkdir(DataListFolder)
end

% make/get datalist - %provide datalist.mat or otherwise need to DIR##data
if exist(fullfile(DataListFolder, [rootname '.mat']), 'file') ~= 2
    datalist_TDTomStacks_stimFinal = uipickfiles('FilterSpec', '/Users/galileo/Dropbox (HMS)/Data/*stack*.mat');
else
    load(fullfile(DataListFolder, [rootname '.mat']));     %DIR#-required list only
end
datalist = eval(rootname); %address of scanimage stacks .tif
clear(rootname)

NaF = length(datalist);

saveName = fullfile(DataListFolder, 'anatomy_alignment_metadata.mat');
if exist(saveName, 'file') == 2 % if part of this datalist has already been processed
    previousAllStacks = load(saveName); %previous allStacks, the concatenate them
    flies2process = 1 : NaF; % only conider this subset anyways
else
    flies2process = 1 : NaF; %defalut: process all the flies in the datalist. New flies are appended at the end of the datalist
end
%%

% clear allStacks  % do run after making sure allstacks has been save to disk

for f = flies2process       % this f was indexed as z before.
    metaimage = scim_openTif(datalist{f});
%     zooms(f) = metaimage.acq.zoomFactor;
%     acqAnglesStacks(f) = metaimage.acq.scanRotation;
    
    [a,b] = fileparts(datalist{f});
%     stackR = readTiffPP(fullfile(a, ['preprocessed_R_' b '.tif']), 'double');
%     stackG = readTiffPP(fullfile(a, ['preprocessed_G_' b '.tif']), 'double');
stackG = readTiffPP(fullfile(a, ['preprocessed_' b '.tif']), 'double');

%     MIJ.createImage(max(stackR,[],3));
    MIJ.createImage(max(stackG,[],3));
    anglesStck_single(f) = input(sprintf('input angle of manual stack rotation to orthogonality: '));
    stackG = imrotate(stackG, -anglesStck_single(f));
%     stackR = imrotate(stackR, -anglesStck_single(f));
    clear stackG_ad  stackR_ad
    for z = 1 : size(stackG,3)
        stackG_ad(:,:,z) = adapthisteq(stackG(:,:,z));
        stackG_ad(:,:,z) = imadjust(stackG_ad(:,:,z));
%         stackR_ad(:,:,z) = adapthisteq(stackR(:,:,z));
%         stackR_ad(:,:,z) = imadjust(stackR_ad(:,:,z));
    end
    MIJ.closeAllWindows();
    
%     MIJ.createImage(stackR_ad);
    MIJ.createImage(stackG_ad);
    MPrange = input('slice range to keep for Max Projection ( min:max ): ');
    maxProjectionG = max(stackG(:,:,MPrange),[], 3);
%     maxProjectionR = max(stackR(:,:,MPrange),[], 3);
    MIJ.createImage(maxProjectionG)
%     MIJ.createImage(maxProjectionR)
    pause()
    MIJ.closeAllWindows();
    
    
    % to be saved
    allStacks(f).flyNum = str2double(datalist{f}(strfind(datalist{f}, 'fly')+3:strfind(datalist{f}, 'fly')+5));
    allStacks(f).stackG = stackG_ad;
%     allStacks(f).stackR = stackR_ad;
allStacks(f).stackR = [];
    allStacks(f).sourceAddress = datalist{f};
    allStacks(f).metaimage = metaimage.acq;
    allStacks(f).manualRotation = anglesStck_single(f);
    allStacks(f).MPrange = MPrange;
    allStacks(f).maxProjectionG = maxProjectionG;
%     allStacks(f).maxProjectionR = maxProjectionR;
allStacks(f).maxProjectionR = [];
end


% for f  = 1:length(allStacks)
%     figure(f); 
%     imagesc(allStacks(f).maxProjectionR), axis image, colormap(gray);
% end


% HERE made and saved reference_MP from fly 119


%% MODULE 2: align to new_ref_fly119
load('newReferenceStack_fly119_cropped.mat');
movingPoints = [];
fixedPoints = [];  % check if no errors for f==1
figure;
%%
for f = 1: length(allStacks)
    disp(f)
    MP_stack = allStacks(f). maxProjectionG;
    [movingPoints, fixedPoints] = cpselect(MP_stack./max(MP_stack(:)), reference_MP, 'Wait', true);
    tform = fitgeotrans(movingPoints,fixedPoints,'nonreflectivesimilarity');
    tformInv = invert(tform);
    Tinv = tformInv.T;
    ss = Tinv(2,1);
    sc = Tinv(1,1);
    scale_recovered = sqrt(ss*ss + sc*sc)
    theta_recovered = atan2(ss,sc)*180/pi
    Roriginal = imref2d(size(reference_MP));
    recovered = imwarp(MP_stack,tform,'OutputView',Roriginal);
    h_montage=figure; imshowpair(reference_MP,recovered,'montage')
    h_composite=figure; imshowpair(reference_MP,recovered)
    pause
    while 1
        choice = questdlg('Redo?','repeat?','Yes','No','No');
        switch choice
            case 'No'
                break
            case 'Yes'
                [movingPoints, fixedPoints] = cpselect(MP_stack./max(MP_stack(:)), reference_MP, movingPoints, fixedPoints, 'Wait', true);
                tform = fitgeotrans(movingPoints,fixedPoints,'nonreflectivesimilarity');
                tformInv = invert(tform);
                Tinv = tformInv.T;
                ss = Tinv(2,1);
                sc = Tinv(1,1);
                scale_recovered = sqrt(ss*ss + sc*sc)
                theta_recovered = atan2(ss,sc)*180/pi
                Roriginal = imref2d(size(reference_MP));
                recovered = imwarp(MP_stack,tform,'OutputView',Roriginal);
                h_montage=figure; imshowpair(reference_MP,recovered,'montage')
                h_composite=figure; imshowpair(reference_MP,recovered)
                pause
        end
    end
    
    export_fig(h_montage, fullfile(DataListFolder, sprintf('stack2fly119ref_fly%d_montage.png', allStacks(f).flyNum)))
    export_fig(h_composite, fullfile(DataListFolder, sprintf('stack2fly119ref_fly%d_composite.png', allStacks(f).flyNum)),'-m2')
    allStacks(f).tform = tform;
    allStacks(f).recovered = recovered;
    allStacks(f).movingPoints = movingPoints;
    allStacks(f).fixedPoints = fixedPoints;
    allStacks(f).Roriginal = Roriginal;
    
    save(fullfile(DataListFolder, sprintf('stack2fly119ref_fly%d.mat',allStacks(f).flyNum)), ...
        'tform', 'recovered', 'movingPoints', 'fixedPoints', 'Roriginal')
    

    %% apply the same transformation to the entire stack
    alignedstackG = zeros([size(recovered), size(allStacks(f).stackG,3)]);
%     alignedstackR = zeros([size(recovered), size(allStacks(f).stackG,3)]);
    for z = 1: size(allStacks(f).stackG,3)
        alignedstackG(:,:,z) = imwarp(allStacks(f).stackG(:,:,z),tform,'OutputView',Roriginal);
%         alignedstackR(:,:,z) = imwarp(allStacks(f).stackR(:,:,z),tform,'OutputView',Roriginal);
    end
    allStacks(f).alignedstackG = alignedstackG;
%     allStacks(f).alignedstackR = alignedstackR;
allStacks(f).alignedstackR = [];

    %% tdtomato mask
    
%     this has been edited in temporary_smooth_refineTDTMasks.m
edit temporary_smooth_refineTDTMasks.m
    
    
    
%     se = strel('square',7);
%     arbitraryThreshold = 0.6; %between 0 and 1
%     sizeMAP = size(reference_MP);
%     masked_alignedstackR = zeros([sizeMAP, size(allStacks(f).stackG,3)]);
%     boundaries = cell([1, size(allStacks(f).stackG,3)]);
%     for z = allStacks(f).MPrange
%         tdTomPositivePxls = false(sizeMAP);
%         tdTomPositivePxls(alignedstackR(:,:,z)>=arbitraryThreshold) = 1;
%         tdTomPositivePxls = logical(tdTomPositivePxls);
%         tdTomPositivePxls = imdilate(tdTomPositivePxls, se);
%         tdTomPositivePxls = imclose(tdTomPositivePxls, strel('square',9));
% %         figure; imagesc(tdTomPositivePxls); axis image; colormap(gray);
%         masked_alignedstackR(:,:,z) = tdTomPositivePxls;
%        %         figure; imagesc(masked_alignedstackR(:,:,z)); axis image; colormap(gray);
%         [B,L] = bwboundaries(tdTomPositivePxls,'noholes');
%         boundaries{z} = B;
% %         boundary = B{1};
% %         plot(boundary(:,2), boundary(:,1), 'LineWidth', 1)
%     end
%     allStacks(f).masked_alignedstackR = masked_alignedstackR;
%     allStacks(f).aligned_Rboundaries = boundaries;
allStacks(f).masked_alignedstackR = [];
allStacks(f).aligned_Rboundaries = [];
    
end        
% save(saveName, 'allStacks', '-v7.3')


if exist('previousAllStacks', 'var')
    allStacks = cat(2, previousAllStacks.allStacks, allStacks);
end
save(saveName, 'allStacks', '-v7.3')



%% add functional rectangles, add tdT boundaries, if any. SAVE


%% align Gal80 flies to fly 119 (need to estimate MAXPROJ range first).






%% only serial stack with rectangles
for f = flies2process       % this f was indexed as z before.
    countfig = 1;
    flyfolder = flyfoldersUnique{f};
    [~, flyname] = fileparts(flyfolder);
    flyname = flyname(1:6);
    stack2fly150_info = matfile(fullfile(flyfolder, ['stack2fly150sStack_' flyname '.mat']));
    
    stackName = fullfile(flyfolder, 'orto_preprocessed_*.tif'); 
    importlist = dir(stackName);
    if strcmp(flyname, 'fly150')
        alignStackName = fullfile(flyfolder, importlist(1).name);
    else
        alignStackName = fullfile(flyfolder, ['stack2fly150sStack_' importlist(1).name]);
    end
    alignedstack = readTiffPP(alignStackName, 'double');
    [fig_alStackSer, ax_alStackSer] = plotSerialStack_imagesc(alignedstack,[],gray);
    
    ZZ = strfind(datalist, flyfolder);
    zs_fly = zeros(size(ZZ));
    for i = 1:length(zs_fly)
        if ~isempty(ZZ{i})
            zs_fly(i) = 1; % indices into datalist relative to fly f
        end
    end
    zs_fly = logical(zs_fly);
    %% sort runs from ventral to dorsal
    zetas = find(zs_fly);
    sliceNums = regGeneralData.alignedSliceNumbers(1,zetas);
    [sliceNums,b] = sort(sliceNums);
    zetas = zetas(b);
    for iz = 1:sum(zs_fly)
        z = zetas(iz);
        sliceNum = sliceNums(iz);
        mapmask = ones(size(aZ{z}.pixels2keep));
        % load specific run's transformation
        regData = matfile(fullfile(sets.analysislist{z}, ['stackRegistration_' basenames{z} '.mat']));
        mapmask = imrotate(mapmask, regGeneralData.totAng(1,z));
        mapmask = imwarp(mapmask,regData.tform,'OutputView',regData.Roriginal);
        if ~strcmp(flyname, 'fly150')
            mapmask = imwarp(mapmask,stack2fly150_info.tform,'OutputView',stack2fly150_info.Roriginal);
        end
        mapmask = logical(mapmask);
        [B,L] = bwboundaries(mapmask,'noholes');
        boundary = B{1};
        
        
        axes(ax_alStackSer(sliceNum)); hold on
        plot(boundary(:,2), boundary(:,1), 'LineWidth', 1)
        ax_alStackSer(sliceNum).YLabel.String = cat(1, ax_alStackSer(regGeneralData.alignedSliceNumbers(1,z)).YLabel.String, basenames(z));
        ax_alStackSer(sliceNum).YLabel.Color = [0, 0.75, 1];
        ax_alStackSer(regGeneralData.alignedSliceNumbers(1,z)).YLabel.Interpreter = 'none';
        
%         % there was an error in saving the full stack for fly 82
%         figure;  imagesc(alignedstack(:,:,sliceNum)); axis image; axis off; colormap(gray); hold on
%         plot(boundary(:,2), boundary(:,1), 'LineWidth', 1)
%         sStName = fullfile(flyfolder, sprintf('serial_stack2fly150sStack_%s_z%d', flyname,iz));
%         export_fig(sStName, '-eps');
    end
    
    disp('---------------------------------------------------------------------------')
    disp('move figures to the left, then press any key here to continue:')
    pause
    disp('Resuming')
    fig_alStackSer.Position = [-900           5        2280         801];
    sStName = fullfile(flyfolder, ['serial_stack2fly150sStack_' flyname '_StackWithRectangles']);
    export_fig(fig_alStackSer, sStName, '-eps');
    close(fig_alStackSer)
end


%% only serial stack RED with rectangles for tdTom+
for f = [7 8]       % this f was indexed as z before.
    countfig = 1;
    flyfolder = flyfoldersUnique{f};
    [~, flyname] = fileparts(flyfolder);
    flyname = flyname(1:6);
    stack2fly150_info = matfile(fullfile(flyfolder, ['stack2fly150sStack_' flyname '.mat']));
    
    stackName = fullfile(flyfolder, 'preprocessed_R_*.tif'); 
    importlist = dir(stackName);
    stackName = fullfile(flyfolder, importlist(1).name);
    redStack = readTiffPP(stackName, 'double');
    % rotate:
    redStackRot = imrotate(redStack, -1*regGeneralData.anglesStck_single(1, f));
    al2fly150_redStack = imwarp(redStackRot,stack2fly150_info.tform,'OutputView',stack2fly150_info.Roriginal); %,'FillValues',[1,1,1]);
    alignStackName = fullfile(flyfolder, ['stack2fly150sStack_' importlist(1).name]);
 writetiff(al2fly150_redStack, alignStackName);
 
 
    [fig_alStackSer, ax_alStackSer] = plotSerialStack_imagesc(al2fly150_redStack,[],gray);
    
    ZZ = strfind(datalist, flyfolder);
    zs_fly = zeros(size(ZZ));
    for i = 1:length(zs_fly)
        if ~isempty(ZZ{i})
            zs_fly(i) = 1; % indices into datalist relative to fly f
        end
    end
    zs_fly = logical(zs_fly);
    %% sort runs from ventral to dorsal
    zetas = find(zs_fly);
    sliceNums = regGeneralData.alignedSliceNumbers(1,zetas);
    [sliceNums,b] = sort(sliceNums);
    zetas = zetas(b);
    for iz = 1:sum(zs_fly)
        z = zetas(iz);
        sliceNum = sliceNums(iz);
        mapmask = ones(size(aZ{z}.pixels2keep));
        % load specific run's transformation
        regData = matfile(fullfile(sets.analysislist{z}, ['stackRegistration_' basenames{z} '.mat']));
        mapmask = imrotate(mapmask, regGeneralData.totAng(1,z));
        mapmask = imwarp(mapmask,regData.tform,'OutputView',regData.Roriginal);
        mapmask = imwarp(mapmask,stack2fly150_info.tform,'OutputView',stack2fly150_info.Roriginal);
        mapmask = logical(mapmask);
        [B,L] = bwboundaries(mapmask,'noholes');
        boundary = B{1};
        axes(ax_alStackSer(sliceNum)); hold on
        plot(boundary(:,2), boundary(:,1), 'LineWidth', 1)
        ax_alStackSer(sliceNum).YLabel.String = cat(1, ax_alStackSer(regGeneralData.alignedSliceNumbers(1,z)).YLabel.String, basenames(z));
        ax_alStackSer(sliceNum).YLabel.Color = [0, 0.75, 1];
        ax_alStackSer(regGeneralData.alignedSliceNumbers(1,z)).YLabel.Interpreter = 'none';
    end
    disp('---------------------------------------------------------------------------')
    disp('move figures to the left, then press any key here to continue:')
    pause
    disp('Resuming')
    fig_alStackSer.Position = [-900           5        2280         801];
    sStName = fullfile(flyfolder, ['serial_stack2fly150sStack_' flyname '_RED_StackWithRectangles']);
    export_fig(fig_alStackSer, sStName, '-eps');
    close(fig_alStackSer)
end

 
%% extra to detect step size in each stack
for f = flies2process
    flyfolder = flyfoldersUnique{f};
    [~, flyname] = fileparts(flyfolder);
    flyname = flyname(1:6);
    disp(flyname)
    %% get each stack's zoom and any acquisition angle FAREI MEGLIO A SCEGLIERLO MANUALMENTE
    %find stack
    clear metaimage
    countStacks = 0;
    stackRef = dir(fullfile(flyfolder, 'orto_preprocessed_*stack*.tif'));
    stackl = dir(fullfile(flyfolder, '*stack*.tif'));
    i = 1;
    while i <= length(stackl) 
        res = strfind(stackRef.name(length('orto_preprocessed_')+1:end), stackl(i).name);
        if ~isempty(res)
            disp(stackl(i).name)
            stackslist{f} = fullfile(flyfolder,stackl(i).name);
            try 
                metaimage = scim_openTif(stackslist{f});
                countStacks = countStacks+1;
            catch ME
                disp('nope')
                if ~strcmp(ME.identifier, 'MATLAB:imagesci:Tiff:unableToOpenFile')
                    rethrow(ME)
                end
            end 
        end
        i=i+1;
    end  
    if countStacks~=1
        disp('----------------')
        disp('reference stack name:')
        disp(stackRef.name)
        stackl.name
        disp('----------------')
        error('check stacks:')
    end
    % get data
    if metaimage.motor.motorZEnable
        zfactor = 1;
    else
        zfactor = 0.25;
    end
    regGeneralData.zStepSize_um(1,f) = abs(metaimage.acq.zStepSize*zfactor);    
end




%% addendum to purge a datalist
% % delete from datalist and regGenData - ONLY WORKS FOR ONE ITERATION OF Z
% origLengt = length(datalist);
% z = 40;
% datalist = datalist([1:z-1,z+1:origLengt]);
% %%
% regGeneralData.alignedSliceNumbers  = alignedSliceNumbers([1:z-1,z+1:origLengt]);
% regGeneralData.allLandmarks         = allLandmarks([1:z-1,z+1:origLengt]);
% regGeneralData.anglesRuns           = anglesRuns([1:z-1,z+1:origLengt]);
% regGeneralData.anglesStck           = anglesStck([1:z-1,z+1:origLengt]);
% regGeneralData.flyfolders           = flyfolders([1:z-1,z+1:origLengt]);
% regGeneralData.runfolders           = runfolders([1:z-1,z+1:origLengt]);
% regGeneralData.totAng               = totAng([1:z-1,z+1:origLengt]);
% 
% %% save datalist
% eval(rootname) = datalist;
% DataSaveFolder = fullfile(ParentFolder, rootfolder);
% save(fullfile(DataSaveFolder, [rootname '.mat']), 'datalist_WINDS_86', '-v7.3');

%%
for i = 2:40
    close(figure(i))
end
