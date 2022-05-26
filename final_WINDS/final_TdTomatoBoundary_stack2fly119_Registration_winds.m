%these are all Gal80
%% settings - OK
% start from the 'new' to-be-aligned-stacks (100). Then add info from flies 86
load_winds_downsampled; clear R;
rootfolder = 'WINDS_100/anatomyImages';
ParentFolder = '/Users/galileo/Dropbox (HMS)/Data';
DataListFolder = fullfile(ParentFolder,rootfolder);
runfolders = sets.analysislist;

saveName = fullfile(DataListFolder, 'anatomy_alignment_metadata.mat');

if exist(saveName, 'file') == 2 % if part of this datalist has already been processed
    IntegrateOnly = 1;
    regGeneralData = matfile(saveName,'Writable',true); %single z planes to their own stack, general infos
    flies2process = 1:NaF;
%     flies2process = length(regGeneralData.flyfoldersUnique)+1 : NaF; %default: process all flies which haven't been yet
    allStacks = regGeneralData.allStacks;
else
    IntegrateOnly = 0;
    flies2process = 1 : NaF; %defalut: process all the flies in the datalist. New flies are appended at the end of the datalist
    % save(saveName, 'allStacks', 'anglesRuns', 'anglesStck', 'totAng', 'zooms', 'anglesStck_single','flyfoldersUnique','flyfolders','runfolders', '-v7.3')
end

[flies2process,~] = listdlg('PromptString','Select flies to process:',...
             'ListString',flyfoldersUnique, 'InitialValue', flies2process, 'ListSize', [300 340] );

% it's dangerous to keep only the relevant zetas. I will have to use he indices of processed runs within new datalist
%%
anglesRuns = zeros(1, NaZ);
stackslist = cell(1,NaF);
acqAnglesStacks = zeros(1, NaF);
anglesStck_single = zeros(1, NaF); %this is the angle of the manual rotation, and is eneterd manually
anglesStck = nan(1, NaZ);
alignedSliceNumbers = zeros(1, NaZ);
rightIpsiSide = nan(1, NaZ);

%% una tantum prelims
% % (una tantum make/get stack datalist - %provide datalist.mat or otherwise need to DIR##data
% datalist_WINDS = datalist;
% rootname = 'datalist_allWINDS_stacks';
% if exist(DataListFolder, 'dir') ~= 7
%     mkdir(DataListFolder)
% end
% 
% if exist(fullfile(DataListFolder, [rootname '.mat']), 'file') ~= 2
%     for f = 1:NaF
%         datalist_allWINDS_stacks{f} = uipickfiles('FilterSpec', [flyfoldersUnique{f} '/*stack*.tif'], 'Output', 'char');
%     end
%     save(fullfile(DataListFolder, [rootname '.mat']), 'datalist_allWINDS_stacks')
% else
%     load(fullfile(DataListFolder, [rootname '.mat']));     %DIR#-required list only
% end
% datalist = eval(rootname); %address of scanimage stacks .tif
% clear(rootname) 
% 
% % prelimins (100)
% warning off
% for f = 1:3       % this f was indexed as z before.
%     disp(f)
%     metaimage = scim_openTif(datalist{f});
% %     zooms(f) = metaimage.acq.zoomFactor;
% %     acqAnglesStacks(f) = metaimage.acq.scanRotation;
%     
%     [a,b] = fileparts(datalist{f});
%     try
%         stackG = readTiffPP(fullfile(a, ['preprocessed_' b '.tif']), 'double');
%     catch err
%         d = dir(fullfile(a, 'preprocessed_*.tif'));
%         if length(d) == 1
%             stackG = readTiffPP(fullfile(a, d.name), 'double');
%         else
%             throw(err)
%         end
%     end
%     MIJ.createImage(max(stackG,[],3));
%     anglesStck_single(f) = input(sprintf('input angle of manual stack rotation to orthogonality: '));
%     stackG = imrotate(stackG, -anglesStck_single(f));
% 
%     clear stackG_ad 
%     for z = 1 : size(stackG,3)
%         stackG_ad(:,:,z) = adapthisteq(stackG(:,:,z));
%         stackG_ad(:,:,z) = imadjust(stackG_ad(:,:,z));
%     end
%     MIJ.closeAllWindows();
%     
%     MIJ.createImage(stackG_ad);
%     MPrange = input('slice range to keep for Max Projection ( min:max ): ');
%     maxProjectionG = max(stackG(:,:,MPrange),[], 3);
% 
%     MIJ.createImage(maxProjectionG)
%     pause()
%     MIJ.closeAllWindows();
%     
%     
%     % to be saved
%     allStacks(f).flyNum = str2double(datalist{f}(strfind(datalist{f}, 'fly')+3:strfind(datalist{f}, 'fly')+5));
%     allStacks(f).stackG = stackG_ad;
%     allStacks(f).sourceAddress = datalist{f};
%     allStacks(f).metaimage = metaimage.acq;
%     allStacks(f).manualRotation = anglesStck_single(f);
%     allStacks(f).MPrange = MPrange;
%     allStacks(f).maxProjectionG = maxProjectionG;
% end
% save(saveName, 'allStacks', '-v7.3')
% 
% 
% % prelimins (86)
% oldAlign86 = matfile('/Users/galileo/Dropbox (HMS)/Data/WINDS_86/anatomyImages/anatomy_alignment_metadata.mat');
% for f = 4:7       % this f was indexed as z before.
%     disp(f)
%     metaimage = scim_openTif(datalist{f});
% %     zooms(f) = metaimage.acq.zoomFactor;
% %     acqAnglesStacks(f) = metaimage.acq.scanRotation;
%     
%     [a,b] = fileparts(datalist{f});
%     try
%         stackG = readTiffPP(fullfile(a, ['preprocessed_' b '.tif']), 'double');
%     catch err
%         d = dir(fullfile(a, 'preprocessed_*.tif'));
%         if length(d) == 1
%             stackG = readTiffPP(fullfile(a, d.name), 'double');
%         else
%             throw(err)
%         end
%     end
% %     MIJ.createImage(max(stackG,[],3));
%     anglesStck_single(f) = oldAlign86.anglesStck_single(1,f-3); %input(sprintf('input angle of manual stack rotation to orthogonality: '));
%     stackG = imrotate(stackG, -anglesStck_single(f));
% 
%     clear stackG_ad 
%     for z = 1 : size(stackG,3)
%         stackG_ad(:,:,z) = adapthisteq(stackG(:,:,z));
%         stackG_ad(:,:,z) = imadjust(stackG_ad(:,:,z));
%     end
%     MIJ.closeAllWindows();
%     
%     MIJ.createImage(stackG_ad);
%     MPrange = input('slice range to keep for Max Projection ( min:max ): ');
%     maxProjectionG = max(stackG(:,:,MPrange),[], 3);
% 
%     MIJ.createImage(maxProjectionG)
%     pause()
%     MIJ.closeAllWindows();
%     
%     
%     % to be saved
%     allStacks(f).flyNum = str2double(datalist{f}(strfind(datalist{f}, 'fly')+3:strfind(datalist{f}, 'fly')+5));
%     allStacks(f).stackG = stackG_ad;
%     allStacks(f).sourceAddress = datalist{f};
%     allStacks(f).metaimage = metaimage.acq;
%     allStacks(f).manualRotation = anglesStck_single(f);
%     allStacks(f).MPrange = MPrange;
%     allStacks(f).maxProjectionG = maxProjectionG;
% end
% save(saveName, 'allStacks', '-v7.3')

anglesStck_single = cat(2,allStacks.manualRotation); % bug fix

%% MODULE 1 (100): functional runs alignment
figure(1)
cab(1)
for f = flies2process
    movingPoints=[];
    fixedPoints=[];
    flyfolder = flyfoldersUnique{f};
    [~, flyname] = fileparts(flyfolder);
    flyname = flyname(1:6);
    ZZ = strfind(datalist, flyfolder);
    zs_fly = zeros(size(ZZ));
    for i = 1:length(zs_fly)
        if ~isempty(ZZ{i})
            zs_fly(i) = 1; % indices into datalist relative to fly f
        end
    end
    zs_fly = logical(zs_fly);
    
    %%  get each funcional run's acquisition angle
    for z = find(zs_fly)
        metaimage = aZ{z}.metaim;
        anglesRuns(z) = metaimage.scanRotation;
    end
    
    %% match rotation angle to stack
    acqAnglesStacks(f) = allStacks(f).metaimage.scanRotation;
    anglesStck(zs_fly) = acqAnglesStacks(f) + anglesStck_single(f);  %this has been used incorrectly for flies 100.
    totAng = anglesRuns-anglesStck; %this will be updated for each fly
    
    %% problem 1: align each z plane in its own stack
    % fly118 gets aligned to stack of fly 126
    % stackName = fullfile('/Users/galileo/Dropbox/Data/fly126_PP', 'orto_preprocessed_*.tif');
    % importlist = dir(stackName);
    % stackName2 = fullfile('/Users/galileo/Dropbox/Data/fly126_PP', importlist(1).name);
    % stack = readTiffPP(stackName2, 'double');
    
    stack = allStacks(f).stackG;
    for z = find(zs_fly)
        disp(basenames{z})
        runfolder = runfolders{z};
        
        %% manually determine slice and control points. Save all relevant data for future use
        avgName = fullfile(runfolder, ['AVG_*' basenames{z} '*.tif']);
        importlist = dir(avgName);
        avgName = fullfile(runfolder, importlist(1).name);
        avgrun = readTiffPP(avgName, 'double');
        avgrun(end,:) = [];         % same cropping
        avgrun(:,end-5:end) = [];   % same cropping
        avgrun = imrotate(avgrun, totAng(z));
        avgrun = mat2gray(avgrun);
        avgrun = imadjust(avgrun);
        avgrun = adapthisteq(avgrun);
        
        if alignedSliceNumbers(z) == 0
            MIJ.createImage(avgrun);
            MIJ.run('In [+]')
            MIJ.run('In [+]')
            MIJ.run('In [+]')
            MIJ.createImage(stack);
            MIJ.run('In [+]')
            %% 
            alignedSliceNumbers(z) = input('stack slice number:');
            rightIpsiSide(z) = input('right-left. Right? (1-0):');
            %% 
        end
        MIJ.closeAllWindows();
        target = stack(:,:,alignedSliceNumbers(z)); %so far manually determined         % original
        
        %% try to recover scale and translation within slice
        % manually pick reference points
        %     target = imadjust(target);
        %     target = adapthisteq(target);
        if isempty(movingPoints) || isempty(fixedPoints)
            [movingPoints, fixedPoints] = cpselect(3, 2, avgrun, target, 'Wait', true);
        else
            [movingPoints, fixedPoints] = cpselect(3, 2, avgrun, target, movingPoints, fixedPoints, 'Wait', true);
        end
%         [movingPoints, fixedPoints] = cpselect_PP_zoomHijacked(avgrun, target, movingPoints, fixedPoints, 'Wait', true);
        %
        moving_pts_adj= cpcorr(movingPoints, fixedPoints, avgrun, target);
        tform = fitgeotrans(moving_pts_adj,fixedPoints,'nonreflectivesimilarity');
        tformInv = invert(tform);
        Tinv = tformInv.T;
        ss = Tinv(2,1);
        sc = Tinv(1,1);
        scale_recovered = sqrt(ss*ss + sc*sc)
        theta_recovered = atan2(ss,sc)*180/pi
        Roriginal = imref2d(size(target));
        recovered = imwarp(avgrun,tform,'OutputView',Roriginal);
        
        h_montage=figure; imshowpair(target,recovered,'montage')
        h_composite=figure; imshowpair(target,recovered)
        title(sprintf('%s - z: %d',basenames{z}, z), 'Interpreter', 'none')
        while 1
            choice = questdlg('Redo?','repeat?','Yes','No','No');
            switch choice
                case 'No'
                    break
                case 'Yes'
                    [movingPoints, fixedPoints] = cpselect(3, 2, avgrun, target, movingPoints, fixedPoints, 'Wait', true);
                    %
                    allLandmarks(z).movingPoints = movingPoints;        % To be saved
                    allLandmarks(z).fixedPoints = fixedPoints;        % To be saved
                    moving_pts_adj= cpcorr(movingPoints, fixedPoints, avgrun, target);
                    tform = fitgeotrans(moving_pts_adj,fixedPoints,'nonreflectivesimilarity');
                    tformInv = invert(tform);
                    Tinv = tformInv.T;
                    ss = Tinv(2,1);
                    sc = Tinv(1,1);
                    scale_recovered = sqrt(ss*ss + sc*sc)
                    theta_recovered = atan2(ss,sc)*180/pi
                    Roriginal = imref2d(size(target));
                    recovered = imwarp(avgrun,tform,'OutputView',Roriginal);
                    
                    h_montage=figure; imshowpair(target,recovered,'montage')
                    h_composite=figure; imshowpair(target,recovered)
                    title(sprintf('%s - z: %d',basenames{z}, z), 'Interpreter', 'none')
            end
        end
        allLandmarks(z).movingPoints = moving_pts_adj;      % To be saved
        allLandmarks(z).fixedPoints = fixedPoints;        % To be saved
        
        export_fig(h_montage, fullfile(runfolder, ['stackRegistration_' basenames{z} '_montage']), '-png', '-m2.5')
        export_fig(h_composite,fullfile(runfolder, ['stackRegistration_' basenames{z} '_composite']),'-tiff', '-CMYK', '-m1')
        save(fullfile(runfolder, ['stackRegistration_' basenames{z} '.mat']), ...
            'tform', 'recovered', 'target', 'avgName', 'Roriginal')
    end
    %% save all data up to the current fly
    if IntegrateOnly % if part of this datalist has already been processed
        regGeneralData.alignedSliceNumbers(1,find(zs_fly)) = alignedSliceNumbers(1,find(zs_fly)); % it works...
        regGeneralData.rightIpsiSide(1,find(zs_fly)) = rightIpsiSide(1,find(zs_fly)); 
        regGeneralData.allLandmarks(1,find(zs_fly)) = allLandmarks(1,find(zs_fly));
        regGeneralData.anglesRuns(1,find(zs_fly)) = anglesRuns(1,find(zs_fly));
        regGeneralData.anglesStck(1,find(zs_fly)) = anglesStck(1,find(zs_fly));
        regGeneralData.flyfolders(1,find(zs_fly)) = flyfolders(1,find(zs_fly));
        regGeneralData.runfolders(1,find(zs_fly)) = runfolders(1,find(zs_fly));
        regGeneralData.totAng(1,find(zs_fly)) = totAng(1,find(zs_fly));
        regGeneralData.anglesStck_single(1,f) = anglesStck_single(1,f);
        regGeneralData.flyfoldersUnique(1,f) = flyfoldersUnique(1,f);
    else
        save(saveName, 'allStacks', 'anglesRuns', 'anglesStck', 'totAng', 'rightIpsiSide', 'anglesStck_single','flyfoldersUnique','flyfolders','runfolders', 'allLandmarks', 'alignedSliceNumbers', '-v7.3')
        IntegrateOnly = 1;
        regGeneralData = matfile(saveName,'Writable',true); %single z planes to their own stack, general infos
    end
regGeneralData

    %% second pass?
    runfolders{zs_fly}
    fprintf('\nZetas this fly:\t %s\n\n', num2str(find(zs_fly)))
    figure(1)
    cab(1)
    for z = find(zs_fly)
        runfolder = runfolders{z};
        target = stack(:,:,alignedSliceNumbers(z)); %so far manually determined         % original
        load(fullfile(runfolder, ['stackRegistration_' basenames{z} '.mat']),'recovered')
        h_composite=figure('WindowStyle', 'normal', 'Position', [1921         -40        1440         805]); imshowpair(target,recovered)
        title(sprintf('%s - z: %d',basenames{z}, z), 'Interpreter', 'none')
        reapplyTransform = input(sprintf('Current z = %d. Reapply transformation from z # (0 for no):  ', z));
        if reapplyTransform ~= 0
            avgName = fullfile(runfolder, ['AVG_*' basenames{z} '*.tif']);
            importlist = dir(avgName);
            avgName = fullfile(runfolder, importlist(1).name);
            avgrun = readTiffPP(avgName, 'double');
            avgrun(end,:) = [];         % same cropping
            avgrun(:,end-5:end) = [];   % same cropping
            avgrun = imrotate(avgrun, totAng(z));
            avgrun = mat2gray(avgrun);
            avgrun = imadjust(avgrun);
            avgrun = adapthisteq(avgrun);
            movingPoints = allLandmarks(reapplyTransform).movingPoints;
            fixedPoints = allLandmarks(reapplyTransform).fixedPoints;
%             [movingPoints, fixedPoints] = cpselect(3, 2, avgrun, target, movingPoints, fixedPoints, 'Wait', true);  
            moving_pts_adj= cpcorr(movingPoints, fixedPoints, avgrun, target);
            tform = fitgeotrans(moving_pts_adj,fixedPoints,'nonreflectivesimilarity');
            tformInv = invert(tform);
            Tinv = tformInv.T;
            ss = Tinv(2,1);
            sc = Tinv(1,1);
            scale_recovered = sqrt(ss*ss + sc*sc)
            theta_recovered = atan2(ss,sc)*180/pi
            Roriginal = imref2d(size(target));
            recovered = imwarp(avgrun,tform,'OutputView',Roriginal);
            
            h_composite=figure; imshowpair(target,recovered)
            title(sprintf('%s - z: %d',basenames{z}, z), 'Interpreter', 'none')
            while 1
                choice = questdlg('Redo?','repeat?','Yes','No', 'KeepOriginal', 'No');
                switch choice
                    case 'No'
                        allLandmarks(z).movingPoints = moving_pts_adj;      % To be saved
                        allLandmarks(z).fixedPoints = fixedPoints;        % To be saved
                        
                        export_fig(h_montage, fullfile(runfolder, ['stackRegistration_' basenames{z} '_montage']), '-png', '-m2.5')
                        export_fig(h_composite,fullfile(runfolder, ['stackRegistration_' basenames{z} '_composite']),'-tiff', '-CMYK', '-m1')
                        save(fullfile(runfolder, ['stackRegistration_' basenames{z} '.mat']), ...
                            'tform', 'recovered', 'target', 'avgName', 'Roriginal')
                        break
                        
                    case 'Yes'
                        [movingPoints, fixedPoints] = cpselect(3, 2, avgrun, target, movingPoints, fixedPoints, 'Wait', true);
                        %
                        allLandmarks(z).movingPoints = movingPoints;        % To be saved
                        allLandmarks(z).fixedPoints = fixedPoints;        % To be saved
                        moving_pts_adj= cpcorr(movingPoints, fixedPoints, avgrun, target);
                        tform = fitgeotrans(moving_pts_adj,fixedPoints,'nonreflectivesimilarity');
                        tformInv = invert(tform);
                        Tinv = tformInv.T;
                        ss = Tinv(2,1);
                        sc = Tinv(1,1);
                        scale_recovered = sqrt(ss*ss + sc*sc)
                        theta_recovered = atan2(ss,sc)*180/pi
                        Roriginal = imref2d(size(target));
                        recovered = imwarp(avgrun,tform,'OutputView',Roriginal);
                        
                        h_montage=figure; imshowpair(target,recovered,'montage')
                        h_composite=figure; imshowpair(target,recovered)
                        title(sprintf('%s - z: %d',basenames{z}, z), 'Interpreter', 'none')
                    case 'KeepOriginal'
                        break
                end
            end
        end
    end
    
    %% save all data up to the current fly
    if IntegrateOnly % if part of this datalist has already been processed
        regGeneralData.alignedSliceNumbers(1,find(zs_fly)) = alignedSliceNumbers(1,find(zs_fly)); % it works...
        regGeneralData.rightIpsiSide(1,find(zs_fly)) = rightIpsiSide(1,find(zs_fly)); 
        regGeneralData.allLandmarks(1,find(zs_fly)) = allLandmarks(1,find(zs_fly));
        regGeneralData.anglesRuns(1,find(zs_fly)) = anglesRuns(1,find(zs_fly));
        regGeneralData.anglesStck(1,find(zs_fly)) = anglesStck(1,find(zs_fly));
        regGeneralData.flyfolders(1,find(zs_fly)) = flyfolders(1,find(zs_fly));
        regGeneralData.runfolders(1,find(zs_fly)) = runfolders(1,find(zs_fly));
        regGeneralData.totAng(1,find(zs_fly)) = totAng(1,find(zs_fly));
        regGeneralData.anglesStck_single(1,f) = anglesStck_single(1,f);
        regGeneralData.flyfoldersUnique(1,f) = flyfoldersUnique(1,f);
    else
        save(saveName, 'allStacks', 'anglesRuns', 'anglesStck', 'totAng', 'rightIpsiSide', 'anglesStck_single','flyfoldersUnique','flyfolders','runfolders', 'allLandmarks', 'alignedSliceNumbers', '-v7.3')
        IntegrateOnly = 1;
        regGeneralData = matfile(saveName,'Writable',true); %single z planes to their own stack, general infos
    end
    
end


%% MODULE 1 (86): functional runs alignment
oldAlign86 = matfile('/Users/galileo/Dropbox (HMS)/Data/WINDS_86/anatomyImages/anatomy_alignment_metadata.mat');
alignedSliceNumbers(NaZ_100+1:NaZ) = oldAlign86.alignedSliceNumbers;
cab(1)
flies2process = 6:7; %86 flies
for f = flies2process
    movingPoints=[];
    fixedPoints=[];
    flyfolder = flyfoldersUnique{f};
    [~, flyname] = fileparts(flyfolder);
    flyname = flyname(1:6);
    ZZ = strfind(datalist, flyfolder);
    zs_fly = zeros(size(ZZ));
    for i = 1:length(zs_fly)
        if ~isempty(ZZ{i})
            zs_fly(i) = 1; % indices into datalist relative to fly f
        end
    end
    zs_fly = logical(zs_fly);
    
    %%  get each funcional run's acquisition angle
    for z = find(zs_fly)
        metaimage = aZ{z}.metaim;
        anglesRuns(z) = metaimage.scanRotation;
    end
    
    %% match rotation angle to stack
    acqAnglesStacks(f) = allStacks(f).metaimage.scanRotation;
    anglesStck(zs_fly) = acqAnglesStacks(f) + anglesStck_single(f); %this has been used incorrectly for flies 100.
    totAng = anglesRuns-anglesStck; %this will be updated for each fly.
    
    %% problem 1: align each z plane in its own stack
    % fly118 gets aligned to stack of fly 126
    % stackName = fullfile('/Users/galileo/Dropbox/Data/fly126_PP', 'orto_preprocessed_*.tif');
    % importlist = dir(stackName);
    % stackName2 = fullfile('/Users/galileo/Dropbox/Data/fly126_PP', importlist(1).name);
    % stack = readTiffPP(stackName2, 'double');
    
    stack = allStacks(f).stackG;
    for z = find(zs_fly)
        disp(basenames{z})
        runfolder = runfolders{z};
        
        %% manually determine slice and control points. Save all relevant data for future use
        avgName = fullfile(runfolder, ['AVG_*' basenames{z} '*.tif']);
        importlist = dir(avgName);
        avgName = fullfile(runfolder, importlist(1).name);
        avgrun = readTiffPP(avgName, 'double');
        avgrun(end,:) = [];         % same cropping
        avgrun(:,end-5:end) = [];   % same cropping
        avgrun = imrotate(avgrun, totAng(z));
        avgrun = mat2gray(avgrun);
        avgrun = imadjust(avgrun);
        avgrun = adapthisteq(avgrun);
        
        if alignedSliceNumbers(z) == 0 || isnan(rightIpsiSide(z))
            MIJ.createImage(avgrun);
            MIJ.run('In [+]')
            MIJ.run('In [+]')
            MIJ.run('In [+]')
            %% 
            if alignedSliceNumbers(z) == 0
                MIJ.createImage(stack);
                MIJ.run('In [+]')
                alignedSliceNumbers(z) = input('stack slice number:');
                rightIpsiSide(z) = input('right-left. Right? (1-0):');
            else
                MIJ.createImage(stack(:,:,alignedSliceNumbers(z)));
                MIJ.run('In [+]')
                rightIpsiSide(z) = input('right-left. Right? (1-0):');
            end
            %% 
        end
        MIJ.closeAllWindows();
        target = stack(:,:,alignedSliceNumbers(z)); %so far manually determined         % original
        
        %% try to recover scale and translation within slice
        % manually pick reference points
        %     target = imadjust(target);
        %     target = adapthisteq(target);
        if isempty(movingPoints) || isempty(fixedPoints)
            [movingPoints, fixedPoints] = cpselect(3, 2, avgrun, target, 'Wait', true);
        else
            [movingPoints, fixedPoints] = cpselect(3, 2, avgrun, target, movingPoints, fixedPoints, 'Wait', true);
        end
%         [movingPoints, fixedPoints] = cpselect_PP_zoomHijacked(avgrun, target, movingPoints, fixedPoints, 'Wait', true);
        %
        moving_pts_adj= cpcorr(movingPoints, fixedPoints, avgrun, target);
        tform = fitgeotrans(moving_pts_adj,fixedPoints,'nonreflectivesimilarity');
        tformInv = invert(tform);
        Tinv = tformInv.T;
        ss = Tinv(2,1);
        sc = Tinv(1,1);
        scale_recovered = sqrt(ss*ss + sc*sc)
        theta_recovered = atan2(ss,sc)*180/pi
        Roriginal = imref2d(size(target));
        recovered = imwarp(avgrun,tform,'OutputView',Roriginal);
        
        h_montage=figure; imshowpair(target,recovered,'montage')
        h_composite=figure; imshowpair(target,recovered)
        title(sprintf('%s - z: %d',basenames{z}, z), 'Interpreter', 'none')
        
        while 1
            choice = questdlg('Redo?','repeat?','Yes','No','No');
            switch choice
                case 'No'
                    break
                case 'Yes'
                    [movingPoints, fixedPoints] = cpselect(3, 2, avgrun, target, movingPoints, fixedPoints, 'Wait', true);
                    %
                    allLandmarks(z).movingPoints = movingPoints;        % To be saved
                    allLandmarks(z).fixedPoints = fixedPoints;        % To be saved
                    moving_pts_adj= cpcorr(movingPoints, fixedPoints, avgrun, target);
                    tform = fitgeotrans(moving_pts_adj,fixedPoints,'nonreflectivesimilarity');
                    tformInv = invert(tform);
                    Tinv = tformInv.T;
                    ss = Tinv(2,1);
                    sc = Tinv(1,1);
                    scale_recovered = sqrt(ss*ss + sc*sc)
                    theta_recovered = atan2(ss,sc)*180/pi
                    Roriginal = imref2d(size(target));
                    recovered = imwarp(avgrun,tform,'OutputView',Roriginal);
                    
                    h_montage=figure; imshowpair(target,recovered,'montage')
                    h_composite=figure; imshowpair(target,recovered)
                    title(sprintf('%s - z: %d',basenames{z}, z), 'Interpreter', 'none')
            end
        end
        allLandmarks(z).movingPoints = moving_pts_adj;      % To be saved
        allLandmarks(z).fixedPoints = fixedPoints;        % To be saved
        
        export_fig(h_montage, fullfile(runfolder, ['stackRegistration_' basenames{z} '_montage']), '-png', '-m2.5')
        export_fig(h_composite,fullfile(runfolder, ['stackRegistration_' basenames{z} '_composite']),'-tiff', '-CMYK', '-m1')
        save(fullfile(runfolder, ['stackRegistration_' basenames{z} '.mat']), ...
            'tform', 'recovered', 'target', 'avgName', 'Roriginal')
    end
    %% save all data up to the current fly
    if IntegrateOnly % if part of this datalist has already been processed
        regGeneralData.alignedSliceNumbers(1,find(zs_fly)) = alignedSliceNumbers(1,find(zs_fly)); % it works...
        regGeneralData.rightIpsiSide(1,find(zs_fly)) = rightIpsiSide(1,find(zs_fly)); 
        regGeneralData.allLandmarks(1,find(zs_fly)) = allLandmarks(1,find(zs_fly));
        regGeneralData.anglesRuns(1,find(zs_fly)) = anglesRuns(1,find(zs_fly));
        regGeneralData.anglesStck(1,find(zs_fly)) = anglesStck(1,find(zs_fly));
        regGeneralData.flyfolders(1,find(zs_fly)) = flyfolders(1,find(zs_fly));
        regGeneralData.runfolders(1,find(zs_fly)) = runfolders(1,find(zs_fly));
        regGeneralData.totAng(1,find(zs_fly)) = totAng(1,find(zs_fly));
        regGeneralData.anglesStck_single(1,f) = anglesStck_single(1,f);
        regGeneralData.flyfoldersUnique(1,f) = flyfoldersUnique(1,f);
    else
        save(saveName, 'allStacks', 'anglesRuns', 'anglesStck', 'totAng', 'rightIpsiSide', 'anglesStck_single','flyfoldersUnique','flyfolders','runfolders', 'allLandmarks', 'alignedSliceNumbers', '-v7.3')
        IntegrateOnly = 1;
        regGeneralData = matfile(saveName,'Writable',true); %single z planes to their own stack, general infos
    end
regGeneralData
    
    %% second pass?
    runfolders{zs_fly}
    fprintf('\nZetas this fly:\t %s\n\n', num2str(find(zs_fly)))
    figure(1)
    cab(1)
    for z = find(zs_fly)
        runfolder = runfolders{z};
        target = stack(:,:,alignedSliceNumbers(z)); %so far manually determined         % original
        load(fullfile(runfolder, ['stackRegistration_' basenames{z} '.mat']),'recovered')
        h_composite=figure('WindowStyle', 'normal', 'Position', [1921         -40        1440         805]); imshowpair(target,recovered)
        title(sprintf('%s - z: %d',basenames{z}, z), 'Interpreter', 'none')
        reapplyTransform = input(sprintf('Current z = %d. Reapply transformation from z # (0 for no):  ', z));
        if reapplyTransform ~= 0
            avgName = fullfile(runfolder, ['AVG_*' basenames{z} '*.tif']);
            importlist = dir(avgName);
            avgName = fullfile(runfolder, importlist(1).name);
            avgrun = readTiffPP(avgName, 'double');
            avgrun(end,:) = [];         % same cropping
            avgrun(:,end-5:end) = [];   % same cropping
            avgrun = imrotate(avgrun, totAng(z));
            avgrun = mat2gray(avgrun);
            avgrun = imadjust(avgrun);
            avgrun = adapthisteq(avgrun);
            movingPoints = allLandmarks(reapplyTransform).movingPoints;
            fixedPoints = allLandmarks(reapplyTransform).fixedPoints;
%             [movingPoints, fixedPoints] = cpselect(3, 2, avgrun, target, movingPoints, fixedPoints, 'Wait', true);  
            moving_pts_adj= cpcorr(movingPoints, fixedPoints, avgrun, target);
            tform = fitgeotrans(moving_pts_adj,fixedPoints,'nonreflectivesimilarity');
            tformInv = invert(tform);
            Tinv = tformInv.T;
            ss = Tinv(2,1);
            sc = Tinv(1,1);
            scale_recovered = sqrt(ss*ss + sc*sc)
            theta_recovered = atan2(ss,sc)*180/pi
            Roriginal = imref2d(size(target));
            recovered = imwarp(avgrun,tform,'OutputView',Roriginal);
            
            h_composite=figure; imshowpair(target,recovered)
            title(sprintf('%s - z: %d',basenames{z}, z), 'Interpreter', 'none')
            while 1
                choice = questdlg('Redo?','repeat?','Yes','No', 'KeepOriginal', 'No');
                switch choice
                    case 'No'
                        allLandmarks(z).movingPoints = moving_pts_adj;      % To be saved
                        allLandmarks(z).fixedPoints = fixedPoints;        % To be saved
                        
                        export_fig(h_montage, fullfile(runfolder, ['stackRegistration_' basenames{z} '_montage']), '-png', '-m2.5')
                        export_fig(h_composite,fullfile(runfolder, ['stackRegistration_' basenames{z} '_composite']),'-tiff', '-CMYK', '-m1')
                        save(fullfile(runfolder, ['stackRegistration_' basenames{z} '.mat']), ...
                            'tform', 'recovered', 'target', 'avgName', 'Roriginal')
                        break
                        
                    case 'Yes'
                        [movingPoints, fixedPoints] = cpselect(3, 2, avgrun, target, movingPoints, fixedPoints, 'Wait', true);
                        %
                        allLandmarks(z).movingPoints = movingPoints;        % To be saved
                        allLandmarks(z).fixedPoints = fixedPoints;        % To be saved
                        moving_pts_adj= cpcorr(movingPoints, fixedPoints, avgrun, target);
                        tform = fitgeotrans(moving_pts_adj,fixedPoints,'nonreflectivesimilarity');
                        tformInv = invert(tform);
                        Tinv = tformInv.T;
                        ss = Tinv(2,1);
                        sc = Tinv(1,1);
                        scale_recovered = sqrt(ss*ss + sc*sc)
                        theta_recovered = atan2(ss,sc)*180/pi
                        Roriginal = imref2d(size(target));
                        recovered = imwarp(avgrun,tform,'OutputView',Roriginal);
                        
                        h_montage=figure; imshowpair(target,recovered,'montage')
                        h_composite=figure; imshowpair(target,recovered)
                        title(sprintf('%s - z: %d',basenames{z}, z), 'Interpreter', 'none')
                        
                    case 'KeepOriginal'
                        break
                end
            end
        end
    end
    
    %% save all data up to the current fly
    if IntegrateOnly % if part of this datalist has already been processed
        regGeneralData.alignedSliceNumbers(1,find(zs_fly)) = alignedSliceNumbers(1,find(zs_fly)); % it works...
        regGeneralData.rightIpsiSide(1,find(zs_fly)) = rightIpsiSide(1,find(zs_fly)); 
        regGeneralData.allLandmarks(1,find(zs_fly)) = allLandmarks(1,find(zs_fly));
        regGeneralData.anglesRuns(1,find(zs_fly)) = anglesRuns(1,find(zs_fly));
        regGeneralData.anglesStck(1,find(zs_fly)) = anglesStck(1,find(zs_fly));
        regGeneralData.flyfolders(1,find(zs_fly)) = flyfolders(1,find(zs_fly));
        regGeneralData.runfolders(1,find(zs_fly)) = runfolders(1,find(zs_fly));
        regGeneralData.totAng(1,find(zs_fly)) = totAng(1,find(zs_fly));
        regGeneralData.anglesStck_single(1,f) = anglesStck_single(1,f);
        regGeneralData.flyfoldersUnique(1,f) = flyfoldersUnique(1,f);
    else
        save(saveName, 'allStacks', 'anglesRuns', 'anglesStck', 'totAng', 'rightIpsiSide', 'anglesStck_single','flyfoldersUnique','flyfolders','runfolders', 'allLandmarks', 'alignedSliceNumbers', '-v7.3')
        IntegrateOnly = 1;
        regGeneralData = matfile(saveName,'Writable',true); %single z planes to their own stack, general infos
    end
regGeneralData
end

% %% fix
% save(saveName, 'allStacks', 'anglesRuns', 'anglesStck', 'totAng', 'rightIpsiSide', 'anglesStck_single','flyfoldersUnique','flyfolders','runfolders', 'allLandmarks', 'alignedSliceNumbers', '-v7.3')


%% MODULE 2: align to new_ref_fly119
load('/Users/galileo/Dropbox (HMS)/Data/TDTstacks/newReferenceStack_fly119_cropped.mat'); %reference reference_MP
movingPoints = [];
fixedPoints = [];  % check if no errors for f==1
figure;


%
for f = 1 : length(allStacks)
    cab(1)
    disp(f)
    MP_stack = allStacks(f). maxProjectionG;
    [movingPoints, fixedPoints] = cpselect(1.25,1,MP_stack./max(MP_stack(:)), reference_MP, 'Wait', true);
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
                [movingPoints, fixedPoints] = cpselect(1.25,1,MP_stack./max(MP_stack(:)), reference_MP, movingPoints, fixedPoints, 'Wait', true);
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

%% second pass
% - compare
for f = 1:length(allStacks)
   figure(f), 
   subplot(1,2,1)
   imshowpair(reference_MP, allStacks(f).recovered)
   subplot(1,2,2)
   imshow(allStacks(f).recovered)
end
%
for f = 2
    disp(f)
    MP_stack = allStacks(f). maxProjectionG;
    [movingPoints, fixedPoints] = cpselect(1.25,1,MP_stack./max(MP_stack(:)), reference_MP, allStacks(f).movingPoints, allStacks(f).fixedPoints,'Wait', true);
    tform = fitgeotrans(movingPoints,fixedPoints,'nonreflectivesimilarity');
    tformInv = invert(tform);
    Tinv = tformInv.T;
    ss = Tinv(2,1);
    sc = Tinv(1,1);
    scale_recovered = sqrt(ss*ss + sc*sc)
    theta_recovered = atan2(ss,sc)*180/pi
    Roriginal = imref2d(size(reference_MP));
    recovered = imwarp(MP_stack,tform,'OutputView',Roriginal);
    figure,
    subplot(1,2,1)
    imshowpair(reference_MP, recovered)
    subplot(1,2,2)
    imshow(recovered)
    pause
    pause
    while 1
        choice = questdlg('Redo?','repeat?','Yes','No','No');
        switch choice
            case 'No'
                break
            case 'Yes'
                [movingPoints, fixedPoints] = cpselect(1.25,1,MP_stack./max(MP_stack(:)), reference_MP, movingPoints, fixedPoints, 'Wait', true);
                tform = fitgeotrans(movingPoints,fixedPoints,'nonreflectivesimilarity');
                tformInv = invert(tform);
                Tinv = tformInv.T;
                ss = Tinv(2,1);
                sc = Tinv(1,1);
                scale_recovered = sqrt(ss*ss + sc*sc)
                theta_recovered = atan2(ss,sc)*180/pi
                Roriginal = imref2d(size(reference_MP));
                recovered = imwarp(MP_stack,tform,'OutputView',Roriginal);
                figure, 
                subplot(1,2,1)
                imshowpair(reference_MP, recovered)
                subplot(1,2,2)
                imshow(recovered)
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
allStacks(f).masked_alignedstackR = [];
allStacks(f).aligned_Rboundaries = [];
    
end        
save(saveName, 'allStacks', '-v7.3')


%% plot serialStacks - one row - no rectangles
anatomyfolder = '/Users/galileo/Dropbox (HMS)/Data/WINDS_100/anatomyImages';
allStacks_flyNums = cat(1,allStacks.flyNum);
for f = 2 : length(allStacks)
    disp(f)
    [fig_alStackSer, ax_alStackSer, fPosition] = plotSerialStack_imagesc_oneRowPxls(allStacks(f).alignedstackG,[],gray);
    disp('move figures to the left, then press any key here to continue:')
    pause
    disp('Resuming')
    fig_alStackSer.Position = fPosition;
    sStName = fullfile(anatomyfolder, sprintf('serial_stack2fly119_onerow_fly%d', allStacks_flyNums(f)));
    export_fig(fig_alStackSer, sStName, '-eps');
    close(fig_alStackSer)
end


%% only serial stack with rectangles
anatomyfolder = '/Users/galileo/Dropbox (HMS)/Data/WINDS_100/anatomyImages';

for f = flies2process       % this f was indexed as z before.

    flyfolder = flyfoldersUnique{f};
    [~, flyname] = fileparts(flyfolder);
    flyname = flyname(1:6);
    stack2fly150_info = matfile(fullfile(anatomyfolder, ['stack2fly119ref_' flyname '.mat']));
    
    [fig_alStackSer, ax_alStackSer] = plotSerialStack_imagesc(allStacks(f).alignedstackG,[],gray);
    
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
    stack2fly150_info = matfile(fullfile(flyfolder, ['stack2fly119ref_fly' flyname '.mat']));
    
    stackName = fullfile(flyfolder, 'preprocessed_R_*.tif'); 
    importlist = dir(stackName);
    stackName = fullfile(flyfolder, importlist(1).name);
    redStack = readTiffPP(stackName, 'double');
    % rotate:
    redStackRot = imrotate(redStack, -1*regGeneralData.anglesStck_single(1, f));
    al2fly150_redStack = imwarp(redStackRot,stack2fly150_info.tform,'OutputView',stack2fly150_info.Roriginal); %,'FillValues',[1,1,1]);
    alignStackName = fullfile(flyfolder, ['stack2fly119ref_fly' importlist(1).name]);
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
    sStName = fullfile(flyfolder, ['serial_stack2fly150sStack_' flyname '_RED_StackWithRectangles']);   %afterward moved to the common WINDS_100/anatomyImages folder
    export_fig(fig_alStackSer, sStName, '-eps');
    close(fig_alStackSer)
end

 
%% extra to detect step size in each stack
zfactor = 0.25;
for f = flies2process
    zStepSize_um(f) = abs(allStacks(f).metaimage.zStepSize) * zfactor;    
end
regGeneralData.zStepSize_um = zStepSize_um;



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
