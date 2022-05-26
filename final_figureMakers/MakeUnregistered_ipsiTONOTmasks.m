%%
load_wedAmmc_downsampled;

% load('/Users/galileo/Dropbox (HMS)/Data/chunks_stimfinal.mat', 'chunks')
% iRun = 14;
% aZ = matfile(datalist{iRun}); %174_run01
% T = aZ.T;
% fastStimMeta = aZ.fastStimMeta;
% clear aZ
% for z = 1:NaZ
%     aZ{z} = matfile(datalist{z});
% end
load(fullfile(Folder2Save, 'R_matrix_Downsampled_wedammc.mat'), 'pxKeep')
load(fullfile(Folder2Save, 'R_matrix_sortedTones_Downsampled_smallWindow.mat'), 'R', 'pixel2keep', 'iKeep')
R = R(pxKeep,:);

clusterFolder = fullfile(Folder2Save, '/rawDataLinkage/ward_ZSCpix');
cd(clusterFolder)

%%
alMaps = load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/19clusters/clusterHeatMaps/alignedData.mat');

keepZetas_TONOTipsi = ...
    [ 19    20    21    22    24    25    26    16     5     6    12    14    36    41    42];
A = mean(alMaps.singleKmaps_px2keep(:,:,keepZetas_TONOTipsi), 3);
figure; imshow(A, []);
h = imfreehand;
iTONOTmask = h.createMask;



%% make iTONOT on the registered stack and invert back.

% prelims
for z = 1:length(datalist)
    [runfolders{z}, ~] = fileparts(datalist{z});
end
regGeneralData = matfile('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/anatomyImages/anatomy_alignment_metadata.mat');
load('/Users/galileo/Dropbox (HMS)/Data/TDTstacks/anatomy_alignment_metadata.mat') %allStacks
allStacks_flyNums = cat(1,allStacks.flyNum);



t0_size = sizeMap;
t1_size = sizeMap + [1,4]; %after t1 has been applied
t2_size = [86 128];
t3_size = [86 128] - [1,6];
t5_size = [700, 700];
t6_size = [486, 589]; 

for f = 1 : NaF
    flyNum = flyNumUnique(f);
    i_allSt = find(allStacks_flyNums == flyNum);
    flyfolder = flyfoldersUnique{f};
    zetas = find(dataTable.fly == flyNum);
    
    if f == 5 %118, skip : it's all WED
        
        warning('this code is actually ewrong, and has been fixed post hoc -- not here')
        for z = zetas'
            thisMZ = dataTable.matrixsZ(z) + 11;
            masksregIPSITONOT(:,:,z) = iTONOTmask & alignedWEDmask & ~alignedAMMCmasks(:,:,thisMZ); %WED
        end
        continue 
        
        
    end
    
    for z = zetas'
        disp(z)
        omap = reshape(pixel2keep(:,z), sizeMap); %for sanity check
        regData = matfile(fullfile(runfolders{z}, ['stackRegistration_' basenames{z} '.mat']));
        
        %goal is to inverse-transform masks for AMMC and WED.
        thisMZ = dataTable.matrixsZ(z) + 11;
        masksregIPSITONOT(:,:,z) = iTONOTmask & alignedWEDmask & ~alignedAMMCmasks(:,:,thisMZ); %WED
        %go on inverting the masks
        
        % 6 
        inv6 = invert(allStacks(i_allSt).tform);
        inv6_Roriginal = imref2d(t5_size);
%         mAMMC = imwarp(masksregAMMC(:,:,z), inv6, 'OutputView', inv6_Roriginal);
        mWED = imwarp(masksregIPSITONOT(:,:,z), inv6, 'OutputView', inv6_Roriginal);     
% figure; imshowpair(recovered_mapClean,mAMMC); title('inv t6') % OK - this could gets rid of lateral cell bodies already, if I define the mask in t6.
%correct

        % 5
        t4_size = size(imrotate(zeros(t3_size), regGeneralData.totAng(1,z))); %need to load it
        inv5 = invert(regData.tform);
        inv5_Roriginal = imref2d(t4_size);
%         mAMMC = imwarp(mAMMC, inv5, 'OutputView', inv5_Roriginal);
        mWED = imwarp(mWED, inv5, 'OutputView', inv5_Roriginal);
% figure; imshowpair(mapClean,mAMMC); title('inv t5') % OK - this could gets rid of lateral cell bodies already, if I define the mask in t6. % OK
% % correct
           
        % 4 first reconstruct the actual rotation
        theta = regGeneralData.totAng(1,z);
        inv4 = affine2d([cosd(theta) sind(theta) 0;...
            -sind(theta) cosd(theta) 0; 0 0 1]);
        inv4_Roriginal = imref2d(t3_size);
%         mAMMC = imwarp(mAMMC, inv4);
%         mAMMC = imresize(mAMMC, size(mAMMC)+1);
        mWED = imwarp(mWED, inv4);
        mWED = imresize(mWED, size(mWED)+1);
        
        % assume it's SORT-OF-perfectly centered, need to crop the size you need (t3_size) at the center
        crop_pre = ceil((size(mWED) - t3_size) ./ 2);
%         mAMMC(1:crop_pre(1), :) = [];
%         mAMMC(:, 1:crop_pre(2)) = [];
%         mAMMC(t3_size(1)+1:end,:) = [];
%         mAMMC(:,t3_size(2)+1:end) = [];
        mWED(1:crop_pre(1), :) = [];
        mWED(:, 1:crop_pre(2)) = [];
        mWED(t3_size(1)+1:end,:) = [];
        mWED(:,t3_size(2)+1:end) = [];        
% figure; imshowpair(mapClean,mAMMC); title('inv t4') % OK - this could gets rid of lateral cell bodies already, if I define the mask in t6. % COOL
% % correct (i.e. -- just as wrong)
    
        % 3 - uncrop
%         i4_mapClean = mAMMC;
%         mAMMC = zeros(t2_size);
%         mAMMC(1:t3_size(1),1:t3_size(2)) = i4_mapClean;
        i4_mapClean = mWED;
        mWED = zeros(t2_size);
        mWED(1:t3_size(1),1:t3_size(2)) = i4_mapClean;
        % 2 - downsize
%         mAMMC = imresize(mAMMC, t1_size);
        mWED = imresize(mWED, t1_size);
        % 1 - recrop
%         mAMMC = mAMMC(1:t0_size(1),1:t0_size(2));
        mWED = mWED(1:t0_size(1),1:t0_size(2));   
% figure; imshowpair(omap,mAMMC); title('full inversion') % OK - this could gets rid of lateral cell bodies already, if I define the mask in t6. % COOL
% % correct -- it works now (with fly 120)
             
        % finally, re-binarize
%         mAMMC(mAMMC<=0.5) = 0;
%         mAMMC(mAMMC>0.5) = 1;
        mWED(mWED<=0.5) = 0;
        mWED(mWED>0.5) = 1;
%         figure; imshowpair(omap,mAMMC); title('full inversion')
%         figure; imshowpair(omap,mWED); title('full inversion') % COOL

        % store
%         masksOrAMMC(:,:,z) = mAMMC;
        masksOrWED(:,:,z) = mWED;
    end
end
% masksOrAMMC = logical(masksOrAMMC);
masksOrWED = logical(masksOrWED);
masksOrITONOT = masksOrWED;

%% use inverted masks to select pixels and assign them to specific conditions
% intersect pixels2keep and masksOr
% pxAMMC = [];
pxITONOT = [];
for z = 1:length(datalist)
    p2k = pixel2keep(:,z);
%     pxAMMC_z = reshape(masksOrAMMC(:,:,z), [],1);
%     pxAMMC_z = pxAMMC_z(p2k);
    pxITONOT_z = reshape(masksOrITONOT(:,:,z), [],1);
    pxITONOT_z = pxITONOT_z(p2k);
%     pxAMMC = cat(1, pxAMMC, pxAMMC_z);
    pxITONOT = cat(1, pxITONOT, pxITONOT_z);
end
clear pxITONOT_z
% pxAMMC = logical(pxAMMC);
pxITONOT = logical(pxITONOT);

% %% double-check
% figure
% for z = 1:49
%     figure;
%     omap = reshape(pixel2keep(:,z), sizeMap); %for sanity check
%     imshowpair(omap,masksOrAMMC(:,:,z)); title(sprintf('AMMC %d',z))
%     pause
%     close
% end
% now YES


% save some stuff
save('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/pixelwisePeakFT.mat', 'pxITONOT', '-append')
save('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/pixelwisePeakFT_mean.mat', 'pxITONOT', '-append')
