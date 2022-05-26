% what do I need to load here?
% not sure I need R
% I need fTC_data for cdfs at single pixel level % load(fullfile(clusterFolder,'fTC_data_01_06secAfterOnset.mat'))
% I need AMMCmasks and a quick separation of the relevant

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
% assert(size(R,1) == sum(pxKeep))
% 
% folder2figure = '/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/figure5n';
% if ~exist(folder2figure, 'dir')
%     mkdir(folder2figure)
% end
% 
% load('klust_maxClust19.mat'), nK = length(unique(Tparent)); 


%% tuning at single pixel level - (old): peak in a short window
% output: dff_peak (all possible pixels, not just significant ones);
[xvar, ia, ic] = unique([chunks.pips.carrierFreqLevels; chunks.tones.carrierFreqLevels(1)]);
useChunks = [3,5];
for i = 1:2
    mClass(i).matfile = matfile(sprintf('FTdata%d.mat', useChunks(i)));
end %baselines and responses for pips_high and tones_high
for i = 1 : 2
    disp(i)
    baselines = mClass(i).matfile.baselines;
    responses = mClass(i).matfile.responses;
    %     if ~isempty(responseWindowPoints)
    %         responses = responses(:,:,1:responseWindowPoints,:);
    responses = responses(:,:,3:7,:); %4:9 for 02 07 sec
    %     end
    baselines = squeeze(reshape(baselines, [], 1, size(baselines,3), size(baselines,4) ) ) ;
    responses = squeeze(reshape(responses, [], 1, size(responses,3), size(responses,4) ) ) ;
    baselines = squeeze(mean(baselines, 2));
    responses = squeeze(max(responses, [],2));
    size(baselines)
    size(responses) %still all possible pixels, by stims
    if i == 1
        dff_peak = 100 * (responses-baselines) ./ baselines;
    else
        dff_peak_add = 100 * (responses-baselines) ./ baselines;
        dff_peak = cat(2, dff_peak_add(:,1), dff_peak);
    end
end
clear dff_peak_add baselines responses

% reduce dff_peak to pxKeep size
dff_peak = dff_peak(pixel2keep(:),:);
save('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/pixelwisePeakFT.mat', 'dff_peak');

%% tuning at single pixel level - averaging in a bigger window to reduce noise
i_resp = 4:16;
[xvar, ia, ic] = unique([chunks.pips.carrierFreqLevels; chunks.tones.carrierFreqLevels(1)]);
useChunks = [3,5];
for i = 1:2
    mClass(i).matfile = matfile(sprintf('FTdata%d.mat', useChunks(i)));
end %baselines and responses for pips_high and tones_high
for i = 1 : 2
    disp(i)
    baselines = mClass(i).matfile.baselines;
    responses = mClass(i).matfile.responses;
    %     if ~isempty(responseWindowPoints)
    %         responses = responses(:,:,1:responseWindowPoints,:);
    responses = responses(:,:,i_resp,:); %4:9 for 02 07 sec
    %     end
    baselines = squeeze(reshape(baselines, [], 1, size(baselines,3), size(baselines,4) ) ) ;
    responses = squeeze(reshape(responses, [], 1, size(responses,3), size(responses,4) ) ) ;
    baselines = squeeze(mean(baselines, 2));
    responses = squeeze(mean(responses,2));
    
    size(baselines)
    size(responses) %still all possible pixels, by stims
    if i == 1
        dff_peak = 100 * (responses-baselines) ./ baselines;
    else
        dff_peak_add = 100 * (responses-baselines) ./ baselines;
        dff_peak = cat(2, dff_peak_add(:,1), dff_peak);
    end
end
clear dff_peak_add baselines responses

% reduce dff_peak to pxKeep size
dff_peak = dff_peak(pixel2keep(:),:);
save('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/pixelwisePeakFT_mean.mat', 'dff_peak');

%% first separate data between AMMC and WED, possibly in the original (unregistered) pixel space
% prelims
for z = 1:length(datalist)
    [runfolders{z}, ~] = fileparts(datalist{z});
end
regGeneralData = matfile('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/anatomyImages/anatomy_alignment_metadata.mat');
load('/Users/galileo/Dropbox (HMS)/Data/TDTstacks/anatomy_alignment_metadata.mat') %allStacks
allStacks_flyNums = cat(1,allStacks.flyNum);

% % make a simple WED mask - UNA TANTUM
% figure; imshow(allStacks(12).recovered, []);
% h = imrect(gca, h.getPosition);
% alignedWEDmask = createMask(h);
% save('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/datalist_WEDAMMC_piezo.mat', 'alignedWEDmask', '-append')

% %% benchmark forward
% %Take a sigpixel map, its transformations, its transformed map, and check that you can invert all the transormations.
% z = 1;
% flyNum = dataTable.fly(z);
% i_allSt = find(allStacks_flyNums == flyNum);
% omap = reshape(pixel2keep(:,z), sizeMap); figure; imshow(omap,[])
% % retrieve other forward transformations
% regData = matfile(fullfile(runfolders{z}, ['stackRegistration_' basenames{z} '.mat']));
% 
% 
% % apply forward transformations - TEST
% % 1 
% % %convert from downsampled to fullsize
% mapClean = omap;
% mapClean = cat(1, omap, zeros(1,size(mapClean, 2))); %61 x 86
% mapClean = cat(2, mapClean, zeros(size(mapClean, 1), 4)); %61 x 90
% % 2
% mapClean = imresize(mapClean, [86 128]); %should become 86 x 128. The quality sucks, but so be it
% % 3
% % %re-crop
% mapClean(end,:) = [];
% mapClean(:,end-5:end) = [];
% mapClean(mapClean<0.3) = 0;
% mapClean(mapClean>=0.1) = 1;                        %note 1 and not Ki!
% % 4
% mapClean = imrotate(mapClean, regGeneralData.totAng(1,z)); %preliminary transformation
% % 5
% recovered_mapClean = imwarp(mapClean,regData.tform,'OutputView',regData.Roriginal,'FillValues',0); %indexed image
% % 6
% al2fly119_mapClean = imwarp(recovered_mapClean,allStacks(i_allSt).tform,'OutputView',allStacks(i_allSt).Roriginal,'FillValues',0);
% figure; imshow(al2fly119_mapClean,[]);
% 
% 
% %% benchmark inverse
% t0_size = sizeMap;
% t1_size = sizeMap + [1,4]; %after t1 has been applied
% t2_size = [86 128];
% t3_size = [86 128] - [1,6];
% t4_size = size(imrotate(zeros(t3_size), regGeneralData.totAng(1,z))); %need to load it
% t5_size = [700, 700]; %need to load it
% t6_size = allStacks(i_allSt).Roriginal.ImageSize; %need to load it
% 
% % 6 
% inv6 = invert(allStacks(i_allSt).tform);
% inv6_Roriginal = imref2d(t5_size);
% calc_recovered_mapClean = imwarp(al2fly119_mapClean, inv6, 'OutputView', inv6_Roriginal);
% figure; imshowpair(recovered_mapClean,calc_recovered_mapClean); title('inv t6') % OK - this could gets rid of lateral cell bodies already, if I define the mask in t6.
% 
% % 5 
% inv5 = invert(regData.tform);
% inv5_Roriginal = imref2d(t4_size);
% calc_mapClean = imwarp(calc_recovered_mapClean, inv5, 'OutputView', inv5_Roriginal);
% figure; imshowpair(mapClean,calc_mapClean); title('inv t5') % OK - this could gets rid of lateral cell bodies already, if I define the mask in t6. % OK
% 
% % 4 first recontruct the actual rotation
% theta = regGeneralData.totAng(1,z);
% inv4 = affine2d([cosd(theta) sind(theta) 0;...
%     -sind(theta) cosd(theta) 0; 0 0 1]);
% inv4_Roriginal = imref2d(t3_size);
% i4_mapClean = imwarp(calc_mapClean, inv4);
% i4_mapClean = imresize(i4_mapClean, size(i4_mapClean)+1);
% % assume it's SORT-OF-perfectly centered, need to crop the size you need (t3_size) at the center
% crop_pre = ceil((size(i4_mapClean) - t3_size) ./ 2);
% i4_mapClean(1:crop_pre(1), :) = [];
% i4_mapClean(:, 1:crop_pre(2)) = [];
% i4_mapClean(t3_size(1)+1:end,:) = [];
% i4_mapClean(:,t3_size(2)+1:end) = [];
% % figure; imshow(i4_mapClean,[]);
% figure; imshowpair(mapClean,i4_mapClean); title('inv t4') % OK - this could gets rid of lateral cell bodies already, if I define the mask in t6. % COOL
% 
% % 3 - uncrop
% i3_mapClean = zeros(t2_size);
% i3_mapClean(1:t3_size(1),1:t3_size(2)) = i4_mapClean;
% % 2 - downsize
% i2_mapClean = imresize(i3_mapClean, t1_size);
% % 1 - recrop
% i_mapClean = i2_mapClean(1:t0_size(1),1:t0_size(2));
% figure; imshow(i_mapClean,[]);
% figure; imshowpair(omap,i_mapClean); title('full inversion') % OK - this could gets rid of lateral cell bodies already, if I define the mask in t6. % COOL

%% now for all: compute all the relevant masks on the registered stack and invert back.
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
        for z = zetas'
            thisMZ = dataTable.matrixsZ(z) + 11;
            masksregAMMC(:,:,z) = alignedAMMCmasks(:,:,thisMZ); %AMMC
            masksregWED(:,:,z) = alignedWEDmask & ~alignedAMMCmasks(:,:,thisMZ); %WED
            masksOrAMMC(:,:,z) = zeros(size(mAMMC));
            masksOrWED(:,:,z) = ones(size(mWED));
        end
        continue 
    end
    
    for z = zetas'
        disp(z)
        omap = reshape(pixel2keep(:,z), sizeMap); %for sanity check
        regData = matfile(fullfile(runfolders{z}, ['stackRegistration_' basenames{z} '.mat']));
        
        %goal is to inverse-transform masks for AMMC and WED.
        thisMZ = dataTable.matrixsZ(z) + 11;
        masksregAMMC(:,:,z) = alignedAMMCmasks(:,:,thisMZ); %AMMC
        masksregWED(:,:,z) = alignedWEDmask & ~alignedAMMCmasks(:,:,thisMZ); %WED
        %go on inverting the masks
        
        % 6 
        inv6 = invert(allStacks(i_allSt).tform);
        inv6_Roriginal = imref2d(t5_size);
        mAMMC = imwarp(masksregAMMC(:,:,z), inv6, 'OutputView', inv6_Roriginal);
        mWED = imwarp(masksregWED(:,:,z), inv6, 'OutputView', inv6_Roriginal);     
% figure; imshowpair(recovered_mapClean,mAMMC); title('inv t6') % OK - this could gets rid of lateral cell bodies already, if I define the mask in t6.
%correct

        % 5
        t4_size = size(imrotate(zeros(t3_size), regGeneralData.totAng(1,z))); %need to load it
        inv5 = invert(regData.tform);
        inv5_Roriginal = imref2d(t4_size);
        mAMMC = imwarp(mAMMC, inv5, 'OutputView', inv5_Roriginal);
        mWED = imwarp(mWED, inv5, 'OutputView', inv5_Roriginal);
% figure; imshowpair(mapClean,mAMMC); title('inv t5') % OK - this could gets rid of lateral cell bodies already, if I define the mask in t6. % OK
% % correct
           
        % 4 first reconstruct the actual rotation
        theta = regGeneralData.totAng(1,z);
        inv4 = affine2d([cosd(theta) sind(theta) 0;...
            -sind(theta) cosd(theta) 0; 0 0 1]);
        inv4_Roriginal = imref2d(t3_size);
        mAMMC = imwarp(mAMMC, inv4);
        mAMMC = imresize(mAMMC, size(mAMMC)+1);
        mWED = imwarp(mWED, inv4);
        mWED = imresize(mWED, size(mWED)+1);
        
        % assume it's SORT-OF-perfectly centered, need to crop the size you need (t3_size) at the center
        crop_pre = ceil((size(mAMMC) - t3_size) ./ 2);
        mAMMC(1:crop_pre(1), :) = [];
        mAMMC(:, 1:crop_pre(2)) = [];
        mAMMC(t3_size(1)+1:end,:) = [];
        mAMMC(:,t3_size(2)+1:end) = [];
        mWED(1:crop_pre(1), :) = [];
        mWED(:, 1:crop_pre(2)) = [];
        mWED(t3_size(1)+1:end,:) = [];
        mWED(:,t3_size(2)+1:end) = [];        
% figure; imshowpair(mapClean,mAMMC); title('inv t4') % OK - this could gets rid of lateral cell bodies already, if I define the mask in t6. % COOL
% % correct (i.e. -- just as wrong)
    
        % 3 - uncrop
        i4_mapClean = mAMMC;
        mAMMC = zeros(t2_size);
        mAMMC(1:t3_size(1),1:t3_size(2)) = i4_mapClean;
        i4_mapClean = mWED;
        mWED = zeros(t2_size);
        mWED(1:t3_size(1),1:t3_size(2)) = i4_mapClean;
        % 2 - downsize
        mAMMC = imresize(mAMMC, t1_size);
        mWED = imresize(mWED, t1_size);
        % 1 - recrop
        mAMMC = mAMMC(1:t0_size(1),1:t0_size(2));
        mWED = mWED(1:t0_size(1),1:t0_size(2));   
% figure; imshowpair(omap,mAMMC); title('full inversion') % OK - this could gets rid of lateral cell bodies already, if I define the mask in t6. % COOL
% % correct -- it works now (with fly 120)
             
        % finally, re-binarize
        mAMMC(mAMMC<=0.5) = 0;
        mAMMC(mAMMC>0.5) = 1;
        mWED(mWED<=0.5) = 0;
        mWED(mWED>0.5) = 1;
%         figure; imshowpair(omap,mAMMC); title('full inversion')
%         figure; imshowpair(omap,mWED); title('full inversion') % COOL

        % store
        masksOrAMMC(:,:,z) = mAMMC;
        masksOrWED(:,:,z) = mWED;
    end
end
masksOrAMMC = logical(masksOrAMMC);
masksOrWED = logical(masksOrWED);

%% use inverted masks to select pixels and assign them to specific conditions
% intersect pixels2keep and masksOr
pxAMMC = [];
pxWED = [];
for z = 1:length(datalist)
    p2k = pixel2keep(:,z);
    pxAMMC_z = reshape(masksOrAMMC(:,:,z), [],1);
    pxAMMC_z = pxAMMC_z(p2k);
    pxWED_z = reshape(masksOrWED(:,:,z), [],1);
    pxWED_z = pxWED_z(p2k);
    pxAMMC = cat(1, pxAMMC, pxAMMC_z);
    pxWED = cat(1, pxWED, pxWED_z);
end
clear pxAMMC_z pxWED_z
pxAMMC = logical(pxAMMC);
pxWED = logical(pxWED);

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


%% save some stuff
save('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/pixelwisePeakFT.mat', 'pxAMMC', 'pxWED', '-append')


%%
load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/pixelwisePeakFT.mat')


%% recalculate tuning curves in region, clusters, and split by fly (small window peak) 
sizeMap = [60, 86];
load('klust_maxClust19.mat'), nK = length(unique(Tparent)); 
clust = remapIndexedPointsToRun_RoisMaps(Tparent, pixel2keep, sizeMap);
ROIs = clust.ROIs;
[xvar, ia, ic] = unique([chunks.pips.carrierFreqLevels; chunks.tones.carrierFreqLevels(1)]);
useChunks = [3,5];
saveNameFTData = fullfile(clusterFolder, 'FTdata.mat');

for i = 1:2
    mClass(i).matfile = matfile(sprintf('FTdata%d.mat', useChunks(i)));
end %baselines and responses for pips_high and tones_high

for i = 1 : 2
    disp(i)
    baselines = mClass(i).matfile.baselines;
    responses = mClass(i).matfile.responses;
    responses = responses(:,:,6:8,:); %for peak in small onset window, new version possibly slightly less noisy
    %     responses = responses(:,:,i_resp,:); %for averaging in bigger window
    baselines = squeeze(reshape(baselines, [], 1, size(baselines,3), size(baselines,4) ) ) ;
    responses = squeeze(reshape(responses, [], 1, size(responses,3), size(responses,4) ) ) ;
    %     size(baselines)
    %     size(responses)
    baselines = squeeze(nanmean(baselines, 2));
    responses = squeeze(nanmean(responses, 2));
    tcs(i).dff_peak = 100 * (responses-baselines) ./ baselines;
end

%% AllPixelIndices
fliesAllPixelIndices = [];
for f = 1:length(flyNumUnique) % split flies indices\
    zetas = dataTable.fly == flyNumUnique(f);
    builder = false(prod(sizeMap), length(zetas));
    builder(:,zetas) = true;
    fliesAllPixelIndices = cat(2,fliesAllPixelIndices, builder(:));
end %ok
fliesAllPixelIndices = logical(fliesAllPixelIndices);

% also do a ROIAllPixelIndices
for k = 1 : size(ROIs,2)
    ROIAllPixelIndices(:,k) = logical(reshape(squeeze(ROIs(:,k,:)), [], 1));
end

% just concatenate the 4Hz point once and for all
tctp = cat(2, tcs(2).dff_peak(:,1), tcs(1).dff_peak);

%% bring to pxKeep-reference all indices
ROIIndices = ROIAllPixelIndices(pixel2keep(:),:);
fliesIndices = fliesAllPixelIndices(pixel2keep(:),:);
tctp = tctp(pixel2keep(:),:);

load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/pixelwisePeakFT.mat', 'pxAMMC', 'pxWED')

% AMMC
for f = 1 : length(flyNumUnique)
    for k = 1 : size(ROIs,2)
        idx = fliesIndices(:,f) & ROIIndices(:,k) & pxAMMC;
        tcAMMC_clusterTuningFly(k,:,f) = nanmean(tctp(idx,:),1);
    end
end

% WED
for f = 1 : length(flyNumUnique)
    for k = 1 : size(ROIs,2)
        idx = fliesIndices(:,f) & ROIIndices(:,k) & pxWED;
        tcWED_clusterTuningFly(k,:,f) = nanmean(tctp(idx,:),1);
    end
end

% check
colors = distinguishable_colors(11);
[hfig, hax] = figureTracesI_PP_pixels( 19, 160*ones(1,2), 160 );
for k = 1:19
    axes(hax(k,1)); hold on
    for f = 1:11
        plot(xvar, squeeze(tcAMMC_clusterTuningFly(k,:,f)), 'Color', colors(f,:));
    end
end

for k = 1:19
    axes(hax(k,2)); hold on
    for f = 1:11
        plot(xvar, squeeze(tcWED_clusterTuningFly(k,:,f)), 'Color', colors(f,:));
    end
end

% save it
save('ftc_AMMC_WED_clusterTuningFly.mat', 'tcAMMC_clusterTuningFly', 'tcWED_clusterTuningFly') %still 10 clusters here





%% in sorted subsets of pixels, compute cdfs and relative parameters 
onset = 4;
offset = 600; %cropping before, makes things worse
cutoff = 0.5;
x2 = logspace(log10(onset), log10(offset), 2000);

yvar2AMMC = interp1(xvar, dff_peak(pxAMMC,:)', x2 );
yvar2WED = interp1(xvar, dff_peak(pxWED,:)', x2 );

meanAMMC = mean(yvar2AMMC,2);
meanWED = mean(yvar2WED,2);
figure; hold on
plot(x2, meanAMMC)
plot(x2, meanWED)


cdfsAMMC = bsxfun(@rdivide, cumsum(yvar2AMMC), max(cumsum(yvar2AMMC)) );
% figure; hold on
% plot(x2, mean(cdfsAMMC,2));

cdfsWED = bsxfun(@rdivide, cumsum(yvar2WED), max(cumsum(yvar2WED)) );
% figure; hold on
% plot(x2, mean(cdfsWED,2));

% The means do not really mean anything. Not all the pixels even respond to
% frequency.



%% find best frequency (half maximum point, pixel wise)
% figure; hold on; plot(x2, cdfsWED(:,12000:12010))
%simply:
i_cdfWED = cdfsWED < 0.5;
s_cdfWED = sum(i_cdfWED);
[BWED, IWED] = sort(s_cdfWED);
figure; subplot(2,1,1)
imshow(i_cdfWED(:,IWED));
ax = gca;
ax.YAxis.Visible = 'on';
ax.YDir = 'normal';

i_cdfAMMC = cdfsAMMC < 0.5;
s_cdfAMMC = sum(i_cdfAMMC);
[BAMMC, IAMMC] = sort(s_cdfAMMC);
subplot(2,1,2)
imshow(i_cdfAMMC(:,IAMMC));
ax = gca;
ax.YAxis.Visible = 'on';
ax.YDir = 'normal';

x2Zero = cat(2, 0, x2); 
%%
figure; hold on
line(1:length(BAMMC), x2Zero(BAMMC+1),'Color', 'r'); 
ax1 = gca;
ax1.XLim(2) = length(BAMMC);
ax1.YLim = [x2(1), x2Zero(end)];
ax1.XColor = 'r';
ax1.YColor = 'r';
ax1.XLabel.String = 'AMMC pixels';
ax1.YLabel.String = 'half maximum cdf';
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');

line(1:length(BWED), x2Zero(BWED+1),'Color', 'k'); 
ax2.XLim(2) = length(BWED);
ylim([x2(1), x2Zero(end)]);
ax2.XLabel.String = 'WED pixels';
ax2.YGrid = 'on';%%
ax1.YScale = 'log';
ax2.YScale = 'log';

% export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/preferredDistrib_WEDvsAMMC.eps')

%% CHANGE TO MAKE IT AS #PIXELS AS A FUNCTION OF FREQUENCY (AND NOT THE OTHER WAY AROUND)
figure; hold on
line(x2Zero(BAMMC+1), 1:length(BAMMC),'Color', 'r'); 
ax1 = gca;
ax1.YLim(2) = length(BAMMC);
ax1.XLim = [x2(1), x2Zero(end)];
ax1.XColor = 'r';
ax1.YColor = 'r';
ax1.YLabel.String = 'AMMC pixels';
ax1.XLabel.String = 'half maximum cdf (best freq)';
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');

line(x2Zero(BWED+1), 1:length(BWED),'Color', 'k'); 
ax2.YLim(2) = length(BWED);
xlim([x2(1), x2Zero(end)]);
ax2.YLabel.String = 'WED pixels';
ax2.XGrid = 'on';%%
ax1.XScale = 'log';
ax2.XScale = 'log';
ax1.YScale = 'linear';
ax2.YScale = 'linear';

% export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/preferredDistrib_WEDvsAMMC.eps')












