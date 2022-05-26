% what do I need to load here?
% not sure I need R
% I need fTC_data for cdfs at single pixel level % load(fullfile(clusterFolder,'fTC_data_01_06secAfterOnset.mat'))
% I need AMMCmasks and a quick separation of the relevant

load_panSpec_downsampled;

% load('/Users/galileo/Dropbox (HMS)/Data/chunks_stimfinal.mat', 'chunks')
% iRun = 14;
% aZ = matfile(datalist{iRun}); %174_run01
% T = aZ.T;
% fastStimMeta = aZ.fastStimMeta;
% clear aZ
% for z = 1:NaZ
%     aZ{z} = matfile(datalist{z});
% end
load('R_matrix_Downsampled_smallWindow.mat', 'pixel2keep', 'iKeep', 'pxKeep')
clusterFolder = '/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix';
cd(clusterFolder)

load('/Users/galileo/Dropbox (HMS)/Data/chunks_stimfinal.mat')
% cd(clusterFolder)
% assert(size(R,1) == sum(pxKeep))
% 
% folder2figure = '/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/figure5n';
% if ~exist(folder2figure, 'dir')
%     mkdir(folder2figure)
% end
% 


%% tuning at single pixel level - (old): peak in a short window % UNA TANTUM - ALREADY DONE
% % output: dff_peak (all possible pixels, not just significant ones);
% [xvar, ia, ic] = unique([chunks.pips.carrierFreqLevels; chunks.tones.carrierFreqLevels(1)]);
% useChunks = [3,5];
% for i = 1:2
%     mClass(i).matfile = matfile(sprintf('FTdata%d.mat', useChunks(i)));
% end %baselines and responses for pips_high and tones_high
% for i = 1 : 2
%     disp(i)
%     baselines = mClass(i).matfile.baselines;
%     responses = mClass(i).matfile.responses;
%     %     if ~isempty(responseWindowPoints)
%     %         responses = responses(:,:,1:responseWindowPoints,:);
%     responses = responses(:,:,3:7,:); %4:9 for 02 07 sec
%     %     end
%     baselines = squeeze(reshape(baselines, [], 1, size(baselines,3), size(baselines,4) ) ) ;
%     responses = squeeze(reshape(responses, [], 1, size(responses,3), size(responses,4) ) ) ;
%     baselines = squeeze(mean(baselines, 2));
%     responses = squeeze(max(responses, [],2));
%     size(baselines)
%     size(responses) %still all possible pixels, by stims
%     if i == 1
%         dff_peak = 100 * (responses-baselines) ./ baselines;
%     else
%         dff_peak_add = 100 * (responses-baselines) ./ baselines;
%         dff_peak = cat(2, dff_peak_add(:,1), dff_peak);
%     end
% end
% clear dff_peak_add baselines responses
% 
% % reduce dff_peak to pxKeep size
% dff_peak = dff_peak(pixel2keep(:),:);
% save('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix/pixelwisePeakFT.mat', 'dff_peak');

%% tuning at single pixel level - averaging in a bigger window to reduce noise %UNA TANTUM - ALREADY DONE
% iRun = 14;
% aZ = matfile(datalist{iRun}); %174_run01
% T = aZ.T;
% ts_dec_pip = T(7).ts_dec;
% i_resp = (find(ts_dec_pip>0.2 & ts_dec_pip<1.3) )'; %4:16
% 
% 
% % output: dff_peak (all possible pixels, not just significant ones);
% [xvar, ia, ic] = unique([chunks.pips.carrierFreqLevels; chunks.tones.carrierFreqLevels(1)]);
% useChunks = [3,5];
% for i = 1:2
%     mClass(i).matfile = matfile(sprintf('FTdata%d.mat', useChunks(i)));
% end %baselines and responses for pips_high and tones_high
% for i = 1 : 2
%     disp(i)
%     baselines = mClass(i).matfile.baselines;
%     responses = mClass(i).matfile.responses;
%     %     if ~isempty(responseWindowPoints)
%     %         responses = responses(:,:,1:responseWindowPoints,:);
%     responses = responses(:,:,i_resp,:); %4:9 for 02 07 sec
%     %     end
%     baselines = squeeze(reshape(baselines, [], 1, size(baselines,3), size(baselines,4) ) ) ;
%     responses = squeeze(reshape(responses, [], 1, size(responses,3), size(responses,4) ) ) ;
%     baselines = squeeze(mean(baselines, 2));
%     responses = squeeze(mean(responses,2));
%     
%     size(baselines)
%     size(responses) %still all possible pixels, by stims
%     if i == 1
%         dff_peak = 100 * (responses-baselines) ./ baselines;
%     else
%         dff_peak_add = 100 * (responses-baselines) ./ baselines;
%         dff_peak = cat(2, dff_peak_add(:,1), dff_peak);
%     end
% end
% clear dff_peak_add baselines responses
% 
% % reduce dff_peak to pxKeep size
% dff_peak = dff_peak(pixel2keep(:),:);
% save('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix/pixelwisePeakFT_mean.mat', 'dff_peak');


%% 171206
load('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix/pixelwisePeakFT_mean.mat', 'dff_peak');
sizeMap = [60, 86];
load('klust_maxClust10.mat'), nK = length(unique(Tparent)); 
clust = remapIndexedPointsToRun_RoisMaps(Tparent, pixel2keep, sizeMap);
ROIs = clust.ROIs;
[xvar, ia, ic] = unique([chunks.pips.carrierFreqLevels; chunks.tones.carrierFreqLevels(1)]);

for k = 1 : size(ROIs,2)
    ROIAllPixelIndices(:,k) = logical(reshape(squeeze(ROIs(:,k,:)), [], 1));
end
ROIIndices = ROIAllPixelIndices(pixel2keep(:),:);
clear ROIAllPixelIndices


flyNumUnique = flyNumUnique(1:6);
tcJON_clusterTuningFly = nan( size(ROIIndices,2), size(dff_peak,2) , length(flyNumUnique) );
tcJON_sumPixels_clusterFly = nan(size(ROIIndices,2),length(flyNumUnique) );
for f = 1:length(flyNumUnique)
    
    zetas = dataTable.fly == flyNumUnique(f);
    builder = false(size(pixel2keep));
    builder(:,zetas) = true;
    idx_thisFly = builder(pixel2keep(:)); %all significant pixels belonging to this fly
    
    indices_K_thisFly = bsxfun(@and, ROIIndices, idx_thisFly);
    
    tcJON_sumPixels_clusterFly(:,f) = sum(indices_K_thisFly);
    for k = 1 : size(ROIs,2)
        tcJON_clusterTuningFly(k,:,f) = nanmean(dff_peak( indices_K_thisFly(:, k) ,:), 1); % OK NOW
    end
end


% % double check the cluster-averaged tuning curves % THIS DOES WORK  %okPP#171206
% for k = 1 : size(ROIs,2)
%     clustAvgTC_JON(k,:) = squeeze(nanmean(tcJON_clusterTuningFly(k,:,:), 3));
% end
% [hfig, hax] = figureTracesI_PP_pixels( size(ROIs,2), 160, 200);
% for k = 1:size(ROIs,2)
%     axes(hax(k,1)); hold on
%     plot(xvar, clustAvgTC_JON(k,:));
%     ylabel([k])
%     hax(k,1).FontSize = 4;
% end 

save('ftcJON_clusterTuningFly.mat', 'tcJON_clusterTuningFly', 'tcJON_sumPixels_clusterFly') %still 10 clusters here


