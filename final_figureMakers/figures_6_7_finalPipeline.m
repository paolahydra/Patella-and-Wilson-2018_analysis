%% 6 zero:
% may want to tweak the ipsiTONOT inclusion criteria, along with a more
% restrictive mask of tehe lateral wedge than currently using.
% use to select a mask

edit MakeUnregistered_ipsiTONOTmasks.m
% saved in:
% save('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/pixelwisePeakFT_mean.mat', 'pxITONOT', '-append')
% also saved 'positionIPSITONOMASK' if you want to reconstruct the ask that
% was used

% % older notes
% % alMaps = load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/19clusters/clusterHeatMaps/alignedData.mat');
% % keepZetas_TONOTipsi = ...
% %     [ 19    20    21    22    24    25    26    16     5     6    12    14    36    41    42];
% % A = mean(alMaps.singleKmaps_px2keep(:,:,keepZetas_TONOTipsi), 3);
% % figure; imshow(A, []);
% % h = imfreehand;
% % iTONOTmask = h.createMask;
% %% check
% figure; imshow(alMaps.singleKmaps_px2keep(:,:,40), [])

%% 6B rev -FINISHED
edit JON_pixelwise_distributions
edit JON_pixelwise_distributions_bestFreq
%revised as for CNS below - ALTHOUGH I HAVE NOT USED A CROP MASK HERE... NO BIG DEAL FOR PAN, but you may fix it
% (find code in JONvsWEDvsAMMC_clustFlyWise_distributions_selectZPlanes.m )

edit WEDvsAMMC_pixelwise_distributions_selectZPlanesAndRegions.m %con extra divisione TONOTboth.
% this file also has the code for reverse transformations!

%FIXED ERROR THAT SCREWED ALL TUNING CURVES!!
% OUTPUT:
% save('ftcJON_clusterTuningFly.mat', 'tcJON_clusterTuningFly', 'tcJON_sumPixels_clusterFly') %still 10 clusters here
%
% save('ftc_AMMCi_TONOTi_TONOTc_clusterTuningFly_mean.mat', 'tcAMMCipsi_clusterTuningFly', 'tcTONOTipsi_clusterTuningFly', 'tcTONOTcontra_clusterTuningFly', 'tcTONOTboth_clusterTuningFly', ...
%     'tcAMMCipsi_sumPixels_clusterFly', 'tcTONOTipsi_sumPixels_clusterFly', 'tcTONOTcontra_sumPixels_clusterFly', 'tcTONOTboth_sumPixels_clusterFly') % all 19 clusters

% put things together and plot:
edit JONvsWEDvsAMMC_clustFlyWise_distributions_selectZPlanes.m

%% 6Bbis
% replot and save pixelwise cdfs by regions to show that WED is more
% uniformily covering frequency space
edit JONvsWEDvsAMMC_pixelwise_distributions_selectZPlanes.m

%% 6A: tonotopy index
edit save_iJON_iCNS_indexList.m %UNA TANTUM, after running WEDvsAMMC_pixelwise_distributions_selectZPlanesAndRegions.m
% --DONE


edit JONvsWEDvsAMMC_clustFlyWise_tonotopyIndex.m % also includes pixelwise PDFs
edit tonotopicAxiMaps.m %currenly in figure folder

edit tonotopicAxiMaps_peakFreq.m    % revision_2018

% data saved and reloaded/plotted from:
% corrTonotopies(1) = load('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix/tonotopyIndex_JONS.mat', 'correlations', 'rotations');
% corrTonotopies(2) = load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/tonotopyIndex_AMMC.mat', 'correlations', 'rotations');
% corrTonotopies(3) = load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/tonotopyIndex_iWED.mat', 'correlations', 'rotations');
% corrTonotopies(4) = load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/tonotopyIndex_cWED.mat', 'correlations', 'rotations');
% export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/tonotopicIndex_selectedZPlanes.eps')


%mostly done.
% still need to do:
% figure out whether we want to exclude a couple of runs from current iWED
% (119_5-6), also in lifetime sparseness calculation

% display (some) maps with relative axis and index - AND CHECK ALL
% ROTATIONS (+-pi)

% strip chart of correlation values (redo)

%% 6C adaptation index, and 6D
edit JONvsWEDvsAMMC_adaptationIndex.m


%% 2 
% winds figure: I had to do something which I do not even remember
% e mo vi avita st citt?