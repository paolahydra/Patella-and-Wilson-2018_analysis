% CNS: all relevant accessory information about single pixel tuning, keeping z planes, and all possible indices
% grouped and saved together


% updated iLATWED mask and included zs for ipsiWED on 12/12
% fixed fly118 post-hoc 12/12



%run after running: 
WEDvsAMMC_pixelwise_distributions_selectZPlanesAndRegions




%% manually select AMMC planes (excluding tdTomato flies, contralateral planes, and planes where just a few pixels should be targeted
%(not actually happening))
% run 
% WEDvsAMMC_pixelwise_distributions_selectZPlanesAndRegions

%%
iCNS.excludeKs = CNS.excludeKs;
iCNS.pixel2keep = pixel2keep;
iCNS.ROIIndices = ROIIndices;
iCNS.dataTable = dataTable;
iCNS.Tparent = Tparent;

iCNS.AMMC.keepZetas = keepZetas_AMMCipsi;
iCNS.AMMC.keepPixelsFromZetas = keepPixels_AMMCipsi;
iCNS.AMMC.pxAMMC = pxAMMC;
iCNS.AMMC.AMMCipsi_Ki = AMMCipsi_Ki;
iCNS.AMMC.AMMCipsi_all = logical(sum(AMMCipsi_Ki,2));

iCNS.WEDi.keepZetas = keepZetas_TONOTipsi;
iCNS.WEDi.keepPixelsFromZetas = keepPixels_TONOTipsi;
iCNS.WEDi.pxWED = pxITONOT;     % pxWED;
iCNS.WEDi.TONOTipsi_Ki = TONOTipsi_Ki;
iCNS.WEDi.TONOTipsi_all = logical(sum(TONOTipsi_Ki,2));

iCNS.WEDb.keepZetas = keepZetas_TONOTboth;
iCNS.WEDb.keepPixelsFromZetas = keepPixels_TONOTboth;
iCNS.WEDb.pxWED = pxWED;
iCNS.WEDb.TONOTboth_Ki = TONOTboth_Ki;
iCNS.WEDb.TONOTboth_all = logical(sum(TONOTboth_Ki,2));

iCNS.WEDc.keepZetas = keepZetas_TONOTcontra;
iCNS.WEDc.keepPixelsFromZetas = keepPixels_TONOTcontra;
iCNS.WEDc.pxWED = pxWED;
iCNS.WEDc.TONOTcontra_Ki = TONOTcontra_Ki;
iCNS.WEDc.TONOTcontra_all = logical(sum(TONOTcontra_Ki,2));

save('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/iCNS.mat', 'iCNS')


%% JONs


%%
iJON.excludeKs = [8,9,10]; % brown , pull, push
iJON.pixel2keep = pixel2keep;
iJON.ROIIndices = ROIIndices;
iJON.dataTable = dataTable;
iJON.Tparent = Tparent;

iJON.PAN.keepZetas = 1:31;
iJON.PAN.keepPixelsFromZetas = keepPixels_AMMCipsi; %temporarily used this name
% load('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/anatomyImages/JON_cropMask.mat')
% iJON.PAN.cropMask = cropMask; %this is in the registered frame

save('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix/iJON.mat', 'iJON')

