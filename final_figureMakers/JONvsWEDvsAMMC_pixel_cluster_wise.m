% here pixel wise FT and separation between AMMC and WED has already been
% calculated.
% all tcs are log interpolated in 200 points here

% high SNR pixels only (percentage within each cluster, from selected runs/flies)

%% load data and relevant indices
load('/Users/galileo/Dropbox (HMS)/Data/chunks_stimfinal.mat')

% load JON
JON = load('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix/pixelwisePeakFT_mean.mat');%NEW - mean of 13 points, rather than peak
% load AMMC + WED
CNS = load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/pixelwisePeakFT_mean.mat');%NEW - mean of 13 points, rather than peak
% CNS.datalist = load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/datalist_WEDAMMC_piezo.mat');

CNS.clusterFolder = '/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix';
JON.clusterFolder = '/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix';

CNS.clusters = load(fullfile(CNS.clusterFolder, 'klust_maxClust19.mat'));
JON.clusters = load(fullfile(JON.clusterFolder, 'klust_maxClust10.mat'));
JON.excludeKs = [8,9,10]; % brown , pull, push
CNS.excludeKs = [1,2]; % pull, push
JON.includedK_i = ~ismember(JON.clusters.Tparent, JON.excludeKs);
CNS.includedK_i = ~ismember(CNS.clusters.Tparent, CNS.excludeKs);

% (log) interpolation parameters here
[xvar, ~, ~] = unique([chunks.pips.carrierFreqLevels; chunks.tones.carrierFreqLevels(1)]);
onset = 4;
offset = 600; %cropping before, makes things worse
cutoff = 0.5;
x2 = logspace(log10(onset), log10(offset), 200);
x2Zero = cat(2, 0, x2); 

%% define z inclusion criteria (not crucial here, just a reminder)
% tdtFlies = unique(CNS.dataTable.fly(CNS.datalist.dataTable.tdTomato==1));
% allCNSflies = unique(CNS.dataTable.fly, 'stable');
% keepCNSflies = ~(ismember(allCNSflies, tdtFlies));
%% manually select AMMC planes (excluding tdTomato flies, contralateral planes, and planes where just a few pixels should be targeted (not actually happening))
keepZetas_AMMCipsi = [ 1     2     4     5     7     8     9    10    11    13    14    35    37    38    39    40    41];
keepPixels_AMMCipsi = false(size(CNS.pixel2keep));
keepPixels_AMMCipsi(:,keepZetas_AMMCipsi) = true;
keepPixels_AMMCipsi = keepPixels_AMMCipsi(:);
keepPixels_AMMCipsi = keepPixels_AMMCipsi(CNS.pixel2keep(:)); %OK

%% manually select ipsi WED planes (excl: contralateral planes, ventral planes Z<=21, and planes where just a few pixels should be targeted) - incl tdT+
keepZetas_WEDipsi = [2     5     6    11    12    14    15    16    20    21    22    25    26    36    41    42 ];
keepPixels_WEDipsi = false(size(CNS.pixel2keep));
keepPixels_WEDipsi(:,keepZetas_WEDipsi) = true;
keepPixels_WEDipsi = keepPixels_WEDipsi(:);
keepPixels_WEDipsi = keepPixels_WEDipsi(CNS.pixel2keep(:)); %OK

%% manually select ipsi WED planes (only include planes where ipsi tonotopy can be separated) - incl tdT+
keepZetas_TONOTipsi = [19    20    21    22    24    25    26    16     5     6    12    14    36    41    42 ];
keepPixels_TONOTipsi = false(size(CNS.pixel2keep));
keepPixels_TONOTipsi(:,keepZetas_TONOTipsi) = true;
keepPixels_TONOTipsi = keepPixels_TONOTipsi(:);
keepPixels_TONOTipsi = keepPixels_TONOTipsi(CNS.pixel2keep(:)); %OK

%% manually select contra WED planes (only include planes where ipsi tonotopy can be separated) - incl tdT+
keepZetas_TONOTcontra = [28    27    30     3    43    44    46    47    48    49 ];
keepPixels_TONOTcontra = false(size(CNS.pixel2keep));
keepPixels_TONOTcontra(:,keepZetas_TONOTcontra) = true;
keepPixels_TONOTcontra = keepPixels_TONOTcontra(:);
keepPixels_TONOTcontra = keepPixels_TONOTcontra(CNS.pixel2keep(:)); %OK

%% keep JONs??


%% SNR CNS - load
load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/Downsampled/R_matrix_Downsampled_wedammc.mat', 'pxKeep')
load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/Downsampled/R_matrix_sortedTones_Downsampled_smallWindow.mat', 'R');
CNS.R = R(pxKeep,:);

saveName = fullfile('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/Downsampled/', 'R_chunks.mat');
load(saveName, 'R_chirps')
baseline_idxs = [5:14, 172:181];
percentileESClude = 33;

clear R pxKeep saveName

% separate by clusters, for all clusters
for k = 1:length(CNS.clusters.klust)
    CNS_Ki(:,k) = CNS.clusters.Tparent == k;
end




%% AMMCipsi
AMMCipsi_Ki = bsxfun(@and, CNS_Ki, (CNS.pxAMMC & keepPixels_AMMCipsi) );
% remove the excluded clusters (not freq-responsive)
AMMCipsi_Ki(:,CNS.excludeKs) = [];

% calculate pixel-wise SNRs
for k = 1:size(AMMCipsi_Ki,2)   
    % SNR 1: peak DF/F
    Rsel = CNS.R(AMMCipsi_Ki(:,k),:);
    clust(k).snrValues(:,1) = max(Rsel,[],2);
    
    % SNR 2: correlation with cluster mean
    clustMean = (mean(Rsel))';
    clust(k).snrValues(:,2) = corr(Rsel', clustMean);
    
    % SNR 3: score by snr and select. Only consider chirps for now
    Rsel = R_chirps(AMMCipsi_Ki(:,k),:);
    clust(k).snrValues(:,3) = std(Rsel, [], 2) ./ std(Rsel(:,baseline_idxs), [], 2);
end


AMMCipsi_Ki_highSNR = false(size(AMMCipsi_Ki));
for k = 1:size(AMMCipsi_Ki,2)
    excl_I3 = clust(k).snrValues(:,3) < prctile( clust(k).snrValues(:,3), percentileESClude);
    excl_I2 = clust(k).snrValues(:,2) < prctile( clust(k).snrValues(:,2), percentileESClude);
    excl_I1 = clust(k).snrValues(:,1) < prctile( clust(k).snrValues(:,1), percentileESClude);
    % sum((excl_I2 | excl_I3) & excl_I1)
    lowSNRpixels = (excl_I2 | excl_I3) & excl_I1;
    indices = find(AMMCipsi_Ki(:,k));
    AMMCipsi_Ki_highSNR(indices(~lowSNRpixels), k) = 1;
end
clear clust

keepPixels_AMMCipsi_highSNR = logical(sum(AMMCipsi_Ki_highSNR,2));

%% TONOTipsi
TONOTipsi_Ki = bsxfun(@and, CNS_Ki, (CNS.pxWED & keepPixels_TONOTipsi) );
% remove the excluded clusters (not freq-responsive)
TONOTipsi_Ki(:,CNS.excludeKs) = [];

% calculate pixel-wise SNRs
for k = 1:size(TONOTipsi_Ki,2)   
    % SNR 1: peak DF/F
    Rsel = CNS.R(TONOTipsi_Ki(:,k),:);
    clust(k).snrValues(:,1) = max(Rsel,[],2);
    
    % SNR 2: correlation with cluster mean
    clustMean = (mean(Rsel))';
    clust(k).snrValues(:,2) = corr(Rsel', clustMean);
    
    % SNR 3: score by snr and select. Only consider chirps for now
    Rsel = R_chirps(TONOTipsi_Ki(:,k),:);
    clust(k).snrValues(:,3) = std(Rsel, [], 2) ./ std(Rsel(:,baseline_idxs), [], 2);
end


TONOTipsi_Ki_highSNR = false(size(TONOTipsi_Ki));
for k = 1:size(TONOTipsi_Ki,2)
    excl_I3 = clust(k).snrValues(:,3) < prctile( clust(k).snrValues(:,3), percentileESClude);
    excl_I2 = clust(k).snrValues(:,2) < prctile( clust(k).snrValues(:,2), percentileESClude);
    excl_I1 = clust(k).snrValues(:,1) < prctile( clust(k).snrValues(:,1), percentileESClude);
    % sum((excl_I2 | excl_I3) & excl_I1)
    lowSNRpixels = (excl_I2 | excl_I3) & excl_I1;
    indices = find(TONOTipsi_Ki(:,k));
    TONOTipsi_Ki_highSNR(indices(~lowSNRpixels), k) = 1;
end
clear clust

keepPixels_TONOTipsi_highSNR = logical(sum(TONOTipsi_Ki_highSNR,2));

%% TONOTcontra
TONOTcontra_Ki = bsxfun(@and, CNS_Ki, (CNS.pxWED & keepPixels_TONOTcontra) );
% remove the excluded clusters (not freq-responsive)
TONOTcontra_Ki(:,CNS.excludeKs) = [];

% calculate pixel-wise SNRs
for k = find(sum(TONOTcontra_Ki))
    % SNR 1: peak DF/F
    Rsel = CNS.R(TONOTcontra_Ki(:,k),:);
    clust(k).snrValues(:,1) = max(Rsel,[],2);
    
    % SNR 2: correlation with cluster mean
    clustMean = (mean(Rsel))';
    clust(k).snrValues(:,2) = corr(Rsel', clustMean);
    
    % SNR 3: score by snr and select. Only consider chirps for now
    Rsel = R_chirps(TONOTcontra_Ki(:,k),:);
    clust(k).snrValues(:,3) = std(Rsel, [], 2) ./ std(Rsel(:,baseline_idxs), [], 2);
end


TONOTcontra_Ki_highSNR = false(size(TONOTcontra_Ki));
for k = find(sum(TONOTcontra_Ki))
    excl_I3 = clust(k).snrValues(:,3) < prctile( clust(k).snrValues(:,3), percentileESClude);
    excl_I2 = clust(k).snrValues(:,2) < prctile( clust(k).snrValues(:,2), percentileESClude);
    excl_I1 = clust(k).snrValues(:,1) < prctile( clust(k).snrValues(:,1), percentileESClude);
    % sum((excl_I2 | excl_I3) & excl_I1)
    lowSNRpixels = (excl_I2 | excl_I3) & excl_I1;
    indices = find(TONOTcontra_Ki(:,k));
    TONOTcontra_Ki_highSNR(indices(~lowSNRpixels), k) = 1;
end
clear clust

keepPixels_TONOTcontra_highSNR = logical(sum(TONOTcontra_Ki_highSNR,2));

%% JONs freq-responsive pixels cdfs
tc.JON = interp1(xvar, JON.dff_peak(JON.includedK_i,:)', x2 );
cdf.JON = bsxfun(@rdivide, cumsum(tc.JON), max(cumsum(tc.JON)) );

%% AMMC freq-responsive pixels cdfs
tc.AMMC = interp1(xvar, CNS.dff_peak(keepPixels_AMMCipsi_highSNR,:)', x2 );
cdf.AMMC = bsxfun(@rdivide, cumsum(tc.AMMC), max(cumsum(tc.AMMC)) );

%% WED freq-responsive pixels cdfs (ipsi)
tc.WED = interp1(xvar, CNS.dff_peak(CNS.includedK_i & CNS.pxWED & keepPixels_WEDipsi,:)', x2 );
cdf.WED = bsxfun(@rdivide, cumsum(tc.WED), max(cumsum(tc.WED)) );

%% ipsi TONOTOPY freq-responsive pixels cdfs
tc.TONOTipsi = interp1(xvar, CNS.dff_peak(keepPixels_TONOTipsi_highSNR,:)', x2 );
cdf.TONOTipsi = bsxfun(@rdivide, cumsum(tc.TONOTipsi), max(cumsum(tc.TONOTipsi)) );

%% contra TONOTOPY freq-responsive pixels cdfs
tc.TONOTcontra = interp1(xvar, CNS.dff_peak(keepPixels_TONOTcontra_highSNR,:)', x2 );
cdf.TONOTcontra = bsxfun(@rdivide, cumsum(tc.TONOTcontra), max(cumsum(tc.TONOTcontra)) );

%% find best frequency (half maximum point in the cdf, pixel wise)
[cdf.half_i_WED, ~] = sort(sum(cdf.WED < 0.5));
[cdf.half_i_TONOTipsi, ~] = sort(sum(cdf.TONOTipsi < 0.5));
[cdf.half_i_TONOTcontra, ~] = sort(sum(cdf.TONOTcontra < 0.5));
[cdf.half_i_AMMC, ~] = sort(sum(cdf.AMMC < 0.5));
[cdf.half_i_JON, ~] = sort(sum(cdf.JON < 0.5));


%% plot #PIXELS AS A FUNCTION OF FREQUENCY (AND NOT THE OTHER WAY AROUND)
cols = brewermap(5, 'Dark2');
figure; hold on
setN = 1; %JONS;
l(setN) = line(x2Zero(cdf.half_i_JON+1), 1:length(cdf.half_i_JON),'Color', cols(setN, :), 'LineWidth', 2); pause;

ax(setN) = gca;
ax(setN).YLim(2) = length(cdf.half_i_JON);
ax(setN).XLim = [x2(1), x2Zero(end)];
ax(setN).XColor = cols(setN, :);
ax(setN).YColor = cols(setN, :);
ax(setN).YTick = [0, 1e4];
ax(setN).YMinorTick	 = 'on';
ax(setN).YAxis.MinorTickValues = round(length(cdf.half_i_JON)/2);
ax(setN).YMinorGrid = 'on';
ax(setN).XGrid = 'on';
ax(setN).YLabel.String = 'pixels sorted by best frequency';
ax(setN).XLabel.String = 'best frequency (half maximum cdf)';

pos = ax(setN).Position; % position of first axes

setN = 2;
ax(setN) = axes('Position', pos,...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none');

l(setN) = line(x2Zero(cdf.half_i_AMMC+1), 1:length(cdf.half_i_AMMC),'Color', cols(setN, :), 'LineWidth', 2); pause;

ax(setN).XAxis.Visible = 'off';
ax(setN).YLim(2) = length(cdf.half_i_AMMC);
ax(setN).XLim = [x2(1), x2Zero(end)];
ax(setN).XColor = cols(setN, :);
ax(setN).YColor = cols(setN, :);
ax(setN).YTick = [0, 1e4];


setN = 3;
ax(setN) = axes('Position', pos,...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none');

l(setN) = line(x2Zero(cdf.half_i_WED+1), 1:length(cdf.half_i_WED),'Color', cols(setN, :), 'LineWidth', 2); pause;

ax(setN).XAxis.Visible = 'off';
ax(setN).YLim(2) = length(cdf.half_i_WED);
ax(setN).XLim = [x2(1), x2Zero(end)];
ax(setN).XColor = cols(setN, :);
ax(setN).YColor = cols(setN, :);
ax(setN).YTick = [0, 1e4];


setN = 4;
ax(setN) = axes('Position', pos,...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none');

l(setN) = line(x2Zero(cdf.half_i_TONOTipsi+1), 1:length(cdf.half_i_TONOTipsi),'Color', cols(setN, :), 'LineWidth', 2); pause;

ax(setN).XAxis.Visible = 'off';
ax(setN).YLim(2) = length(cdf.half_i_TONOTipsi);
ax(setN).XLim = [x2(1), x2Zero(end)];
ax(setN).XColor = cols(setN, :);
ax(setN).YColor = cols(setN, :);
ax(setN).YTick = [0, 1e4];


setN = 5;
ax(setN) = axes('Position', pos,...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none');

l(setN) = line(x2Zero(cdf.half_i_TONOTcontra+1), 1:length(cdf.half_i_TONOTcontra),'Color', cols(setN, :), 'LineWidth', 2); pause;

ax(setN).XAxis.Visible = 'off';
ax(setN).YLim(2) = length(cdf.half_i_TONOTcontra);
ax(setN).XLim = [x2(1), x2Zero(end)];
ax(setN).XColor = cols(setN, :);
ax(setN).YColor = cols(setN, :);
ax(setN).YTick = [0, 1e4];


for i = 1:5
    ax(i).XScale = 'log';
end
legend(l, 'JON', 'AMMC', 'WEDipsi', 'TONipsi', 'TONcontra')
legend boxoff




% export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/preferredDistrib_JON_AMMC_WEDi_TONi_TONc_mean.eps')


%% find spread of tuning curve (sd of fitted gaussian over log spaced freq tuning curves)

%% RmeanRmax (tcs are already interpolated - 2000) - normalize between 0 and 1 first
figure; hold on
title('Rmean distributions after normalization')

ttc = tc.JON;
ttc = bsxfun(@minus, ttc, min(ttc));
ttc = bsxfun(@rdivide, ttc, max(ttc));
x = linspace(1, size(tc.JON,2), size(tc.JON,2));
plot(x, sort(nanmean(ttc)), 'Color', cols(1,:), 'LineWidth', 2);
RmeanRmaxJON = nanmean(ttc);

ttc = tc.AMMC;
ttc = bsxfun(@minus, ttc, min(ttc));
ttc = bsxfun(@rdivide, ttc, max(ttc));
x = linspace(1, size(tc.JON,2), size(tc.AMMC,2));
plot(x, sort(nanmean(ttc)), '-', 'Color', cols(2,:), 'LineWidth', 2);
RmeanRmaxAMMC = nanmean(ttc);

ttc = tc.WED;
ttc = bsxfun(@minus, ttc, min(ttc));
ttc = bsxfun(@rdivide, ttc, max(ttc));
x = linspace(1, size(tc.JON,2), size(tc.WED,2));
plot(x, sort(nanmean(ttc)), '-', 'Color', cols(3,:), 'LineWidth', 2);
RmeanRmaxWED = nanmean(ttc);

ttc = tc.TONOTipsi;
ttc = bsxfun(@minus, ttc, min(ttc));
ttc = bsxfun(@rdivide, ttc, max(ttc));
x = linspace(1, size(tc.JON,2), size(tc.TONOTipsi,2));
plot(x, sort(nanmean(ttc)), '-', 'Color', cols(4,:), 'LineWidth', 2);
RmeanRmaxTONOTipsi = nanmean(ttc);

ttc = tc.TONOTcontra;
ttc = bsxfun(@minus, ttc, min(ttc));
ttc = bsxfun(@rdivide, ttc, max(ttc));
x = linspace(1, size(tc.JON,2), size(tc.TONOTcontra,2));
plot(x, sort(nanmean(ttc)), '-', 'Color', cols(5,:), 'LineWidth', 2);
RmeanRmaxTONOTcontra = nanmean(ttc);

legend( 'JON', 'AMMC', 'WEDipsi', 'TONipsi', 'TONcontra')
legend boxoff


% export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/normalizedRmeans_JON_AMMC_WEDi_TONi_TONc_.eps')



allSpreads = cat(2, RmeanRmaxJON, RmeanRmaxAMMC, RmeanRmaxWED, RmeanRmaxTONOTipsi, RmeanRmaxTONOTcontra);
origin = cat(2, ones(size(RmeanRmaxJON)), 2*ones(size(RmeanRmaxAMMC)), 3*ones(size(RmeanRmaxWED)), 4*ones(size(RmeanRmaxTONOTipsi)), 5*ones(size(RmeanRmaxTONOTcontra)) );

figure; hold on
boxplot(allSpreads, origin);

ax = gca;
ax.XTickLabel = {'JON', 'AMMC', 'WEDi', 'TONi', 'TONc'};
ax.YLabel.String = {'R_m_e_a_n/R_m_a_x, interpolated tcs'};
ax.YLim = [0,1];

%% lifetime sparseness
nSt = length(x2);

ttc = tc.JON;
nums = nanmean(ttc).^2;
denoms = nanmean(ttc.^2);
lftsJON = 1 - nums./denoms;

ttc = tc.AMMC;
nums = nanmean(ttc).^2;
denoms = nanmean(ttc.^2);
lftsAMMC = 1 - nums./denoms;

ttc = tc.WED;
nums = nanmean(ttc).^2;
denoms = nanmean(ttc.^2);
lftsWED = 1 - nums./denoms;

ttc = tc.TONOTipsi;
nums = nanmean(ttc).^2;
denoms = nanmean(ttc.^2);
lftsTONOTipsi = 1 - nums./denoms;

ttc = tc.TONOTcontra;
nums = nanmean(ttc).^2;
denoms = nanmean(ttc.^2);
lftsTONOTcontra = 1 - nums./denoms;




allSpreads = cat(2, lftsJON, lftsAMMC, lftsWED, lftsTONOTipsi, lftsTONOTcontra);
% origin = cat(2, ones(size(RmeanRmaxJON)), 2*ones(size(RmeanRmaxAMMC)), 3*ones(size(RmeanRmaxWED)) );

figure; hold on
boxplot(allSpreads, origin);

ax = gca;
ax.XTickLabel = {'JON', 'AMMC', 'WEDi', 'TONi', 'TONc'};
ax.YLabel.String = {'lifetime sparseness, interpolated tcs'};

%% compare mu with halfcdf (--> overall looks fine)
% % 1. single pixel delta: boxplot 
% allHalfCdfs = cat(2, cdf.half_i_JON, cdf.half_i_AMMC, cdf.half_i_WED);
% allMus = cat(2, muJON, muAMMC, muWED); %this is in points, and so is allHalfCdfs
% allMus(abs(allMus)>2e3) = nan;
% 
% deltas = abs(allHalfCdfs - allMus);
% figure; hold on
% boxplot(deltas, origin)
% ylabel({'absolute delta (mu vs halfCDF)'; 'log spaced points'});
% ax = gca;
% ax.XTickLabel = {'JON', 'AMMC', 'WED'};
% export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/delta_MuBest_Distrib_JON_WEDvsAMMC_mean.eps')
% 
% 
% % 2. plot mu distribution
% cols = brewermap(3, 'Dark2');
% figure; hold on
% setN = 1; %JONS;
% sortedmu = sort(muJON(isfinite(muJON)));
% sortedmu(abs(sortedmu)>2e3) = nan;
% l(setN) = line(sortedmu, 1:length(sortedmu),'Color', cols(setN, :)); 
% 
% ax(setN) = gca;
% ax(setN).YLim(2) = length(cdf.half_i_JON);
% % ax(setN).XLim = [x2(1), x2Zero(end)];
% ax(setN).XColor = cols(setN, :);
% ax(setN).YColor = cols(setN, :);
% ax(setN).YTick = [0, 1e4];
% ax(setN).YMinorTick	 = 'on';
% ax(setN).YAxis.MinorTickValues = round(length(cdf.half_i_JON)/2);
% ax(setN).YMinorGrid = 'on';
% ax(setN).XGrid = 'on';
% ax(setN).YLabel.String = 'pixels sorted by best frequency';
% % ax(setN).XLabel.String = 'best frequency (half maximum cdf)';
% 
% pos = ax(setN).Position; % position of first axes
% 
% setN = 2;
% ax(setN) = axes('Position', pos,...
%     'XAxisLocation','top',...
%     'YAxisLocation','left',...
%     'Color','none');
% sortedmu = sort(muAMMC(isfinite(muAMMC)));
% sortedmu(abs(sortedmu)>2e3) = nan;
% l(setN) = line(sortedmu, 1:length(sortedmu),'Color', cols(setN, :)); 
% 
% % ax(setN).XAxis.Visible = 'off';
% ax(setN).YLim(2) = length(sortedmu);
% % ax(setN).XLim = [x2(1), x2Zero(end)];
% ax(setN).XColor = cols(setN, :);
% ax(setN).YColor = cols(setN, :);
% ax(setN).YTick = [0, 1e4];
% 
% 
% setN = 3;
% ax(setN) = axes('Position', pos,...
%     'XAxisLocation','top',...
%     'YAxisLocation','left',...
%     'Color','none');
% 
% sortedmu = sort(muWED(isfinite(muWED)));
% % remove outliers
% sortedmu(abs(sortedmu)>2e3) = nan;
% l(setN) = line(sortedmu, 1:length(sortedmu),'Color', cols(setN, :)); 
% 
% % ax(setN).XAxis.Visible = 'off';
% ax(setN).YLim(2) = length(sortedmu);
% % ax(setN).XLim = [x2(1), x2Zero(end)];
% ax(setN).XColor = cols(setN, :);
% ax(setN).YColor = cols(setN, :);
% ax(setN).YTick = [0, 1e4];
% 
% linkaxes(ax, 'x')
% ax(1).XLim = [1,2000];
% for i = 1:3
%     ax(i).XScale = 'linear';
% end
% legend(l, 'JON', 'AMMC', 'WED')
% legend boxoff
% 
% ax(1).XTick = 0:200:2000;
% ax(1).XTick = ax(1).XTick(2:end);
% ax(1).XTickLabel = round(x2(ax(1).XTick));
% 
% % export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/muDistrib_JON_WEDvsAMMC.eps')

%% scatter sigma vs mean
figure;
subplot(3,1,1), hold on
musc = muJON;
sigsc = sigmaJON;
title('JON')
ax(1) = gca;
removeoutliers = ~isfinite(musc) | ~isfinite(sigsc) | musc>4e3 | musc<-2e3 | sigsc>5e3;
scatter(musc(~removeoutliers), sigsc(~removeoutliers));
xlabel('mu')
ylabel('sigma')
axis image

subplot(3,1,2), hold on
musc = muAMMC;
sigsc = sigmaAMMC;
title('AMMC')
ax(2) = gca;
removeoutliers = ~isfinite(musc) | ~isfinite(sigsc) | musc>4e3 | musc<-2e3 | sigsc>4e3;
sum(removeoutliers)
scatter(musc(~removeoutliers), sigsc(~removeoutliers));
xlabel('mu')
ylabel('sigma')
axis image

subplot(3,1,3), hold on
musc = muWED;
sigsc = sigmaWED;
title('WED')
ax(3) = gca;
removeoutliers = ~isfinite(musc) | ~isfinite(sigsc) | musc>4e3 | musc<-2e3 | sigsc>4e3;
sum(removeoutliers)
scatter(musc(~removeoutliers), sigsc(~removeoutliers));
xlabel('mu')
ylabel('sigma')
axis image

ax(1).XLim = [-0.5 2.5].*1e3;
linkaxes(ax, 'xy')

export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/scatter_muVsSigma_JON_WEDvsAMMC.eps')
clear musc sigsc ax ttc


%% spread in octaves
% for each pixel center to its mean and consider the following octaves:
xOctaves = -3:3;
xOctavesFreqMultiplier = 2.^xOctaves;

musc = muJON;
sigsc = sigmaJON;
removeoutliers = ~isfinite(musc) | ~isfinite(sigsc) | musc>2e3 | musc<-0.5 | sigsc>5e3;
musc(removeoutliers) = [];
tcc = tc.JON(:,~removeoutliers);

tic
for p = 1 : sum(~removeoutliers)
    mu1 = x2(round(musc(p)));
    octDistJON(p,:) = interp1(x2, tcc(:,p), mu1.*xOctavesFreqMultiplier);
end
toc %12 sec



musc = muAMMC;
sigsc = sigmaAMMC;
removeoutliers = ~isfinite(musc) | ~isfinite(sigsc) | musc>2e3 | musc<-0.5 | sigsc>5e3;
musc(removeoutliers) = [];
tcc = tc.AMMC(:,~removeoutliers);

tic
for p = 1 : sum(~removeoutliers)
    mu1 = x2(round(musc(p)));
    octDistAMMC(p,:) = interp1(x2, tcc(:,p), mu1.*xOctavesFreqMultiplier);
end
toc


musc = muWED;
sigsc = sigmaWED;
removeoutliers = ~isfinite(musc) | ~isfinite(sigsc) | musc>2e3 | musc<-0.5 | sigsc>5e3;
musc(removeoutliers) = [];
tcc = tc.WED(:,~removeoutliers);

tic
for p = 1 : sum(~removeoutliers)
    mu1 = x2(round(musc(p)));
    octDistWED(p,:) = interp1(x2, tcc(:,p), mu1.*xOctavesFreqMultiplier);
end
toc

%%
% how do I normalize?
octDistJONn = bsxfun(@rdivide, octDistJON, octDistJON(:,4));
octDistAMMCn = bsxfun(@rdivide, octDistAMMC, octDistAMMC(:,4));
octDistWEDn = bsxfun(@rdivide, octDistWED, octDistWED(:,4));

errJON = bsxfun(@rdivide, nanstd(octDistJONn), sqrt(sum(isnan(octDistJONn))) );
errAMMC = bsxfun(@rdivide, nanstd(octDistAMMCn), sqrt(sum(isnan(octDistAMMCn))) );
errWED = bsxfun(@rdivide, nanstd(octDistWEDn), sqrt(sum(isnan(octDistWEDn))) );

figure; hold on
errorbar(xOctaves, nanmean(octDistJONn), errJON, 'Color', cols(1,:))
errorbar(xOctaves, nanmean(octDistAMMCn), errAMMC, 'Color', cols(2,:))
errorbar(xOctaves, nanmean(octDistWEDn), errWED, 'Color', cols(3,:))
legend('JON', 'AMMC', 'WED')
xlabel('octaves')

export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/normResponses2octaves_JON_WEDvsAMMC_mean.eps')

