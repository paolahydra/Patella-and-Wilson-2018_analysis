% here pixel wise FT and separation between AMMC and WED has already been
% calculated.
% all tcs are log interpolated in 2000 points here

%% load data and relevant indices
load('/Users/galileo/Dropbox (HMS)/Data/chunks_stimfinal.mat')

% load JON
JON = load('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix/pixelwisePeakFT_mean.mat');%NEW - mean of 13 points, rather than peak
% load AMMC + WED
CNS = load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/pixelwisePeakFT_mean.mat');%NEW - mean of 13 points, rather than peak
CNS.datalist = load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/datalist_WEDAMMC_piezo.mat');

CNS.clusterFolder = '/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix';
JON.clusterFolder = '/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix';

CNS.clusters = load(fullfile(CNS.clusterFolder, 'klust_maxClust19.mat'));
JON.clusters = load(fullfile(JON.clusterFolder, 'klust_maxClust10.mat'));
JON.excludeKs = [8,9,10]; % brown , pull, push
CNS.excludeKs = [1,2]; % pull, push
JON.includedK_i = ~ismember(JON.clusters.Tparent, JON.excludeKs);
CNS.includedK_i = ~ismember(CNS.clusters.Tparent, CNS.excludeKs);

%
[xvar, ~, ~] = unique([chunks.pips.carrierFreqLevels; chunks.tones.carrierFreqLevels(1)]);
onset = 4;
offset = 600; %cropping before, makes things worse
cutoff = 0.5;
x2 = logspace(log10(onset), log10(offset), 2000);
x2Zero = cat(2, 0, x2); 

%% JONs freq-responsive pixels cdfs
tc.JON = interp1(xvar, JON.dff_peak(JON.includedK_i,:)', x2 );
cdf.JON = bsxfun(@rdivide, cumsum(tc.JON), max(cumsum(tc.JON)) );

%% AMMC freq-responsive pixels cdfs
tc.AMMC = interp1(xvar, CNS.dff_peak(CNS.includedK_i & CNS.pxAMMC,:)', x2 );
cdf.AMMC = bsxfun(@rdivide, cumsum(tc.AMMC), max(cumsum(tc.AMMC)) );

%% WED freq-responsive pixels cdfs
tc.WED = interp1(xvar, CNS.dff_peak(CNS.includedK_i & CNS.pxWED,:)', x2 );
cdf.WED = bsxfun(@rdivide, cumsum(tc.WED), max(cumsum(tc.WED)) );

%% find best frequency (half maximum point in the cdf, pixel wise)
[cdf.half_i_WED, ~] = sort(sum(cdf.WED < 0.5));
[cdf.half_i_AMMC, ~] = sort(sum(cdf.AMMC < 0.5));
[cdf.half_i_JON, ~] = sort(sum(cdf.JON < 0.5));

%% find best frequency (half maximum point in the cdf, pixel wise) - YET NOTE that I still need to select included pixels here
m = movmean(tc.JON, 20);
[~,mi] = max(m);
bestFrJON = x2(mi);

m = movmean(tc.AMMC, 20);
[~,mi] = max(m);
bestFrAMMC = x2(mi);

m = movmean(tc.WED, 20);
[~,mi] = max(m);
bestFrWED = x2(mi);

%% adapt
NBplusplus = 20;
edges = logspace(log10(onset), log10(offset),NBplusplus);
figure; hold on
h1 = histogram(bestFrJON,edges);
h2 = histogram(bestFrAMMC,edges);
h3 = histogram(bestFrWED,edges);


ax = gca;
ax.XAxis.Scale = 'log';
ax.TickDir = 'out';
ax.Box = 'off';
ax.XLim = [onset, offset];

h1.Normalization = 'probability';
h2.Normalization = 'probability';
h3.Normalization = 'probability';

% h1.FaceAlpha = 1;
% h2.FaceAlpha = 1;
% h3.FaceAlpha = 1; %save vectorial, change later
% 

% export_fig(fullfile('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning',...
%                     'JON_AMMC_WED_pdfHistogram_LOG_alpha40.eps'))

centers = zeros(1, NBplusplus-1);
for i = 1:NBplusplus-1
    t = logspace(log10(edges(i)), log10(edges(i+1)),3);
    centers(i) = t(2);
end

figure;
hold on
plot(centers, h1.Values);
plot(centers, h2.Values);
plot(centers, h3.Values);
ax = gca;
ax.XAxis.Scale = 'log';
ax.TickDir = 'out';

ax.XAxis.Limits = [onset, offset];
ax.XAxis.MinorTickValues = [4:10   20:10:100   200:100:600];
ax.XAxis.TickValues = [8,16,25,50,100,225, 350];


ax.XAxis.TickValues = [8,16,25,50,100,225, 350, 600];

% 
% 
% export_fig(fullfile('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning',...
%                     'JON_AMMC_WED_pdfCurves_LOG.eps'))



%% plot #PIXELS AS A FUNCTION OF FREQUENCY (AND NOT THE OTHER WAY AROUND)
cols = brewermap(3, 'Dark2');
figure; hold on
setN = 1; %JONS;
l(setN) = line(x2Zero(cdf.half_i_JON+1), 1:length(cdf.half_i_JON),'Color', cols(setN, :)); 

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

l(setN) = line(x2Zero(cdf.half_i_AMMC+1), 1:length(cdf.half_i_AMMC),'Color', cols(setN, :)); 

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

l(setN) = line(x2Zero(cdf.half_i_WED+1), 1:length(cdf.half_i_WED),'Color', cols(setN, :)); 

ax(setN).XAxis.Visible = 'off';
ax(setN).YLim(2) = length(cdf.half_i_WED);
ax(setN).XLim = [x2(1), x2Zero(end)];
ax(setN).XColor = cols(setN, :);
ax(setN).YColor = cols(setN, :);
ax(setN).YTick = [0, 1e4];

for i = 1:3
    ax(i).XScale = 'log';
end
legend(l, 'JON', 'AMMC', 'WED')
legend boxoff

export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/preferredDistrib_JON_WEDvsAMMC_mean.eps')


%% find spread of tuning curve (sd of fitted gaussian over log spaced freq tuning curves)
% tuning curves are sampled logarithmically, but that is vain if I plot it
% against x2, which is also sampled logarithmically.
% x2
tic
warning off
ttc = tc.WED;
sigmaWED = nan(1, length(ttc));
muWED = nan(1, length(ttc));
ttc = ttc';
for i = 1 : length(ttc)
    [sigmaWED(i), muWED(i)] = gaussfit( 1:2000, ttc(i,:));
end
disp('Done WED')

ttc = tc.AMMC;
sigmaAMMC = nan(1, length(ttc));
muAMMC = nan(1, length(ttc));
ttc = ttc';
for i = 1 : length(ttc)
    [sigmaAMMC(i), muAMMC(i)] = gaussfit( 1:2000, ttc(i,:));
end
disp('Done AMMC')

ttc = tc.JON;
sigmaJON = nan(1, length(ttc));
muJON = nan(1, length(ttc));
ttc = ttc';
for i = 1 : length(ttc)
    [sigmaJON(i), muJON(i)] = gaussfit( 1:2000, ttc(i,:));
end
disp('Done JON')
toc


%%
allSigmas = cat(2, sigmaJON, sigmaAMMC, sigmaWED);
allSigmas(allSigmas>2e3) = nan;
origin = cat(2, ones(size(sigmaJON)), 2*ones(size(sigmaAMMC)), 3*ones(size(sigmaWED)) );
%
cols = brewermap(3, 'Dark2');
figure; hold on
boxplot(allSigmas, origin)
ax = gca;
ax.XTickLabel = {'JON', 'AMMC', 'WED'};
ax.YLabel.String = {'sigma of fitted gaussian'; 'log spaced points'};

export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/sigmaDistrLogspaced_JON_WEDvsAMMC_mean.eps')

%% RmeanRmax (tcs are already interpolated - 2000)
ttc = tc.WED;
RmeanRmaxWED = nanmean(ttc)./ max(ttc);

ttc = tc.AMMC;
RmeanRmaxAMMC = nanmean(ttc)./ max(ttc);

ttc = tc.JON;
RmeanRmaxJON = nanmean(ttc)./ max(ttc);


allSpreads = cat(2, RmeanRmaxJON, RmeanRmaxAMMC, RmeanRmaxWED);
origin = cat(2, ones(size(RmeanRmaxJON)), 2*ones(size(RmeanRmaxAMMC)), 3*ones(size(RmeanRmaxWED)) );

figure; hold on
boxplot(allSpreads, origin);

ax = gca;
ax.XTickLabel = {'JON', 'AMMC', 'WED'};
ax.YLabel.String = {'R_m_e_a_n/R_m_a_x, interpolated tcs'};

% 
figure; hold on
plot(sort(nanmean(tc.JON)), '--', 'Color', cols(1,:))
plot(sort(nanmean(tc.AMMC)), '--', 'Color', cols(2,:))
plot(sort(nanmean(tc.WED)), '--', 'Color', cols(3,:))

plot(sort(max(tc.JON)), '-', 'Color', cols(1,:))
plot(sort(max(tc.AMMC)), '-', 'Color', cols(2,:))
plot(sort(max(tc.WED)), '-', 'Color', cols(3,:))

set(gca, 'XScale', 'log')
set(gca, 'YLim', [-25, 500])
ylabel('delta F / F')

%% RmeanRmax (tcs are already interpolated - 2000) - normalize between 0 and 1 first
figure; hold on
title('Rmean distributions after normalization')

ttc = tc.JON;
ttc = bsxfun(@minus, ttc, min(ttc));
ttc = bsxfun(@rdivide, ttc, max(ttc));
x = linspace(1, size(tc.JON,2), size(tc.JON,2));
plot(x, sort(nanmean(ttc)), '--', 'Color', cols(1,:))
RmeanRmaxJON = nanmean(ttc);

ttc = tc.AMMC;
ttc = bsxfun(@minus, ttc, min(ttc));
ttc = bsxfun(@rdivide, ttc, max(ttc));
x = linspace(1, size(tc.AMMC,2), size(tc.JON,2));
plot(x, sort(nanmean(ttc)), '--', 'Color', cols(2,:))
RmeanRmaxAMMC = nanmean(ttc);

ttc = tc.WED;
ttc = bsxfun(@minus, ttc, min(ttc));
ttc = bsxfun(@rdivide, ttc, max(ttc));
x = linspace(1, size(tc.WED,2), size(tc.JON,2));
plot(x, sort(nanmean(ttc)), '--', 'Color', cols(3,:))
RmeanRmaxWED = nanmean(ttc);


export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/normalizedRmeans.eps')



allSpreads = cat(2, RmeanRmaxJON, RmeanRmaxAMMC, RmeanRmaxWED);
origin = cat(2, ones(size(RmeanRmaxJON)), 2*ones(size(RmeanRmaxAMMC)), 3*ones(size(RmeanRmaxWED)) );

figure; hold on
boxplot(allSpreads, origin);

ax = gca;
ax.XTickLabel = {'JON', 'AMMC', 'WED'};
ax.YLabel.String = {'R_m_e_a_n/R_m_a_x, interpolated tcs'};
ax.YLim = [0,1];

% 

plot(sort(nanmean(tc.JON)), '--', 'Color', cols(1,:))
plot(sort(nanmean(tc.AMMC)), '--', 'Color', cols(2,:))
plot(sort(nanmean(tc.WED)), '--', 'Color', cols(3,:))

plot(sort(max(tc.JON)), '-', 'Color', cols(1,:))
plot(sort(max(tc.AMMC)), '-', 'Color', cols(2,:))
plot(sort(max(tc.WED)), '-', 'Color', cols(3,:))

set(gca, 'XScale', 'log')
set(gca, 'YLim', [-25, 500])
ylabel('delta F / F')

%% lifetime sparseness
nSt = length(x2);
ttc = tc.WED;
nums = nanmean(ttc).^2;
denoms = nanmean(ttc.^2);
lftsWED = 1 - nums./denoms;

ttc = tc.AMMC;
nums = nanmean(ttc).^2;
denoms = nanmean(ttc.^2);
lftsAMMC = 1 - nums./denoms;

ttc = tc.JON;
nums = nanmean(ttc).^2;
denoms = nanmean(ttc.^2);
lftsJON = 1 - nums./denoms;


allSpreads = cat(2, lftsJON, lftsAMMC, lftsWED);
% origin = cat(2, ones(size(RmeanRmaxJON)), 2*ones(size(RmeanRmaxAMMC)), 3*ones(size(RmeanRmaxWED)) );

figure; hold on
boxplot(allSpreads, origin);

ax = gca;
ax.XTickLabel = {'JON', 'AMMC', 'WED'};
ax.YLabel.String = {'lifetime sparseness, interpolated tcs'};

%% compare mu with halfcdf (--> overall looks fine)
% 1. single pixel delta: boxplot 
allHalfCdfs = cat(2, cdf.half_i_JON, cdf.half_i_AMMC, cdf.half_i_WED);
allMus = cat(2, muJON, muAMMC, muWED); %this is in points, and so is allHalfCdfs
allMus(abs(allMus)>2e3) = nan;

deltas = abs(allHalfCdfs - allMus);
figure; hold on
boxplot(deltas, origin)
ylabel({'absolute delta (mu vs halfCDF)'; 'log spaced points'});
ax = gca;
ax.XTickLabel = {'JON', 'AMMC', 'WED'};
export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/delta_MuBest_Distrib_JON_WEDvsAMMC_mean.eps')


% 2. plot mu distribution
cols = brewermap(3, 'Dark2');
figure; hold on
setN = 1; %JONS;
sortedmu = sort(muJON(isfinite(muJON)));
sortedmu(abs(sortedmu)>2e3) = nan;
l(setN) = line(sortedmu, 1:length(sortedmu),'Color', cols(setN, :)); 

ax(setN) = gca;
ax(setN).YLim(2) = length(cdf.half_i_JON);
% ax(setN).XLim = [x2(1), x2Zero(end)];
ax(setN).XColor = cols(setN, :);
ax(setN).YColor = cols(setN, :);
ax(setN).YTick = [0, 1e4];
ax(setN).YMinorTick	 = 'on';
ax(setN).YAxis.MinorTickValues = round(length(cdf.half_i_JON)/2);
ax(setN).YMinorGrid = 'on';
ax(setN).XGrid = 'on';
ax(setN).YLabel.String = 'pixels sorted by best frequency';
% ax(setN).XLabel.String = 'best frequency (half maximum cdf)';

pos = ax(setN).Position; % position of first axes

setN = 2;
ax(setN) = axes('Position', pos,...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none');
sortedmu = sort(muAMMC(isfinite(muAMMC)));
sortedmu(abs(sortedmu)>2e3) = nan;
l(setN) = line(sortedmu, 1:length(sortedmu),'Color', cols(setN, :)); 

% ax(setN).XAxis.Visible = 'off';
ax(setN).YLim(2) = length(sortedmu);
% ax(setN).XLim = [x2(1), x2Zero(end)];
ax(setN).XColor = cols(setN, :);
ax(setN).YColor = cols(setN, :);
ax(setN).YTick = [0, 1e4];


setN = 3;
ax(setN) = axes('Position', pos,...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none');

sortedmu = sort(muWED(isfinite(muWED)));
% remove outliers
sortedmu(abs(sortedmu)>2e3) = nan;
l(setN) = line(sortedmu, 1:length(sortedmu),'Color', cols(setN, :)); 

% ax(setN).XAxis.Visible = 'off';
ax(setN).YLim(2) = length(sortedmu);
% ax(setN).XLim = [x2(1), x2Zero(end)];
ax(setN).XColor = cols(setN, :);
ax(setN).YColor = cols(setN, :);
ax(setN).YTick = [0, 1e4];

linkaxes(ax, 'x')
ax(1).XLim = [1,2000];
for i = 1:3
    ax(i).XScale = 'linear';
end
legend(l, 'JON', 'AMMC', 'WED')
legend boxoff

ax(1).XTick = 0:200:2000;
ax(1).XTick = ax(1).XTick(2:end);
ax(1).XTickLabel = round(x2(ax(1).XTick));

export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/muDistrib_JON_WEDvsAMMC.eps')

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

