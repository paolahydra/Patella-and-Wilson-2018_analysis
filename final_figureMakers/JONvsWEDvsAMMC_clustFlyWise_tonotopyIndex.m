% tuning data previously assembled in JON_pixelwise_distributions and WEDvsAMMC_pixelwise_distributions
% now reassmbled in 
% edit WEDvsAMMC_pixelwise_distributions_selectZPlanesAndRegions.m


% using geom2d toolbox

%% you want to recalculate these based on the selected runs/regions
% load JON
JON = load('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix/ftcJON_clusterTuningFly.mat');
% load AMMC + WED
CNS = load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/ftc_AMMCi_TONOTi_TONOTc_clusterTuningFly_mean.mat');
%CNS.tcAMMC_clusterTuningFly % 19x19x11
% save('ftc_AMMCi_TONOTi_TONOTc_clusterTuningFly.mat', 'tcAMMCipsi_clusterTuningFly', 'tcTONOTipsi_clusterTuningFly', 'tcTONOTcontra_clusterTuningFly', ...
%     'tcAMMCipsi_sumPixels_clusterFly', 'tcTONOTipsi_sumPixels_clusterFly', 'tcTONOTcontra_sumPixels_clusterFly') % all 19 clusters

%% load everything useful

load('/Users/galileo/Dropbox (HMS)/Data/chunks_stimfinal.mat')
load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/iCNS.mat', 'iCNS')
load('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix/iJON.mat', 'iJON')

CNS.clusterFolder = '/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix';
JON.clusterFolder = '/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix';


%the following fileds have now been loaded in iCNS:

% CNS.datalist = load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/datalist_WEDAMMC_piezo.mat');
% CNS.clusters = load(fullfile(CNS.clusterFolder, 'klust_maxClust19.mat'));
% CNS.excludeKs = [1,2]; % pull, push
% CNS.includedK_i = ~ismember(CNS.clusters.Tparent, CNS.excludeKs);
iCNS.includedK_i = ~ismember(iCNS.Tparent, iCNS.excludeKs);

% JON.clusters = load(fullfile(JON.clusterFolder, 'klust_maxClust10.mat'));
% JON.excludeKs = [8,9,10]; % brown , pull, push
% iJON.includedK_i = ~ismember(JON.clusters.Tparent, JON.excludeKs);
iJON.includedK_i = ~ismember(iJON.Tparent, iJON.excludeKs);


% single-pixel tuning curves
load('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix/pixelwisePeakFT_mean.mat', 'dff_peak', 'bestFr_interpolated');
JON.dff_peak = dff_peak;
JON.bestFreq_int = bestFr_interpolated;
load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/pixelwisePeakFT_mean.mat', 'dff_peak', 'bestFr_interpolated'); %171204
CNS.dff_peak = dff_peak;
CNS.bestFreq_int = bestFr_interpolated;
clear dff_peak bestFr_interpolated

%
[xvar, ~, ~] = unique([chunks.pips.carrierFreqLevels; chunks.tones.carrierFreqLevels(1)]);
% we are not interpolating here, and I still think it is not such a good idea.
onset = 4;
offset = 600; %cropping before, makes things worse
cutoff = 0.5;
x2 = logspace(log10(onset), log10(offset), 2000);
x3 = logspace(log10(onset), log10(offset), 120);
x2Zero = cat(2, 0, x2);   


sizeMap = [60, 86];
% keepClustersCNS = 3:19;
xcoord = repmat(1:sizeMap(2), sizeMap(1), 1);
ycoord = repmat((1:sizeMap(1))', 1, sizeMap(2));
xcoord = xcoord(:);
ycoord = ycoord(:);


%% strip chart, with pre-saved correlation values   ------  mapWISE
corrTonotopies(1) = load('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix/tonotopyIndex_JONS.mat', 'correlations', 'rotations');
corrTonotopies(2) = load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/tonotopyIndex_AMMC.mat', 'correlations', 'rotations');
corrTonotopies(3) = load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/tonotopyIndex_iWED.mat', 'correlations', 'rotations');
% corrTonotopies(4) = load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/tonotopyIndex_cWED.mat', 'correlations', 'rotations');





allSpreads = cat(2, cat(2, corrTonotopies(1).correlations.regression2Rot), ...
                    cat(2, corrTonotopies(2).correlations.regression2Rot), ...
                    cat(2, corrTonotopies(3).correlations.regression2Rot) ); %, ...
%                     cat(2, corrTonotopies(3).correlations.regression2Rot), cat(2, corrTonotopies(4).correlations.regression2Rot), ...
%                     cat(2, corrTonotopies(4).correlations.regression2Rot) );
origin = cat(2, 1*ones(1, length(corrTonotopies(1).correlations)), ...
                2*ones(1, length(corrTonotopies(2).correlations)), ...
                3*ones(1, length(corrTonotopies(3).correlations)) ); %, ...
%                 4*ones(1, length(corrTonotopies(3).correlations)+length(corrTonotopies(4).correlations)), ...
%                 5*ones(1, length(corrTonotopies(4).correlations)) );

figure; hold on
% notBoxPlot(allSpreads, origin);
[H] = notBoxPlot(allSpreads, origin, 'style', 'line', 'markMedian', true);

ax = gca;
ax.YLim = [0,1]; %checked that no cropping
ax.XTickLabel =  {'JON', 'AMMCi', 'WEDi'}; %, 'WEDb', 'WEDc'}; 
ax.YLabel.String = {'freq-posit regression'};
ax.FontSize = 7;

export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/tonotopicIndex_selectedZPlanes.eps')

%\% establish significance (pairwise ttest)
figure; normplot(cat(2, corrTonotopies(3).correlations.regression2Rot)); title('iwed')
lillie(3) = lillietest(cat(2, corrTonotopies(3).correlations.regression2Rot));
jbt(3) = jbtest(cat(2, corrTonotopies(3).correlations.regression2Rot));

figure; normplot(cat(2, corrTonotopies(2).correlations.regression2Rot)); title('ammc')
lillie(2) = lillietest(cat(2, corrTonotopies(2).correlations.regression2Rot));
jbt(2) = jbtest(cat(2, corrTonotopies(2).correlations.regression2Rot));

figure; normplot(cat(2, corrTonotopies(1).correlations.regression2Rot)); title('jon')
lillie(1) = lillietest(cat(2, corrTonotopies(1).correlations.regression2Rot));
jbt(1) = jbtest(cat(2, corrTonotopies(1).correlations.regression2Rot));

%ok, normal distribution it is.

% check variances too
[var12.h, var12.p] = vartest2(cat(2, corrTonotopies(1).correlations.regression2Rot) , cat(2, corrTonotopies(2).correlations.regression2Rot))
[var13.h, var13.p] = vartest2(cat(2, corrTonotopies(1).correlations.regression2Rot) , cat(2, corrTonotopies(3).correlations.regression2Rot))
[var23.h, var23.p] = vartest2(cat(2, corrTonotopies(2).correlations.regression2Rot) , cat(2, corrTonotopies(3).correlations.regression2Rot))


% ttest
[tt12.h , tt12.p ] = ttest2(cat(2, corrTonotopies(1).correlations.regression2Rot) , cat(2, corrTonotopies(2).correlations.regression2Rot));
[tt13.h , tt13.p ] = ttest2(cat(2, corrTonotopies(1).correlations.regression2Rot) , cat(2, corrTonotopies(3).correlations.regression2Rot) ,'Vartype','unequal');
[tt23.h , tt23.p ] = ttest2(cat(2, corrTonotopies(2).correlations.regression2Rot) , cat(2, corrTonotopies(3).correlations.regression2Rot) ,'Vartype','unequal'); 

%same overall results if assuming same variances in all cases



%% strip chart, with pre-saved correlation values      -------    flyWISE
corrTonotopies(1) = load('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix/tonotopyIndex_JONS.mat', 'correlations', 'rotations');
corrTonotopies(2) = load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/tonotopyIndex_AMMC.mat', 'correlations', 'rotations');
corrTonotopies(3) = load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/tonotopyIndex_iWED.mat', 'correlations', 'rotations');

% fly grouping
% JONs
cti = 1;
fliesZetas = iJON.dataTable.fly(iJON.PAN.keepZetas);
flies = unique(fliesZetas, 'stable');

for f = 1:length(flies)
    flyNum = flies(f);
    allZetasThisFly = find(fliesZetas==flyNum);
    allFreqSorted = [];
    allFreqRotated = [];
    for z = 1:length(allZetasThisFly)
        zeta = allZetasThisFly(z);
        allFreqSorted = cat(1, allFreqSorted, sort(corrTonotopies(cti).rotations(zeta).frRotated) );
        allFreqRotated = cat(1, allFreqRotated, corrTonotopies(cti).rotations(zeta).frRotated);
    end
    corrTonotFlywise(cti).regression2Rot(f) = corr(allFreqSorted, allFreqRotated);
end

% AMMC
cti = 2;
fliesZetas = iCNS.dataTable.fly(iCNS.AMMC.keepZetas);
flies = unique(fliesZetas, 'stable');

for f = 1:length(flies)
    flyNum = flies(f);
%     allZetasThisFly = iCNS.AMMC.keepZetas(fliesZetas==flyNum); %index in dataTable
    allZetasThisFly = find(fliesZetas==flyNum); %index in keepZetas
    allFreqSorted = [];
    allFreqRotated = [];
    for z = 1:length(allZetasThisFly)
        zeta = allZetasThisFly(z);
        allFreqSorted = cat(1, allFreqSorted, sort(corrTonotopies(cti).rotations(zeta).frRotated) );
        allFreqRotated = cat(1, allFreqRotated, corrTonotopies(cti).rotations(zeta).frRotated);
    end
    corrTonotFlywise(cti).regression2Rot(f) = corr(allFreqSorted, allFreqRotated);
end

% WED
cti = 3;
fliesZetas = iCNS.dataTable.fly(iCNS.WEDi.keepZetas);
flies = unique(fliesZetas, 'stable');

for f = 1:length(flies)
    flyNum = flies(f);
%     allZetasThisFly = iCNS.AMMC.keepZetas(fliesZetas==flyNum); %index in dataTable
    allZetasThisFly = find(fliesZetas==flyNum); %index in keepZetas
    allFreqSorted = [];
    allFreqRotated = [];
    for z = 1:length(allZetasThisFly)
        zeta = allZetasThisFly(z);
        allFreqSorted = cat(1, allFreqSorted, sort(corrTonotopies(cti).rotations(zeta).frRotated) );
        allFreqRotated = cat(1, allFreqRotated, corrTonotopies(cti).rotations(zeta).frRotated);
    end
    corrTonotFlywise(cti).regression2Rot(f) = corr(allFreqSorted, allFreqRotated);
end


save('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/corrTonotFlywise_figure6A.mat', 'corrTonotFlywise')


allSpreads = cat(2, cat(2, corrTonotFlywise(1).regression2Rot), ...
                    cat(2, corrTonotFlywise(2).regression2Rot), ...
                    cat(2, corrTonotFlywise(3).regression2Rot) ); %, ...
%                     cat(2, corrTonotopies(3).correlations.regression2Rot), cat(2, corrTonotopies(4).correlations.regression2Rot), ...
%                     cat(2, corrTonotopies(4).correlations.regression2Rot) );
origin = cat(2, 1*ones(1, length(corrTonotFlywise(1).regression2Rot)), ...
                2*ones(1, length(corrTonotFlywise(2).regression2Rot)), ...
                3*ones(1, length(corrTonotFlywise(3).regression2Rot)) ); %, ...
%                 4*ones(1, length(corrTonotopies(3).correlations)+length(corrTonotopies(4).correlations)), ...
%                 5*ones(1, length(corrTonotopies(4).correlations)) );

%% establish significance (pairwise ttests)
figure; normplot(cat(2, corrTonotFlywise(3).regression2Rot)); title('iwed')
lillie(3) = lillietest(cat(2, corrTonotFlywise(3).regression2Rot));
jbt(3) = jbtest(cat(2, corrTonotFlywise(3).regression2Rot));

figure; normplot(cat(2, corrTonotFlywise(2).regression2Rot)); title('ammc')
lillie(2) = lillietest(cat(2, corrTonotFlywise(2).regression2Rot));
jbt(2) = jbtest(cat(2, corrTonotFlywise(2).regression2Rot));

figure; normplot(cat(2, corrTonotFlywise(1).regression2Rot)); title('jon')
lillie(1) = lillietest(cat(2, corrTonotFlywise(1).regression2Rot));
jbt(1) = jbtest(cat(2, corrTonotFlywise(1).regression2Rot));
%ok, normal distribution it is.

% check variances too
[var12.h, var12.p] = vartest2(cat(2, corrTonotFlywise(1).regression2Rot) , cat(2, corrTonotFlywise(2).regression2Rot))
[var13.h, var13.p] = vartest2(cat(2, corrTonotFlywise(1).regression2Rot) , cat(2, corrTonotFlywise(3).regression2Rot))
[var23.h, var23.p] = vartest2(cat(2, corrTonotFlywise(2).regression2Rot) , cat(2, corrTonotFlywise(3).regression2Rot))


% ttest
[tt12.h , tt12.p ] = ttest2(cat(2, corrTonotFlywise(1).regression2Rot) , cat(2, corrTonotFlywise(2).regression2Rot));
[tt13.h , tt13.p ] = ttest2(cat(2, corrTonotFlywise(1).regression2Rot) , cat(2, corrTonotFlywise(3).regression2Rot));
[tt23.h , tt23.p ] = ttest2(cat(2, corrTonotFlywise(2).regression2Rot) , cat(2, corrTonotFlywise(3).regression2Rot)); 
tt12
tt13
tt23
%same overall results if assuming same variances in all cases

%% establish significance anova1 and tukey post hoc test
clear stats
[p,t,stats] = anova1(allSpreads,origin);
[resultsT] = multcompare(stats,'CType','tukey-kramer');
[resultsB] = multcompare(stats,'CType','bonferroni');


%%


figure; hold on
% notBoxPlot(allSpreads, origin);
[H] = notBoxPlot(allSpreads, origin, 'style', 'line', 'markMedian', true);

ax = gca;
ax.YLim = [0,1]; %checked that no cropping
ax.XTickLabel =  {'JON', 'AMMCi', 'WEDi'}; %, 'WEDb', 'WEDc'}; 
ax.YLabel.String = {'freq-posit regression'};
ax.FontSize = 7;

export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/tonotopicIndex_selectedZPlanes_flyWise.eps')


%% make pdfs center frequency pixelwise
pdf_cf_AMMC = sort(CNS.bestFreq_int(iCNS.AMMC.AMMCipsi_all));
pdf_cf_WEDi = sort(CNS.bestFreq_int(iCNS.WEDi.TONOTipsi_all));
pdf_cf_JON = sort(JON.bestFreq_int(iJON.PAN.keepPixelsFromZetas));
% %%
% figure;
% h1 = histogram(pdf_cf_AMMC);
% hold on
% h2 = histogram(pdf_cf_WEDi);
% h3 = histogram(pdf_cf_JON);
% 
% ax = gca;
% ax.XAxis.Scale = 'log';
% ax.TickDir = 'out';
% ax.Box = 'off';
% ax.XLim = [5, 500]
% 
% h1.Normalization = 'probability';
% h2.Normalization = 'probability';
% h3.Normalization = 'probability';
% 
% BW = 10;
% h1.BinWidth = BW;
% h2.BinWidth = BW;
% h3.BinWidth = BW;

%% log hist (makes more sense to me)
NBplusplus = 40;
edges = logspace(log10(5), log10(500),NBplusplus);
figure; hold on
h1 = histogram(pdf_cf_JON,edges);
h2 = histogram(pdf_cf_AMMC,edges);
h3 = histogram(pdf_cf_WEDi,edges);


ax = gca;
ax.XAxis.Scale = 'log';
ax.TickDir = 'out';
ax.Box = 'off';
ax.XLim = [5, 500];

h1.Normalization = 'probability';
h2.Normalization = 'probability';
h3.Normalization = 'probability';

% h1.FaceAlpha = 1;
% h2.FaceAlpha = 1;
% h3.FaceAlpha = 1; %save vectorial, change later
% 
export_fig(fullfile('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning',...
                    'JON_AMMC_WED_pdfHistogram_LOG_alpha40.eps'))
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

ax.XAxis.Limits = [5, 500];
ax.XAxis.MinorTickValues = [5:10   20:10:100   200:100:500];
ax.XAxis.TickValues = [8,16,25,50,100,225, 350];


ax.XAxis.TickValues = [8,16,25,50,100,225, 350, 600];



export_fig(fullfile('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning',...
                    'JON_AMMC_WED_pdfCurves_LOG.eps'))


%% 2018 01 31
% make pdfs BEST(==peak) frequency pixelwise
[~, mi] = max((CNS.dff_peak)');
bestPeak_Freq_CNS = xvar(mi);

[~, mi] = max((JON.dff_peak)');
bestPeak_Freq_JON = xvar(mi);


pdf_cf_AMMC = sort(bestPeak_Freq_CNS(iCNS.AMMC.AMMCipsi_all));
pdf_cf_WEDi = sort(bestPeak_Freq_CNS(iCNS.WEDi.TONOTipsi_all));
pdf_cf_JON = sort(bestPeak_Freq_JON(iJON.PAN.keepPixelsFromZetas));


% % % peak will not change no matter what - unless you fit a model (e.g.
% % % gaussian

% % interpolate
% tcJON = interp1(xvar, JON.dff_peak', x2 );
% tcCNS = interp1(xvar, CNS.dff_peak', x2 );
% 
% [~, mi] = max(tcCNS);
% bestPeak_Freq_CNS = x2(mi);
% 
% [~, mi] = max(tcJON);
% bestPeak_Freq_JON = x2(mi);
% 
% 
% pdf_cf_AMMC = sort(bestPeak_Freq_CNS(iCNS.AMMC.AMMCipsi_all));
% pdf_cf_WEDi = sort(bestPeak_Freq_CNS(iCNS.WEDi.TONOTipsi_all));
% pdf_cf_JON = sort(bestPeak_Freq_JON(iJON.PAN.keepPixelsFromZetas));



%% log hist (makes more sense to me)
NBplusplus = 10;
edges = logspace(log10(onset), log10(offset),NBplusplus);
figure; hold on
h1 = histogram(pdf_cf_JON,edges);
h2 = histogram(pdf_cf_AMMC,edges);
h3 = histogram(pdf_cf_WEDi,edges);


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
ax.XAxis.TickValues = [4,8,16,25,50,100,225, 350];


ax.XAxis.TickValues = [4,8,16,25,50,100,225, 350, 600];


title('10 bins, pixel tuning curves previously interpolated')
export_fig(fullfile('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning',...
                    'JON_AMMC_WED_pdfCurves_LOG_peakFreq_tc_interpolated.eps'))

                
                
                
                
                
%% TEMP - UNA TANTUM: for all CNS and JON pixels, calculate best freq based on half cdf
% yvar2 = interp1(xvar, CNS.dff_peak', x2 );
% cdfs2 = bsxfun(@rdivide, cumsum(yvar2), max(cumsum(yvar2)) );
% i_cdf = cdfs2 < 0.5;
% s_cdf = sum(i_cdf);       % interpolated
% bestFr_interpolated = nan(size(s_cdf));
% for i = 1:length(s_cdf)
%     if s_cdf(i)>0
%         bestFr_interpolated(i) = x2(s_cdf(i));
%     end
% end
% CNS.bestFreq_int = bestFr_interpolated;
% save('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/pixelwisePeakFT_mean.mat', 'bestFr_interpolated', '-append');
% 
% yvar2 = interp1(xvar, JON.dff_peak', x2 );
% cdfs2 = bsxfun(@rdivide, cumsum(yvar2), max(cumsum(yvar2)) );
% i_cdf = cdfs2 < 0.5;
% s_cdf = sum(i_cdf);       % interpolated
% bestFr_interpolated = nan(size(s_cdf));
% for i = 1:length(s_cdf)
%     if s_cdf(i)>0
%         bestFr_interpolated(i) = x2(s_cdf(i));
%     end
% end
% JON.bestFreq_int = bestFr_interpolated;
% save('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix/pixelwisePeakFT_mean.mat', 'bestFr_interpolated', '-append');


%% contraTONOT
for allZetas = 1: length(iCNS.WEDc.keepZetas)
    
    keepZeta = iCNS.WEDc.keepZetas(allZetas);
    keepPixelsZeta = false(size(iCNS.pixel2keep));
    keepPixelsZeta(:,keepZeta) = true;
    keepPixelsZeta = keepPixelsZeta(:);
    keepPixelsZeta = keepPixelsZeta(iCNS.pixel2keep(:)); % indices of all pixel2keep this run % size: 49608,1 - sum: 1126
    
    % apply other masks... and exclude sust clusters' pixels, then convert back
    keepPixels = keepPixelsZeta & iCNS.includedK_i & iCNS.WEDc.TONOTcontra_all; % size: 49608,1 - sum: 764
    
    % iterative approach for excluding some more spurious pixels
    freqmap = CNS.bestFreq_int(keepPixels);                                     % size: 764,1
    excludePixels = find(isnan(freqmap));
    temp = cumsum(keepPixels);
    for i = 1:length(excludePixels)
        excludePixels_i = find(temp == excludePixels(i), 1);
        keepPixels(excludePixels_i) = 0;
    end
    clear temp excludePixels_i
    keepPixels_thisRunPx2Kp = keepPixels(keepPixelsZeta);                       % size: 1126,1  - sum: 764
    freqmap = CNS.bestFreq_int(keepPixels);
    
    
    thismap = iCNS.pixel2keep(:, keepZeta);                                     % size: 5060,1 - sum: 1126    % all px2KP
    thisMap_keepPixels = false(size(thismap));
    thisMap_keepPixels(thismap) = keepPixels_thisRunPx2Kp;                      % size: 5060,1 - sum: 764     % select px2KP
    
    % reshape given maps:
    freqMAPreshaped = nan(size(thismap));
    freqMAPreshaped(thisMap_keepPixels) = freqmap;
    freqMAPreshaped = reshape(freqMAPreshaped, sizeMap);
%     figure; imshow(freqMAPreshaped, []);
    
    keepPixelsMAPreshaped = thisMap_keepPixels;
    keepPixelsMAPreshaped = reshape(keepPixelsMAPreshaped, sizeMap);
%     figure; imshow(keepPixelsMAPreshaped, []); title('no pushpull clusters')

    
    %% manual rotation of points and projection on x axis (x coord)
    freqREF = (sort(freqmap))';
    %take the ywf matrix
    xyf = [xcoord(thisMap_keepPixels), ycoord(thisMap_keepPixels), freqmap']; 
    
    % % implement rotations by 1 degree
    allAlphasRad = linspace(0,2*pi, 361);
    allAlphasRad(end) = [];
    allAlphasRadPI = allAlphasRad(1:180); 
    
    % [manual ~PCA] -> let's try to implement the variance criterion (with N=1, it looks better)
    binsize = 18; %between 15 and 20 seems to be working fine
    for i = 1:length(allAlphasRad)
        alph = allAlphasRad(i);
        xyRot = rotateVector(xyf(:,1:2), alph);
        [~, I] = sort(xyRot(:,1));
        frRotated = xyf(I,3);
        frRotBinned = frRotated(1 : binsize * floor(length(frRotated)/binsize) );
        frRotBinned = reshape(frRotBinned, binsize,[]);
        frRotRemainder = frRotated(binsize * floor(length(frRotated)/binsize) +1 : end);
        variances(i) = mean(cat(2, var(frRotBinned), var(frRotRemainder) ));
        corrValue(i) = corr(frRotated, freqREF);
    end
    variancesPI = mean(reshape(variances, 180,2),2);
    [M,I] = min(variancesPI);
    alpha_rotation = allAlphasRadPI(I);
    
    
    xyRot = rotateVector(xyf(:,1:2), alpha_rotation);
    [~, I] = sort(xyRot(:,1));
    frRotated = xyf(I,3);
    
    % hopefully check if frequency decreases
    if mean(frRotated(1:8)) > mean(frRotated(end-7:end))
        alpha_rotation = alpha_rotation+pi;
        xyRot = rotateVector(xyf(:,1:2), alpha_rotation);
        [~, I] = sort(xyRot(:,1));
        frRotated = xyf(I,3);
    end
    
    rotations(allZetas).keepZeta = keepZeta;
    rotations(allZetas).variancesRotating = variances;
    rotations(allZetas).corrValuesRotating = corrValue;
    rotations(allZetas).frRotated = frRotated;
    
    rotations(allZetas).alphaMinVarRot = alpha_rotation; %THESE ROTATIONS ARE CLOCK-WISE
    correlations(allZetas).regression2Rot = corr(frRotated, freqREF);
%     figure; hold on
%     subplot(1,2,1); scatter(xyf(:,1),xyf(:,2), 'filled'); title(sprintf('fly%d run%d',iCNS.dataTable.fly(keepZeta), iCNS.dataTable.run(keepZeta))); axis image; set(gca, 'YDir', 'reverse')
%     subplot(1,2,2); scatter(xyRot(:,1),xyRot(:,2), 'filled'); title([alphaMinVarRot(allZetas), corrValue(allZetas)]); axis image; set(gca, 'YDir', 'reverse')
    
    %% let's try bootstrap for significance of correlation
    rng shuffle
    L = length(freqmap);
    for p = 1:1e5
        freqperm = freqmap(randperm(L));
        cPerm(p) = corr(freqperm', freqREF);
    end
    correlations(allZetas).cPerm = cPerm;
    correlations(allZetas).pValue = sum(correlations(allZetas).cPerm >= correlations(allZetas).regression2Rot)/length(correlations(allZetas).cPerm);
%     figure; histogram(cPerm)
%     max(cPerm)
    
    %% figure out the total angle of rotation from the registered map 
    % or I can DO it MANUALLY starting from the edge of the recording window
    
end
% save 
save('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/tonotopyIndex_cWED.mat', 'correlations', 'rotations')


%% ipsiTONOT
clear rotations correlations
for allZetas = 1 : length(iCNS.WEDi.keepZetas)
    
    keepZeta = iCNS.WEDi.keepZetas(allZetas);
    keepPixelsZeta = false(size(iCNS.pixel2keep));
    keepPixelsZeta(:,keepZeta) = true;
    keepPixelsZeta = keepPixelsZeta(:);
    keepPixelsZeta = keepPixelsZeta(iCNS.pixel2keep(:)); % indices of all pixel2keep this run % size: 49608,1 - sum: 1126
    
    % apply other masks... and exclude sust clusters' pixels, then convert back
    keepPixels = keepPixelsZeta & iCNS.includedK_i & iCNS.WEDi.TONOTipsi_all; % size: 49608,1 - sum: 764
    
    % iterative approach for excluding some more spurious pixels
    freqmap = CNS.bestFreq_int(keepPixels);                                     % size: 764,1
    excludePixels = find(isnan(freqmap));
    temp = cumsum(keepPixels);
    for i = 1:length(excludePixels)
        excludePixels_i = find(temp == excludePixels(i), 1);
        keepPixels(excludePixels_i) = 0;
    end
    clear temp excludePixels_i
    keepPixels_thisRunPx2Kp = keepPixels(keepPixelsZeta);                       % size: 1126,1  - sum: 764
    freqmap = CNS.bestFreq_int(keepPixels);
    
    
    thismap = iCNS.pixel2keep(:, keepZeta);                                     % size: 5060,1 - sum: 1126    % all px2KP
    thisMap_keepPixels = false(size(thismap));
    thisMap_keepPixels(thismap) = keepPixels_thisRunPx2Kp;                      % size: 5060,1 - sum: 764     % select px2KP
    
    % reshape given maps:
    freqMAPreshaped = nan(size(thismap));
    freqMAPreshaped(thisMap_keepPixels) = freqmap;
    freqMAPreshaped = reshape(freqMAPreshaped, sizeMap);
%     figure; imshow(freqMAPreshaped, []);
    
    keepPixelsMAPreshaped = thisMap_keepPixels;
    keepPixelsMAPreshaped = reshape(keepPixelsMAPreshaped, sizeMap);
%     figure; imshow(keepPixelsMAPreshaped, []); title('no pushpull clusters')

    
    %% manual rotation of points and projection on x axis (x coord)
    freqREF = (sort(freqmap))';
    %take the ywf matrix
    xyf = [xcoord(thisMap_keepPixels), ycoord(thisMap_keepPixels), freqmap']; 
    
    % % implement rotations by 1 degree
    allAlphasRad = linspace(0,2*pi, 361);
    allAlphasRad(end) = [];
    allAlphasRadPI = allAlphasRad(1:180); 
    
    % [manual ~PCA] -> let's try to implement the variance criterion (with N=1, it looks better)
    binsize = 18; %between 15 and 20 seems to be working fine
    for i = 1:length(allAlphasRad)
        alph = allAlphasRad(i);
        xyRot = rotateVector(xyf(:,1:2), alph);
        [~, I] = sort(xyRot(:,1));
        frRotated = xyf(I,3);
        frRotBinned = frRotated(1 : binsize * floor(length(frRotated)/binsize) );
        frRotBinned = reshape(frRotBinned, binsize,[]);
        frRotRemainder = frRotated(binsize * floor(length(frRotated)/binsize) +1 : end);
        variances(i) = mean(cat(2, var(frRotBinned), var(frRotRemainder) ));
        corrValue(i) = corr(frRotated, freqREF);
    end
    variancesPI = mean(reshape(variances, 180,2),2);
    [M,I] = min(variancesPI);
    alpha_rotation1 = allAlphasRadPI(I);
    
    
    xyRot1 = rotateVector(xyf(:,1:2), alpha_rotation1);
    [~, I] = sort(xyRot1(:,1));
    frRotated1 = xyf(I,3);
    tempCorrV1 = corr(frRotated1, freqREF);
    
    alpha_rotation2 = alpha_rotation1+pi;
    xyRot2 = rotateVector(xyf(:,1:2), alpha_rotation2);
    [~, I] = sort(xyRot2(:,1));
    frRotated2 = xyf(I,3);
    tempCorrV2 = corr(frRotated2, freqREF);
        
    if tempCorrV1 > tempCorrV2
        xyRot = xyRot1;
        alpha_rotation = alpha_rotation1;
        frRotated = frRotated1;
        changeCorrelationsPI = (tempCorrV1 - tempCorrV2)/tempCorrV1;
    else
        xyRot = xyRot2;
        alpha_rotation = alpha_rotation2;
        frRotated = frRotated2;
        changeCorrelationsPI = (tempCorrV2 - tempCorrV1)/tempCorrV2;
    end
        
    
    rotations(allZetas).keepZeta = keepZeta;
    rotations(allZetas).variancesRotating = variances;
    rotations(allZetas).corrValuesRotating = corrValue;
    rotations(allZetas).frRotated = frRotated;
    rotations(allZetas).changeCorrelationsPI = changeCorrelationsPI;
    rotations(allZetas).absdiffMaxRotCorrelations = abs(max(corrValue) - corr(frRotated, freqREF));
    
    
    rotations(allZetas).alphaMinVarRot = alpha_rotation; %THESE ROTATIONS ARE CLOCK-WISE
    correlations(allZetas).regression2Rot = corr(frRotated, freqREF);
    figure; hold on
    subplot(1,2,1); scatter(xyf(:,1),xyf(:,2), 'filled'); title(sprintf('fly%d run%d',iCNS.dataTable.fly(keepZeta), iCNS.dataTable.run(keepZeta))); axis image; set(gca, 'YDir', 'reverse', 'FontSize',7)
    subplot(1,2,2); scatter(xyRot(:,1),xyRot(:,2), 'filled'); title([rotations(allZetas).alphaMinVarRot, correlations(allZetas).regression2Rot]); axis image; set(gca, 'YDir', 'reverse', 'FontSize',7)
    
    %% let's try bootstrap for significance of correlation
    rng shuffle
    L = length(freqmap);
    for p = 1:1e5
        freqperm = freqmap(randperm(L));
        cPerm(p) = corr(freqperm', freqREF);
    end
    correlations(allZetas).cPerm = cPerm;
    correlations(allZetas).pValue = sum(correlations(allZetas).cPerm >= correlations(allZetas).regression2Rot)/length(correlations(allZetas).cPerm);
%     figure; histogram(cPerm)
%     max(cPerm)
    
    %% figure out the total angle of rotation from the registered map 
    % or I can DO it MANUALLY starting from the edge of the recording window
    
end
% save 
save('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/tonotopyIndex_iWED.mat', 'correlations', 'rotations')

% I had to manually change the direction of rotation for fly 125_02


%% AMMC
close all
clear rotations correlations
for allZetas = 1: length(iCNS.AMMC.keepZetas)
    
    keepZeta = iCNS.AMMC.keepZetas(allZetas);
    keepPixelsZeta = false(size(iCNS.pixel2keep));
    keepPixelsZeta(:,keepZeta) = true;
    keepPixelsZeta = keepPixelsZeta(:);
    keepPixelsZeta = keepPixelsZeta(iCNS.pixel2keep(:)); % indices of all pixel2keep this run % size: 49608,1 - sum: 1126
    
    % apply other masks... and exclude sust clusters' pixels, then convert back
    keepPixels = keepPixelsZeta & iCNS.includedK_i & iCNS.AMMC.AMMCipsi_all; % size: 49608,1 - sum: 764
    
    % iterative approach for excluding some more spurious pixels
    freqmap = CNS.bestFreq_int(keepPixels);                                     % size: 764,1
    excludePixels = find(isnan(freqmap));
    temp = cumsum(keepPixels);
    for i = 1:length(excludePixels)
        excludePixels_i = find(temp == excludePixels(i), 1);
        keepPixels(excludePixels_i) = 0;
    end
    clear temp excludePixels_i
    keepPixels_thisRunPx2Kp = keepPixels(keepPixelsZeta);                       % size: 1126,1  - sum: 764
    freqmap = CNS.bestFreq_int(keepPixels);
    
    
    thismap = iCNS.pixel2keep(:, keepZeta);                                     % size: 5060,1 - sum: 1126    % all px2KP
    thisMap_keepPixels = false(size(thismap));
    thisMap_keepPixels(thismap) = keepPixels_thisRunPx2Kp;                      % size: 5060,1 - sum: 764     % select px2KP
    
    % reshape given maps:
    freqMAPreshaped = nan(size(thismap));
    freqMAPreshaped(thisMap_keepPixels) = freqmap;
    freqMAPreshaped = reshape(freqMAPreshaped, sizeMap);
%     figure; imshow(freqMAPreshaped, []);
    
    keepPixelsMAPreshaped = thisMap_keepPixels;
    keepPixelsMAPreshaped = reshape(keepPixelsMAPreshaped, sizeMap);
%     figure; imshow(keepPixelsMAPreshaped, []); title('no pushpull clusters')

    
    %% manual rotation of points and projection on x axis (x coord)
    freqREF = (sort(freqmap))';
    %take the ywf matrix
    xyf = [xcoord(thisMap_keepPixels), ycoord(thisMap_keepPixels), freqmap']; 
    
    % % implement rotations by 1 degree
    allAlphasRad = linspace(0,2*pi, 361);
    allAlphasRad(end) = [];
    allAlphasRadPI = allAlphasRad(1:180); 
    
    % [manual ~PCA] -> let's try to implement the variance criterion (with N=1, it looks better)
    binsize = 18; %between 15 and 20 seems to be working fine
    for i = 1:length(allAlphasRad)
        alph = allAlphasRad(i);
        xyRot = rotateVector(xyf(:,1:2), alph);
        [~, I] = sort(xyRot(:,1));
        frRotated = xyf(I,3);
        frRotBinned = frRotated(1 : binsize * floor(length(frRotated)/binsize) );
        frRotBinned = reshape(frRotBinned, binsize,[]);
        frRotRemainder = frRotated(binsize * floor(length(frRotated)/binsize) +1 : end);
        variances(i) = mean(cat(2, var(frRotBinned), var(frRotRemainder) ));
        corrValue(i) = corr(frRotated, freqREF);
    end
    variancesPI = mean(reshape(variances, 180,2),2);
    [M,I] = min(variancesPI);
    alpha_rotation1 = allAlphasRadPI(I);
    
    
    xyRot1 = rotateVector(xyf(:,1:2), alpha_rotation1);
    [~, I] = sort(xyRot1(:,1));
    frRotated1 = xyf(I,3);
    tempCorrV1 = corr(frRotated1, freqREF);
    
    alpha_rotation2 = alpha_rotation1+pi;
    xyRot2 = rotateVector(xyf(:,1:2), alpha_rotation2);
    [~, I] = sort(xyRot2(:,1));
    frRotated2 = xyf(I,3);
    tempCorrV2 = corr(frRotated2, freqREF);
        
    if tempCorrV1 > tempCorrV2
        xyRot = xyRot1;
        alpha_rotation = alpha_rotation1;
        frRotated = frRotated1;
        changeCorrelationsPI = (tempCorrV1 - tempCorrV2)/tempCorrV1;
    else
        xyRot = xyRot2;
        alpha_rotation = alpha_rotation2;
        frRotated = frRotated2;
        changeCorrelationsPI = (tempCorrV2 - tempCorrV1)/tempCorrV2;
    end
        
    
    rotations(allZetas).keepZeta = keepZeta;
    rotations(allZetas).variancesRotating = variances;
    rotations(allZetas).corrValuesRotating = corrValue;
    rotations(allZetas).frRotated = frRotated;
    rotations(allZetas).changeCorrelationsPI = changeCorrelationsPI;
    rotations(allZetas).absdiffMaxRotCorrelations = abs(max(corrValue) - corr(frRotated, freqREF));
    
    
    rotations(allZetas).alphaMinVarRot = alpha_rotation; %THESE ROTATIONS ARE CLOCK-WISE
    correlations(allZetas).regression2Rot = corr(frRotated, freqREF);
    figure; hold on
    subplot(1,2,1); scatter(xyf(:,1),xyf(:,2), 'filled'); title(sprintf('fly%d run%d',iCNS.dataTable.fly(keepZeta), iCNS.dataTable.run(keepZeta))); axis image; set(gca, 'YDir', 'reverse', 'FontSize',7)
    subplot(1,2,2); scatter(xyRot(:,1),xyRot(:,2), 'filled'); title([rotations(allZetas).alphaMinVarRot, correlations(allZetas).regression2Rot]); axis image; set(gca, 'YDir', 'reverse', 'FontSize',7)
    
    %% let's try bootstrap for significance of correlation
    rng shuffle
    L = length(freqmap);
    for p = 1:1e5
        freqperm = freqmap(randperm(L));
        cPerm(p) = corr(freqperm', freqREF);
    end
    correlations(allZetas).cPerm = cPerm;
    correlations(allZetas).pValue = sum(correlations(allZetas).cPerm >= correlations(allZetas).regression2Rot)/length(correlations(allZetas).cPerm);
%     figure; histogram(cPerm)
%     max(cPerm)
    
    %% figure out the total angle of rotation from the registered map 
    % or I can DO it MANUALLY starting from the edge of the recording window
    
end
% save 
save('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/tonotopyIndex_AMMC.mat', 'correlations', 'rotations')

% need to check maps (if they got the right 180deg rotation) before using
% them




%% JONS
close all
clear rotations correlations
for allZetas = 1: length(iJON.PAN.keepZetas)
    
    keepZeta = iJON.PAN.keepZetas(allZetas);
    keepPixelsZeta = false(size(iJON.pixel2keep));
    keepPixelsZeta(:,keepZeta) = true;
    keepPixelsZeta = keepPixelsZeta(:);
    keepPixelsZeta = keepPixelsZeta(iJON.pixel2keep(:)); % indices of all pixel2keep this run % size: 49608,1 - sum: 1126
    
    % apply other masks... and exclude sust clusters' pixels, then convert back
    keepPixels = keepPixelsZeta & iJON.includedK_i; % size: 49608,1 - sum: 764
    
    % iterative approach for excluding some more spurious pixels
    freqmap = JON.bestFreq_int(keepPixels);                                     % size: 764,1
    excludePixels = find(isnan(freqmap));
    temp = cumsum(keepPixels);
    for i = 1:length(excludePixels)
        excludePixels_i = find(temp == excludePixels(i), 1);
        keepPixels(excludePixels_i) = 0;
    end
    clear temp excludePixels_i
    keepPixels_thisRunPx2Kp = keepPixels(keepPixelsZeta);                       % size: 1126,1  - sum: 764
    freqmap = JON.bestFreq_int(keepPixels);
    
    
    thismap = iJON.pixel2keep(:, keepZeta);                                     % size: 5060,1 - sum: 1126    % all px2KP
    thisMap_keepPixels = false(size(thismap));
    thisMap_keepPixels(thismap) = keepPixels_thisRunPx2Kp;                      % size: 5060,1 - sum: 764     % select px2KP
    
    % reshape given maps:
    freqMAPreshaped = nan(size(thismap));
    freqMAPreshaped(thisMap_keepPixels) = freqmap;
    freqMAPreshaped = reshape(freqMAPreshaped, sizeMap);
%     figure; imshow(freqMAPreshaped, []);
    
    keepPixelsMAPreshaped = thisMap_keepPixels;
    keepPixelsMAPreshaped = reshape(keepPixelsMAPreshaped, sizeMap);
%     figure; imshow(keepPixelsMAPreshaped, []); title('no pushpull clusters')

    
    %% manual rotation of points and projection on x axis (x coord)
    freqREF = (sort(freqmap))';
    %take the ywf matrix
    xyf = [xcoord(thisMap_keepPixels), ycoord(thisMap_keepPixels), freqmap']; 
    
    % % implement rotations by 1 degree
    allAlphasRad = linspace(0,2*pi, 361);
    allAlphasRad(end) = [];
    allAlphasRadPI = allAlphasRad(1:180); 
    
    % [manual ~PCA] -> let's try to implement the variance criterion (with N=1, it looks better)
    binsize = 18; %between 15 and 20 seems to be working fine
    for i = 1:length(allAlphasRad)
        alph = allAlphasRad(i);
        xyRot = rotateVector(xyf(:,1:2), alph);
        [~, I] = sort(xyRot(:,1));
        frRotated = xyf(I,3);
        frRotBinned = frRotated(1 : binsize * floor(length(frRotated)/binsize) );
        frRotBinned = reshape(frRotBinned, binsize,[]);
        frRotRemainder = frRotated(binsize * floor(length(frRotated)/binsize) +1 : end);
        variances(i) = mean(cat(2, var(frRotBinned), var(frRotRemainder) ));
        corrValue(i) = corr(frRotated, freqREF);
    end
    variancesPI = mean(reshape(variances, 180,2),2);
    [M,I] = min(variancesPI);
    alpha_rotation1 = allAlphasRadPI(I);
    
    
    xyRot1 = rotateVector(xyf(:,1:2), alpha_rotation1);
    [~, I] = sort(xyRot1(:,1));
    frRotated1 = xyf(I,3);
    tempCorrV1 = corr(frRotated1, freqREF);
    
    alpha_rotation2 = alpha_rotation1+pi;
    xyRot2 = rotateVector(xyf(:,1:2), alpha_rotation2);
    [~, I] = sort(xyRot2(:,1));
    frRotated2 = xyf(I,3);
    tempCorrV2 = corr(frRotated2, freqREF);
        
    if tempCorrV1 > tempCorrV2
        xyRot = xyRot1;
        alpha_rotation = alpha_rotation1;
        frRotated = frRotated1;
        changeCorrelationsPI = (tempCorrV1 - tempCorrV2)/tempCorrV1;
    else
        xyRot = xyRot2;
        alpha_rotation = alpha_rotation2;
        frRotated = frRotated2;
        changeCorrelationsPI = (tempCorrV2 - tempCorrV1)/tempCorrV2;
    end
        
    
    rotations(allZetas).keepZeta = keepZeta;
    rotations(allZetas).variancesRotating = variances;
    rotations(allZetas).corrValuesRotating = corrValue;
    rotations(allZetas).frRotated = frRotated;
    rotations(allZetas).changeCorrelationsPI = changeCorrelationsPI;
    rotations(allZetas).absdiffMaxRotCorrelations = abs(max(corrValue) - corr(frRotated, freqREF));
    
    
    rotations(allZetas).alphaMinVarRot = alpha_rotation; %THESE ROTATIONS ARE CLOCK-WISE
    correlations(allZetas).regression2Rot = corr(frRotated, freqREF);
    figure; hold on
    subplot(1,2,1); scatter(xyf(:,1),xyf(:,2), 'filled'); title(sprintf('fly%d run%d',iJON.dataTable.fly(keepZeta), iJON.dataTable.run(keepZeta))); axis image; set(gca, 'YDir', 'reverse', 'FontSize',7)
    subplot(1,2,2); scatter(xyRot(:,1),xyRot(:,2), 'filled'); title([rotations(allZetas).alphaMinVarRot, correlations(allZetas).regression2Rot]); axis image; set(gca, 'YDir', 'reverse', 'FontSize',7)
    
    %% let's try bootstrap for significance of correlation
    rng shuffle
    L = length(freqmap);
    for p = 1:1e5
        freqperm = freqmap(randperm(L));
        cPerm(p) = corr(freqperm', freqREF);
    end
    correlations(allZetas).cPerm = cPerm;
    correlations(allZetas).pValue = sum(correlations(allZetas).cPerm >= correlations(allZetas).regression2Rot)/length(correlations(allZetas).cPerm);
%     figure; histogram(cPerm)
%     max(cPerm)
    
    %% figure out the total angle of rotation from the registered map 
    % or I can DO it MANUALLY starting from the edge of the recording window
    
end
% save 
save('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix/tonotopyIndex_JONS.mat', 'correlations', 'rotations')

% need to check maps (if they got the right 180deg rotation) before using
% them





%% strip chart

%% cluster by cluster anatomical axis: need to average and still sucks
% keepClustersCNS = 3:19;
% 
% xcoord = repmat(1:sizeMap(2), sizeMap(1), 1);
% ycoord = repmat((1:sizeMap(1))', 1, sizeMap(2));
% xcoord = xcoord(:);
% ycoord = ycoord(:);
% 
% thismap = iCNS.pixel2keep(:, keepZeta);
% cf = [];
% for k = keepClustersCNS
%     keepPixels = keepPixelsZeta & iCNS.includedK_i & iCNS.WEDc.TONOTcontra_Ki(:,k); % size: 49608,1 
%     if sum(keepPixels)>20
%         keepPixels_thisRunPx2Kp = keepPixels(keepPixelsZeta);
%         thisMap_keepPixels = false(size(thismap));
%         thisMap_keepPixels(thismap) = keepPixels_thisRunPx2Kp;
%         xy = [xcoord(thisMap_keepPixels), ycoord(thisMap_keepPixels)];
%         xy = bsxfun(@minus, xy, mean(xy));
%         [cf(k).coeff, scores, latent] = pca(xy);
%         figure; scatter(xy(:,1), xy(:,2)); hold on; title(k)
%         line(15*[0,cf(k).coeff(1,1)], 15*[0,cf(k).coeff(2,1)], 'LineWidth', 4, 'Color', 'r')
%         line(15*[0,cf(k).coeff(1,2)], 15*[0,cf(k).coeff(2,2)], 'LineWidth', 3, 'Color', 'g')
%         ax = gca;
%         ax.YAxis.Direction = 'reverse';
%         axis image
%     else
%         cf(k).coeff = [];
%     end
% end
% 
% figure; hold on;
% for k = keepClustersCNS
%     if ~isempty(cf(k).coeff)
%         line(15*[0,cf(k).coeff(1,1)], 15*[0,cf(k).coeff(2,1)], 'LineWidth', 2, 'Color', 'r')
%     end
% end
% ax = gca;
% ax.YAxis.Direction = 'reverse';
% axis image
% 
% 
% avgCoef1 = cat(3,cf.coeff);
% avgCoef1 = mean(squeeze(avgCoef1(:,1,:)),2);
% line(15*[0,avgCoef1(1)], 15*[0, avgCoef1(2)],'LineWidth', 2, 'Color', 'g')
% avgCoef2 = cat(3,cf.coeff);
% avgCoef2 = mean(squeeze(avgCoef2(:,2,:)),2);
% line(15*[0,avgCoef2(1)], 15*[0, avgCoef2(2)],'LineWidth', 4, 'Color', 'b')
% % Now, these averages are not good because they are not orthogonal anymore
% % fix it


%
% 
% 
% 
% %% nope
% [Gmag,Gdir] = imgradient(freqMAPreshaped);
% figure; imagesc(Gdir)



% freq01 = mat2gray(freqmap);
% 
% xcoord = repmat(1:sizeMap(2), sizeMap(1), 1);
% ycoord = repmat((1:sizeMap(1))', 1, sizeMap(2));
% xcoord = xcoord(:);
% ycoord = ycoord(:);
% xyf = [xcoord(thisMap_keepPixels), ycoord(thisMap_keepPixels), freq01'];
% figure; scatter3(xyf(:,1), xyf(:,2), xyf(:,3))
% 
% % try simple pca (do this on single-clusters binary maps)
% mucent = mean(xyf);
% xyf = bsxfun(@minus, xyf, mucent);
% [coeff, scores, latent] = pca(xyf);
% 
% 
% hold on
% line(40*[0,coeff(1,1)], 40*[0,coeff(2,1)], 1*[0,coeff(3,1)], 'LineWidth', 4, 'Color', 'r')
% line(40*[0,coeff(1,2)], 40*[0,coeff(2,2)], 1*[0,coeff(3,2)], 'LineWidth', 3, 'Color', 'g')
% line(40*[0,coeff(1,3)], 40*[0,coeff(2,3)], 1*[0,coeff(3,3)], 'LineWidth', 2, 'Color', 'y')
% 
% 
% [coeff, scores, latent] = pca(xyf(:,1:2));
% figure; scatter(xyf(:,1), xyf(:,2))
% hold on
% line(40*[0,coeff(1,1)], 40*[0,coeff(2,1)], 'LineWidth', 4, 'Color', 'r')
% line(40*[0,coeff(1,2)], 40*[0,coeff(2,2)], 'LineWidth', 3, 'Color', 'g')
% line(40*[0,coeff(1,3)], 40*[0,coeff(2,3)], 1*[0,coeff(3,3)], 'LineWidth', 2, 'Color', 'y')

% %% try PLSR
% % predictors = xyf(:,1:2);
% % variable = xyf(:,3);
% 
% predictors = zscore(xyf(:,1:2));
% variable = zscore(xyf(:,3));
% 
% [XL,yl,XS,YS,beta,PCTVAR] = plsregress(predictors,variable,1);
% 
% figure; scatter(predictors(:,1), predictors(:,2))
% hold on
% line(
% % concatenate 










