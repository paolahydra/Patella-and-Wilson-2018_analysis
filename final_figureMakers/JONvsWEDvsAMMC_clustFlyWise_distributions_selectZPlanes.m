% tuning data previously assembled in JON_pixelwise_distributions and WEDvsAMMC_pixelwise_distributions
% now reassmbled in 
% edit WEDvsAMMC_pixelwise_distributions_selectZPlanesAndRegions.m


%% you want to recalculate these based on the selected runs/regions
% load JON
JON = load('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix/ftcJON_clusterTuningFly.mat');
% load AMMC + WED
CNS = load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/ftc_AMMCi_TONOTi_TONOTc_clusterTuningFly_mean.mat');
%CNS.tcAMMC_clusterTuningFly % 19x19x11
% save('ftc_AMMCi_TONOTi_TONOTc_clusterTuningFly.mat', 'tcAMMCipsi_clusterTuningFly', 'tcTONOTipsi_clusterTuningFly', 'tcTONOTcontra_clusterTuningFly', ...
%     'tcAMMCipsi_sumPixels_clusterFly', 'tcTONOTipsi_sumPixels_clusterFly', 'tcTONOTcontra_sumPixels_clusterFly') % all 19 clusters

%% keep loading 

load('/Users/galileo/Dropbox (HMS)/Data/chunks_stimfinal.mat')

CNS.clusterFolder = '/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix';
JON.clusterFolder = '/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix';
CNS.datalist = load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/datalist_WEDAMMC_piezo.mat');

CNS.clusters = load(fullfile(CNS.clusterFolder, 'klust_maxClust19.mat'));
JON.clusters = load(fullfile(JON.clusterFolder, 'klust_maxClust10.mat'));
JON.excludeKs = [8,9,10]; % brown , pull, push
CNS.excludeKs = [1,2]; % pull, push
JON.includedK_i = ~ismember(JON.clusters.Tparent, JON.excludeKs);
CNS.includedK_i = ~ismember(CNS.clusters.Tparent, CNS.excludeKs);

%
[xvar, ~, ~] = unique([chunks.pips.carrierFreqLevels; chunks.tones.carrierFreqLevels(1)]);

% we are not interpolating here, and I still think it is not such a good idea.
onset = 4;
offset = 600; %cropping before, makes things worse
cutoff = 0.5;
x2 = logspace(log10(onset), log10(offset), 2000);
x3 = logspace(log10(onset), log10(offset), 120);
x2Zero = cat(2, 0, x2);   

figure; hold on
plot(linspace(1,length(x3), length(xvar)), xvar, 'xr')
plot(x3, '+g')
xlabel('sampling points')
ylabel('frequency (Hz)')
legend('actual', 'log-interpolated')
export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/x3Interpolation_frequencies.eps')

%% exclude tdt+ flies from CNS
% % pensavo di aver escluso fly 126 once and for all. Invece, sta ancora nel
% % dataset (155 pixels)
% tdtFlies = unique(CNS.datalist.dataTable.fly(CNS.datalist.dataTable.tdTomato==1));
% allCNSflies = unique(CNS.datalist.dataTable.fly, 'stable');
% keepCNSflies = ~(ismember(allCNSflies, tdtFlies));
% 
% figure;
% hold on
% title('average TC AMMC gal80 flies')
% avgAMMCGal80includedKs = CNS.tcAMMC_clusterTuningFly;
% avgAMMCGal80includedKs(CNS.excludeKs,:,:) = [];
% avgAMMCGal80includedKs = avgAMMCGal80includedKs(:,:, keepCNSflies);
% avgAMMCGal80includedKs = squeeze(nanmean(avgAMMCGal80includedKs,1));
% plot(avgAMMCGal80includedKs, 'LineWidth', 2)
% legend(num2str(allCNSflies(keepCNSflies)))
% export
% figure;
% hold on
% title('average TC AMMC tdTomato flies')
% avgAMMCtdTincludedKs = CNS.tcAMMC_clusterTuningFly;
% avgAMMCtdTincludedKs(CNS.excludeKs,:,:) = [];
% avgAMMCtdTincludedKs = avgAMMCtdTincludedKs(:,:, ~keepCNSflies);
% avgAMMCtdTincludedKs = squeeze(nanmean(avgAMMCtdTincludedKs,1));
% plot(avgAMMCtdTincludedKs, 'LineWidth', 2)
% legend(num2str(allCNSflies(~keepCNSflies)))

%% calculate all spread measures within each region ( for each tuning curve )
for region = 1:3 %JON, AMMC, WED
    switch region
        case 1 %JON
            ttc = JON.tcJON_clusterTuningFly;
            ttc(JON.excludeKs,:,:) = [];
        case 2 %AMMCi
            ttc = CNS.tcAMMCipsi_clusterTuningFly;
            ttc(CNS.excludeKs,:,:) = [];
        case 3 %TONOTi
            ttc = CNS.tcTONOTipsi_clusterTuningFly;
            ttc(CNS.excludeKs,:,:) = [];
        case 4 %TONOTc
            ttc = CNS.tcTONOTcontra_clusterTuningFly;
            ttc(CNS.excludeKs,:,:) = [];
        case 5
            ttc = CNS.tcTONOTboth_clusterTuningFly;
            ttc(CNS.excludeKs,:,:) = [];
    end

    %% calculate log-interpolate sigma and mu
%     sigmaREG = nan(size(ttc, 1), size(ttc,3)); %cluster by fly
%     muREG = nan(size(ttc, 1), size(ttc,3)); %cluster by fly
%     for Ki = 1 : size(ttc, 1)
%         for f = 1 : size(ttc, 3)
%             tc = interp1(xvar, ttc(Ki, :, f), x2);
%             [sigmaREG(Ki, f), muREG(Ki, f)] = gaussfit( 1:2000, tc); %, sigma0(Ki), mu0(Ki));
%             %         [sigmaWED(Ki, f), muWED(Ki, f)] = gaussfit( xvar, ttc(Ki, :, f), sigma0(Ki), mu0(Ki));
%         end
%     end
%     
%     muREGHz = round(muREG);
%     for i = 1:numel(muREGHz)
%         if ~isnan(muREGHz(i))
%             muREGHz(i) = x2(muREGHz(i));
%         end
%     end
%     
%     %% calculate rmax/rmean for each tuning curve
%     RmeanRmax = nan(size(ttc, 1), size(ttc,3));
%     for Ki = 1 : size(ttc, 1)
%         for f = 1 : size(ttc, 3)
%             RmeanRmax(Ki, f) = nanmean(ttc(Ki, :, f)) / max(ttc(Ki, :, f));
%         end
%     end
%     
    %% calculate interpolated rmax/rmean for each tuning curve - DO
    RmeanRmax_int = nan(size(ttc, 1), size(ttc,3));
    for Ki = 1 : size(ttc, 1)
        for f = 1 : size(ttc, 3)
            tc = interp1(xvar, ttc(Ki, :, f), x3);
            RmeanRmax_int(Ki, f) = nanmean(tc) / max(tc);
        end
    end
    
    %% calculate lifetime sparseness for each tuning curve
%     nSt = length(xvar);
%     lfts = nan(size(ttc, 1), size(ttc,3));
%     for Ki = 1 : size(ttc, 1)
%         for f = 1 : size(ttc, 3)
%             A = ttc(Ki, :, f);
%             num = sum((A/nSt))^2;
%             den = sum(A.^2/nSt);
%             if sum(A) == 0
%                 lfts(Ki, f) = 1;
%             else
%                 lfts(Ki, f) = 1 - num/den;
%             end
%         end
%     end
    
    %% calculate interpolated lifetime sparseness for each tuning curve - DO
    nSt = length(x3);
    lfts_int = nan(size(ttc, 1), size(ttc,3));
    for Ki = 1 : size(ttc, 1)
        for f = 1 : size(ttc, 3)
            A = interp1(xvar, ttc(Ki, :, f), x3);
            num = sum((A/nSt))^2;
            den = sum(A.^2/nSt);
            if sum(A) == 0
                lfts_int(Ki, f) = 1;
            else
                lfts_int(Ki, f) = 1 - num/den;
            end
        end
    end
    
    %% store region-specific results
    switch region
        case 1 %JON
%             JON.sigma = sigmaREG;
%             JON.mu = muREG;
%             JON.muHz = muREGHz;
%             JON.RmeanRmax = RmeanRmax;
%             JON.lfts = lfts;
            JON.RmeanRmax_int = RmeanRmax_int;
            JON.lfts_int = lfts_int;
        case 2 %AMMC
%             AMMC.sigma = sigmaREG;
%             AMMC.mu = muREG;
%             AMMC.muHz = muREGHz;
%             AMMC.RmeanRmax = RmeanRmax;
%             AMMC.lfts = lfts;
            AMMC.RmeanRmax_int = RmeanRmax_int;
            AMMC.lfts_int = lfts_int;
        case 3 %WED
%             WED.sigma = sigmaREG;
%             WED.mu = muREG;
%             WED.muHz = muREGHz;
%             WED.RmeanRmax = RmeanRmax;
%             WED.lfts = lfts;
            WED.RmeanRmax_int = RmeanRmax_int;
            WED.lfts_int = lfts_int;
        case 4 %WEDc
%             TNC.sigma = sigmaREG;
%             TNC.mu = muREG;
%             TNC.muHz = muREGHz;
%             TNC.RmeanRmax = RmeanRmax;
%             TNC.lfts = lfts;
            TNC.RmeanRmax_int = RmeanRmax_int;
            TNC.lfts_int = lfts_int;
        case 5 %TONOTboth
%             ALT.sigma = sigmaREG;
%             ALT.mu = muREG;
%             ALT.muHz = muREGHz;
%             ALT.RmeanRmax = RmeanRmax;
%             ALT.lfts = lfts;
            ALT.RmeanRmax_int = RmeanRmax_int;
            ALT.lfts_int = lfts_int;
    end
end

%% scatter max and lfts for all regions, log and linear
close all
x3 = logspace(log10(onset), log10(offset), 200);
xl = linspace(onset, offset, 200);
[hax, pos] = tight_subplot(2,3);

% %% JON
% ttc = nanmean(JON.tcJON_clusterTuningFly,3);
% ttc(JON.excludeKs,:) = [];
% % peaks, by clusters
% for Ki = 1 : size(ttc, 1)
%     [~, muI(Ki)] = max(ttc(Ki, :));
% end
% mu0 = xvar(muI);
% 
% % calculate interpolated lifetime sparseness for each tuning curve - linear
% % ttcN = bsxfun(@rdivide, ttc, max(ttc, [], 2) ); % NO CHANGE
% nSt = length(xl);
% lfts_lin = nan(1, size(ttc, 1));
% for Ki = 1 : size(ttc, 1)
%         A = interp1(xvar, ttc(Ki, :), xl);
%         num = sum((A/nSt))^2;
%         den = sum(A.^2/nSt);
%         if sum(A) == 0
%             lfts_lin(Ki) = 1;
%         else
%             lfts_lin(Ki) = 1 - num/den;
%         end
% end
% % calculate interpolated lifetime sparseness for each tuning curve - log
% nSt = length(x3);
% lfts_log = nan(1, size(ttc, 1));
% for Ki = 1 : size(ttc, 1)
%         A = interp1(xvar, ttc(Ki, :), x3);
%         num = sum((A/nSt))^2;
%         den = sum(A.^2/nSt);
%         if sum(A) == 0
%             lfts_log(Ki) = 1;
%         else
%             lfts_log(Ki) = 1 - num/den;
%         end
% end
% %
% scatter(hax(1), mu0, lfts_log, 'filled'), axis image
% hax(1).XLim = [4,600];
% hax(1).XScale = 'log';
% hax(1).XAxisLocation = 'top';
% 
% hold on
% 
% scatter(hax(4), mu0, lfts_lin, 'filled'), axis image
% hax(4).XLim = [4,600];
% 
% %% AMMC
% ttc = nanmean(CNS.tcAMMCipsi_clusterTuningFly,3);
% ttc(CNS.excludeKs,:) = [];
% % peaks, by clusters
% for Ki = 1 : size(ttc, 1)
%     [~, muI(Ki)] = max(ttc(Ki, :));
% end
% mu0 = xvar(muI);
% 
% % calculate interpolated lifetime sparseness for each tuning curve - linear
% % ttcN = bsxfun(@rdivide, ttc, max(ttc, [], 2) ); % NO CHANGE
% nSt = length(xl);
% lfts_lin = nan(1, size(ttc, 1));
% for Ki = 1 : size(ttc, 1)
%         A = interp1(xvar, ttc(Ki, :), xl);
%         num = sum((A/nSt))^2;
%         den = sum(A.^2/nSt);
%         if sum(A) == 0
%             lfts_lin(Ki) = 1;
%         else
%             lfts_lin(Ki) = 1 - num/den;
%         end
% end
% % calculate interpolated lifetime sparseness for each tuning curve - log
% nSt = length(x3);
% lfts_log = nan(1, size(ttc, 1));
% for Ki = 1 : size(ttc, 1)
%         A = interp1(xvar, ttc(Ki, :), x3);
%         num = sum((A/nSt))^2;
%         den = sum(A.^2/nSt);
%         if sum(A) == 0
%             lfts_log(Ki) = 1;
%         else
%             lfts_log(Ki) = 1 - num/den;
%         end
% end
% %
% scatter(hax(2), mu0, lfts_log, 'filled'), axis image
% hax(2).XLim = [4,600];
% hax(2).XScale = 'log';
% hax(2).XAxisLocation = 'top';
% 
% 
% scatter(hax(5), mu0, lfts_lin, 'filled'), axis image
% hax(5).XLim = [4,600];
% 
% %% WED
% ttc = nanmean(CNS.tcWED_clusterTuningFly,3);
% ttc(CNS.excludeKs,:) = [];
% % peaks, by clusters
% for Ki = 1 : size(ttc, 1)
%     [~, muI(Ki)] = max(ttc(Ki, :));
% end
% mu0 = xvar(muI);
% 
% % calculate interpolated lifetime sparseness for each tuning curve - linear
% % ttcN = bsxfun(@rdivide, ttc, max(ttc, [], 2) ); % NO CHANGE
% nSt = length(xl);
% lfts_lin = nan(1, size(ttc, 1));
% for Ki = 1 : size(ttc, 1)
%         A = interp1(xvar, ttc(Ki, :), xl);
%         num = sum((A/nSt))^2;
%         den = sum(A.^2/nSt);
%         if sum(A) == 0
%             lfts_lin(Ki) = 1;
%         else
%             lfts_lin(Ki) = 1 - num/den;
%         end
% end
% % calculate interpolated lifetime sparseness for each tuning curve - log
% nSt = length(x3);
% lfts_log = nan(1, size(ttc, 1));
% for Ki = 1 : size(ttc, 1)
%         A = interp1(xvar, ttc(Ki, :), x3);
%         num = sum((A/nSt))^2;
%         den = sum(A.^2/nSt);
%         if sum(A) == 0
%             lfts_log(Ki) = 1;
%         else
%             lfts_log(Ki) = 1 - num/den;
%         end
% end
% %
% scatter(hax(3), mu0, lfts_log, 'filled'), axis image
% hax(3).XLim = [4,600];
% hax(3).XScale = 'log';
% hax(3).XAxisLocation = 'top';
% 
% 
% scatter(hax(6), mu0, lfts_lin, 'filled'), axis image
% hax(6).XLim = [4,600];
% 
% hax(1).YLim = hax(2).YLim;
% linkaxes(hax,'y')

%% weight (also) JON for cluster prevalence (wathever)
areaJON_K = JON.tcJON_sumPixels_clusterFly;
areaJON_K(JON.excludeKs,:) = [];
areaJON_K = bsxfun(@rdivide, areaJON_K, sum(areaJON_K));


JON.RmeanRmaxINTWeighted = JON.RmeanRmax_int .* areaJON_K;
JON.RmeanRmaxINTWeighted = nansum(JON.RmeanRmaxINTWeighted );

JON.lftsINTWeighted = JON.lfts_int .* areaJON_K;
JON.lftsINTWeighted = nansum(JON.lftsINTWeighted );



%% weight WED and AMMC results by cluster prevalence
% this was assembled in WEDvsAMMC_clusterdistribution_flybyfly.m 
% in fact, no need to split that by fly. Just keep it constant, so that
% there is less noise. Right?

areaAMMC_K = CNS.tcAMMCipsi_sumPixels_clusterFly;
areaAMMC_K(CNS.excludeKs,:) = [];
areaAMMC_K = bsxfun(@rdivide, areaAMMC_K, sum(areaAMMC_K));

areaWED_K = CNS.tcTONOTipsi_sumPixels_clusterFly;
areaWED_K(CNS.excludeKs,:) = [];
areaWED_K = bsxfun(@rdivide, areaWED_K, sum(areaWED_K) );


%
AMMC.RmeanRmaxINTWeighted = AMMC.RmeanRmax_int .* areaAMMC_K;
AMMC.RmeanRmaxINTWeighted = nansum(AMMC.RmeanRmaxINTWeighted );

AMMC.lftsINTWeighted = AMMC.lfts_int .* areaAMMC_K;
AMMC.lftsINTWeighted = nansum(AMMC.lftsINTWeighted );


%
WED.RmeanRmaxINTWeighted = WED.RmeanRmax_int .* areaWED_K;
WED.RmeanRmaxINTWeighted = nansum(WED.RmeanRmaxINTWeighted );

WED.lftsINTWeighted = WED.lfts_int .* areaWED_K;
WED.lftsINTWeighted = nansum(WED.lftsINTWeighted );


%% not used anymore
% areaTNC_K = CNS.tcTONOTcontra_sumPixels_clusterFly;
% areaTNC_K(CNS.excludeKs,:) = [];
% areaTNC_K = bsxfun(@rdivide, areaTNC_K, sum(areaTNC_K) );
% 
% areaALT_K = CNS.tcTONOTboth_sumPixels_clusterFly;
% areaALT_K(CNS.excludeKs,:) = [];
% areaALT_K = bsxfun(@rdivide, areaALT_K, sum(areaALT_K) ); 
% 
% %
% TNC.RmeanRmaxINTWeighted = TNC.RmeanRmax_int .* areaTNC_K;
% TNC.RmeanRmaxINTWeighted = nansum(TNC.RmeanRmaxINTWeighted );
% 
% TNC.lftsINTWeighted = TNC.lfts_int .* areaTNC_K;
% TNC.lftsINTWeighted = nansum(TNC.lftsINTWeighted );
% 
% %
% ALT.RmeanRmaxINTWeighted = ALT.RmeanRmax_int .* areaALT_K;
% ALT.RmeanRmaxINTWeighted = nansum(ALT.RmeanRmaxINTWeighted );
% 
% ALT.lftsINTWeighted = ALT.lfts_int .* areaALT_K;
% ALT.lftsINTWeighted = nansum(ALT.lftsINTWeighted );



%% plot relevant ones - Rmean/Rmx (all for now)
allSpreads = cat(2, JON.RmeanRmaxINTWeighted, AMMC.RmeanRmaxINTWeighted, WED.RmeanRmaxINTWeighted); % ALT.RmeanRmaxINTWeighted, TNC.RmeanRmaxINTWeighted);
origin = cat(2, 1*ones(size(JON.RmeanRmaxINTWeighted)), ...
                2*ones(size(AMMC.RmeanRmaxINTWeighted)), ...
                3*ones(size(WED.RmeanRmaxINTWeighted)) ); %, ...
%                 4*ones(size(ALT.RmeanRmaxINTWeighted)), ...
%                 5*ones(size(TNC.RmeanRmaxINTWeighted)) );

figure; hold on
[H] = notBoxPlot(allSpreads, origin, 'style', 'line', 'markMedian', true);

ax = gca;
ax.YLim = [0,1];
ax.XTickLabel =  {'JON', 'AMMCi', 'WEDi'}; %, 'WEDb', 'WEDc'}; 
ax.YLabel.String = {'R_m_e_a_n/R_m_a_x, weighted by cluster extent - interpolated xs'};
ax.FontSize = 7;

export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/RmeanRmax_selectedZPlanes.eps')


%% plot relevant ones - lifetime sparseness (all for now)
allSpreads = cat(2, JON.lftsINTWeighted, AMMC.lftsINTWeighted, WED.lftsINTWeighted); %, ALT.lftsINTWeighted, TNC.lftsINTWeighted);
origin = cat(2, 1*ones(size(JON.lftsINTWeighted)), ...
                2*ones(size(AMMC.lftsINTWeighted)), ...
                3*ones(size(WED.lftsINTWeighted)) ); %, ...
%                 4*ones(size(ALT.lftsINTWeighted)), ...
%                 5*ones(size(TNC.lftsINTWeighted)) );

figure; hold on
[H] = notBoxPlot(allSpreads, origin, 'style', 'line', 'markMedian', true);

ax = gca;
ax.YLim = [0,1];
ax.XTickLabel =  {'JON', 'AMMCi', 'WEDi'}; %, 'WEDb', 'WEDc'}; 
ax.YLabel.String = {'lifetime sparseness, weighted by cluster extent - interpolated xs'};
ax.FontSize = 7;

export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/lifetimesparseness_selectedZPlanes.eps')

%% establish significance anova1 and tukey post hoc test
clear stats
[p,t,stats] = anova1(allSpreads,origin);
[resultsT] = multcompare(stats,'CType','tukey-kramer');
[resultsB] = multcompare(stats,'CType','bonferroni');


%% establish significance (pairwise ttest)
figure; normplot(WED.lftsINTWeighted); title('iwed')
lillie(3) = lillietest(WED.lftsINTWeighted);
jbt(3) = jbtest(WED.lftsINTWeighted);

figure; normplot(AMMC.lftsINTWeighted); title('ammc')
lillie(2) = lillietest(AMMC.lftsINTWeighted);
jbt(2) = jbtest(AMMC.lftsINTWeighted);

figure; normplot(JON.lftsINTWeighted); title('jon')
lillie(1) = lillietest(JON.lftsINTWeighted);
jbt(1) = jbtest(JON.lftsINTWeighted);

%ok, normal distribution it is.

% check variances too
[var12.h, var12.p] = vartest2(JON.lftsINTWeighted , AMMC.lftsINTWeighted)
[var13.h, var13.p] = vartest2(JON.lftsINTWeighted , WED.lftsINTWeighted)
[var23.h, var23.p] = vartest2(AMMC.lftsINTWeighted , WED.lftsINTWeighted)


% ttest
[tt12.h , tt12.p ] = ttest2(JON.lftsINTWeighted , AMMC.lftsINTWeighted);
[tt13.h , tt13.p ] = ttest2(JON.lftsINTWeighted , WED.lftsINTWeighted);
[tt23.h , tt23.p ] = ttest2(AMMC.lftsINTWeighted , WED.lftsINTWeighted); 



%% need weights for JONs as well. An anatomically underrepresented type will weight less. Are you good with that?
% load('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix/10clusters/clusterHeatMaps/alignedData.mat','singleKmaps')
% load('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/anatomyImages/JON_cropMask.mat') %cropMask
% for z = 1:size(singleKmaps,3)
%     singleKmapsJON(:,:,z,:) = bsxfun( @and, singleKmaps(:,:,z,:), cropMask);
% end
% for Ki = 1:size(singleKmaps,4) 
%     maxJON_K(:,:,Ki) = max(singleKmapsJON(:,:,:,Ki), [], 3);
%     areaJON_K(Ki) = sum(sum(maxJON_K(:,:,Ki)));
% end
% % % sanity check on cluster numbering
% % for Ki = 1:10  
% %     figure; imshow(maxJON_K(:,:,Ki), []); title(sprintf('cluster %d - sumPx: %d', Ki, sum(JON.clusters.Tparent==Ki)))
% % end
% % % ok: they follow this order in responses figures, bottom to up: clusterOrder = [4 1 2 3 7 5 6 8 9 10]; % this is K2use!! bottom-up 
% clear singleKmaps cropMask maxJON_K singleKmapsJON
% areaJON_K(JON.excludeKs) = [];
% areaJON_K = areaJON_K ./ sum(areaJON_K);
% % this leaves areaJON_K


%% try octave fly-cluster-wise
% xOctaves = -3:3;
% xOctavesFreqMultiplier = 2.^xOctaves;
% 
% ttc = JON.tcJON_clusterTuningFly;
% ttc(JON.excludeKs,:,:) = [];
% JON.octDist = nan(size(ttc,1), size(ttc,3), length(xOctaves)); %clusterFlyOctaves
% for Ki = 1 : size(ttc, 1)
%     for f = 1 : size(ttc, 3)
%         [~, muI] = max(ttc(Ki, :, f)); %peak rather than mu
%         mu0 = xvar(muI);
%         JON.octDist(Ki,f,:) = interp1(xvar', ttc(Ki, :, f), mu0.*xOctavesFreqMultiplier);
%     end
% end
% 
% 
% ttc = CNS.tcAMMC_clusterTuningFly;
% ttc(CNS.excludeKs,:,:) = [];
% AMMC.octDist = nan(size(ttc,1), size(ttc,3), length(xOctaves)); %clusterFlyOctaves
% for Ki = 1 : size(ttc, 1)
%     for f = 1 : size(ttc, 3)
%         [~, muI] = max(ttc(Ki, :, f));  %peak rather than mu
%         mu0 = xvar(muI);
%         AMMC.octDist(Ki,f,:) = interp1(xvar', ttc(Ki, :, f), mu0.*xOctavesFreqMultiplier);
%     end
% end
% 
% ttc = CNS.tcWED_clusterTuningFly;
% ttc(CNS.excludeKs,:,:) = [];
% WED.octDist = nan(size(ttc,1), size(ttc,3), length(xOctaves)); %clusterFlyOctaves
% for Ki = 1 : size(ttc, 1)
%     for f = 1 : size(ttc, 3)
%         [~, muI] = max(ttc(Ki, :, f));  %peak rather than mu
%         mu0 = xvar(muI);
%         WED.octDist(Ki,f,:) = interp1(xvar', ttc(Ki, :, f), mu0.*xOctavesFreqMultiplier);
%     end
% end
% 
% %% weigthed average across clusters 
% octDistNull = JON.octDist;
% JON.octavesWeighted = squeeze(nansum(repmat(areaJON_K', 1, size(JON.octDist,2), size(JON.octDist,3)) .* JON.octDist,1)); %flyes-octaves
% JON.octavesCountFlies = squeeze(nanmean(octDistNull));
% JON.octavesWeighted(isnan(JON.octavesCountFlies)) = nan;
% JON.octavesCountFlies = sum(~isnan(JON.octavesCountFlies));
% JON.octavesWeighted = bsxfun(@rdivide, JON.octavesWeighted, JON.octavesWeighted(:,4));
% 
% 
% octDistNull = AMMC.octDist;
% AMMC.octavesWeighted = squeeze(nansum(repmat(areaAMMC_K', 1, size(AMMC.octDist,2), size(AMMC.octDist,3)) .* octDistNull,1)); %flyes-octaves
% AMMC.octavesCountFlies = squeeze(nanmean(octDistNull));
% AMMC.octavesWeighted(isnan(AMMC.octavesCountFlies)) = nan;
% AMMC.octavesCountFlies = sum(~isnan(AMMC.octavesCountFlies));
% AMMC.octavesWeighted = bsxfun(@rdivide, AMMC.octavesWeighted, AMMC.octavesWeighted(:,4));
% 
% 
% octDistNull = WED.octDist;
% WED.octavesWeighted = squeeze(nansum(repmat(areaWED_K', 1, size(WED.octDist,2), size(WED.octDist,3)) .* octDistNull,1)); %flyes-octaves
% WED.octavesCountFlies = squeeze(nanmean(octDistNull));
% WED.octavesWeighted(isnan(WED.octavesCountFlies)) = nan;
% WED.octavesCountFlies = sum(~isnan(WED.octavesCountFlies));
% WED.octavesWeighted = bsxfun(@rdivide, WED.octavesWeighted, WED.octavesWeighted(:,4));
% 
% 
% errFlyJONweighted = nanstd(JON.octavesWeighted)./sqrt(JON.octavesCountFlies);
% errFlyAMMCweighted = nanstd(AMMC.octavesWeighted)./sqrt(AMMC.octavesCountFlies);
% errFlyWEDweighted = nanstd(WED.octavesWeighted)./sqrt(WED.octavesCountFlies);
% 
% figure; hold on
% errorbar(xOctaves, nanmean(JON.octavesWeighted), errFlyJONweighted, 'Color', cols(1,:))
% errorbar(xOctaves, nanmean(AMMC.octavesWeighted), errFlyAMMCweighted, 'Color', cols(2,:))
% errorbar(xOctaves, nanmean(WED.octavesWeighted), errFlyWEDweighted, 'Color', cols(3,:))
% 
% legend('JON', 'AMMC', 'WED')
% xlabel('octaves')
% set(gca, 'XTick', xOctaves)
% ylabel('responses centered and normalized to peak, +- ste across flies')
% title('weigthed average across clusters')
% export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/normResponses2octaves_JON_WEDvsAMMC_weightedAverage.eps')
% 
% 
% %% un-weighted average across clusters 
% octDistNull = JON.octDist;
% JON.octavesClustAverages = squeeze(nanmean(octDistNull)); %flyes-octaves
% JON.octavesClustAverages = bsxfun(@rdivide, JON.octavesClustAverages, JON.octavesClustAverages(:,4));
% 
% 
% octDistNull = AMMC.octDist;
% AMMC.octavesClustAverages = squeeze(nanmean(octDistNull)); %flyes-octaves
% AMMC.octavesClustAverages = bsxfun(@rdivide, AMMC.octavesClustAverages, AMMC.octavesClustAverages(:,4));
% 
% octDistNull = WED.octDist;
% WED.octavesClustAverages = squeeze(nanmean(octDistNull)); %flyes-octaves
% WED.octavesClustAverages = bsxfun(@rdivide, WED.octavesClustAverages, WED.octavesClustAverages(:,4));
% 
% 
% errFlyJON = nanstd(JON.octavesWeighted)./sqrt(JON.octavesCountFlies);
% errFlyAMMC = nanstd(AMMC.octavesWeighted)./sqrt(AMMC.octavesCountFlies);
% errFlyWED = nanstd(WED.octavesWeighted)./sqrt(WED.octavesCountFlies);
% 
% figure; hold on
% errorbar(xOctaves, nanmean(JON.octavesClustAverages), errFlyJON, 'Color', cols(1,:))
% errorbar(xOctaves, nanmean(AMMC.octavesClustAverages), errFlyAMMC, 'Color', cols(2,:))
% errorbar(xOctaves, nanmean(WED.octavesClustAverages), errFlyWED, 'Color', cols(3,:))
% 
% legend('JON', 'AMMC', 'WED')
% xlabel('octaves')
% set(gca, 'XTick', xOctaves)
% ylabel('responses centered and normalized to peak, +- ste across flies')
% title('NON-weighted average across clusters')
% export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/normResponses2octaves_JON_WEDvsAMMC_UNweightedAverage.eps')
% 

%% parked stuff
% ttcALL = nanmean(ttc,3);
% ttcALL(CNS.excludeKs,:) = [];
% for Ki = 1:size(ttcALL,1)
%     tc = interp1(xvar, ttcALL(Ki, :), x2);
%     [sigma0(Ki), mu0(Ki)] = gaussfit( 1:2000, tc);
% %     [sigma0(Ki), mu0(Ki)] = gaussfit( xvar, ttcALL(Ki, :)); 
% end
% figure; scatter(mu0, sigma0, 'filled')
% axis image
% xlabel('mu_0 (Hz)')
% ylabel('sigma_0')
% % title('WED: data not interpolated')
% title('WED: data log-interpolated')



%% retrieve AMMC vs WED vs flies relative cluster contribution, take a weighted average (fly by fly)



