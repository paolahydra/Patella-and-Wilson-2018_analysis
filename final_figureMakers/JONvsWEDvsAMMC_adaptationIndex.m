
%% load everything useful

load('/Users/galileo/Dropbox (HMS)/Data/chunks_stimfinal.mat')
load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/iCNS.mat', 'iCNS', 'cmap') 
load('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix/iJON.mat', 'iJON', 'cmap10') %cmap is sorted already

CNS.clusterFolder = '/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix';
JON.clusterFolder = '/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix';

% save('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix/iJON.mat', 'cmap10', '-append')
% save('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/iCNS.mat', 'cmap', '-append')

%the following fileds have now been loaded in iCNS:

% CNS.datalist = load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/datalist_WEDAMMC_piezo.mat');
% CNS.clusters = load(fullfile(CNS.clusterFolder, 'klust_maxClust19.mat'));
% CNS.excludeKs = [1,2]; % pull, push
% CNS.includedK_i = ~ismember(CNS.clusters.Tparent, CNS.excludeKs);
iCNS.includedK_i = ~ismember(iCNS.Tparent, iCNS.excludeKs);
iCNS.includeKs = 3:19;


% JON.clusters = load(fullfile(JON.clusterFolder, 'klust_maxClust10.mat'));
% JON.excludeKs = [8,9,10]; % brown , pull, push
% JON.includedK_i = ~ismember(JON.clusters.Tparent, JON.excludeKs);
iJON.includedK_i = ~ismember(iJON.Tparent, iJON.excludeKs);
iJON.includeKs = 1:7;

%% chunk's R    -   'chirps_down_up'
% cropbaseline = - 1.7; %first point give errors
load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/Downsampled/R_chunks.mat',  'ts_dec_chirp')

CNS.R_chirps = load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/Downsampled/R_chunks.mat', 'R_chirps');
JON.R_chirps = load('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/R_chunks.mat', 'R_chirps');

CNS.R_chirps = CNS.R_chirps.R_chirps;
JON.R_chirps = JON.R_chirps.R_chirps;

%% split responses by stim and cluster, and find best lag

% make average responses
allChirpResponsesTypes = [];
% JON
for ik = 1 : length(iJON.includeKs) %included clusters
    k = iJON.includeKs(ik);
    allChirpResponsesTypes= cat(1, allChirpResponsesTypes , mean(JON.R_chirps(iJON.Tparent==k, :),1) );
end

% CNS
for ik = 1 : length(iCNS.includeKs) %included clusters
    k = iCNS.includeKs(ik);
    allChirpResponsesTypes= cat(1, allChirpResponsesTypes , mean(CNS.R_chirps(iCNS.Tparent==k, :),1) );
end


% split
allChirpResponsesTypes_down = allChirpResponsesTypes(:,1:167);
allChirpResponsesTypes_up_flipped = fliplr(allChirpResponsesTypes(:,168:end));

%% crosscorrelation, for each type
for k = 1 : size(allChirpResponsesTypes,1)
    [xcf(k,:),lags(k,:),bounds(k,:)] = crosscorr(allChirpResponsesTypes_down(k,:), allChirpResponsesTypes_up_flipped(k,:), 81);
end

% bounds not so informative (all the same)
% lags not useful either
[maxXCF,I] = max(xcf, [], 2);

% yes, DO need to only consider linear ones to find optimal lag:
keepLags = I(maxXCF>0.95)
mean(keepLags)
median(keepLags)
mode(keepLags)

% it's really 63.5
% I will use 64
LAG_i = 64;
% corr_lagged = xcf(:,LAG_i); % not using anymore

%% sort JONS and apply to traces once nd fo all
clusterOrder = fliplr([4 1 2 3 7 5 6]); % this is K2use!! bottom-up
% % double-check --> OK
% figure
% for k = clusterOrder
%     plot(mean(JON.R_chirps(iJON.Tparent==k, :),1) )
%     pause
% end
clusterOrder = cat(2,fliplr([4 1 2 3 7 5 6]), 8:24); % JON and CNS as in figure

allChirpResponsesTypes_up_flipped_sorted = allChirpResponsesTypes_up_flipped(clusterOrder, :);
allChirpResponsesTypes_down_sorted = allChirpResponsesTypes_down(clusterOrder,:);
clear allChirpResponsesTypes_up_flipped allChirpResponsesTypes_down allChirpResponsesTypes

%% need to apply the LAG shift now
% try one direction and check: let's add nans rather than cropping for now
LAG = (lags(1,LAG_i));
LAG = abs(LAG);
allChirpResponsesTypes_up_flipped_sorted_lagged = cat(2, nan(size(allChirpResponsesTypes_up_flipped_sorted,1), LAG), allChirpResponsesTypes_up_flipped_sorted);
allChirpResponsesTypes_down_sorted_lagged = cat(2, allChirpResponsesTypes_down_sorted, nan(size(allChirpResponsesTypes_up_flipped_sorted,1), LAG));

% [hfig, hax] = figureTracesI_PP_pixels( 24, 160, 200);
% for k = 1:24
%     axes(hax(k,1)); hold on
%     plot(allChirpResponsesTypes_up_flipped_sorted_lagged(k,:))
%     plot(allChirpResponsesTypes_down_sorted_lagged(k,:))
%     ylabel(sprintf('%d', k))
% end

% ok! Now clip nans and tails

allChirpResponsesTypes_down_sorted_lagged(:,[1:LAG, end-LAG+1:end]) = [];
allChirpResponsesTypes_up_flipped_sorted_lagged(:,[1:LAG, end-LAG+1:end]) = [];
clear allChirpResponsesTypes_up_flipped_sorted  allChirpResponsesTypes_down_sorted


%% normalize (jointly both up and down) and calculate distance
% just divide by max....
maxvalues = max(cat(2, allChirpResponsesTypes_down_sorted_lagged, allChirpResponsesTypes_up_flipped_sorted_lagged), [], 2);
%bah
zscored = zscore(cat(2, allChirpResponsesTypes_down_sorted_lagged, allChirpResponsesTypes_up_flipped_sorted_lagged), [], 2);
% [hfig, hax] = figureTracesI_PP_pixels( 24, 160, 200);
% for k = 1:24
%     axes(hax(k,1)); hold on
%     plot(zscored(k,150:end))
%     plot(zscored(k,1:149))
%     ylabel(sprintf('%d', k))
% end  %OK

zscored_up = zscored(:,150:end);
zscored_down = zscored(:,1:149);

% calculate distance:
adapt_index = sum((zscored_up - zscored_down).^2, 2);
figure; plot(1:7, adapt_index(1:7), 'x-'); hold on
plot(8:24, adapt_index(8:24), 'x-');

export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/adaptationIndex_sortedClusters.eps')
save('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/adaptationIndex.mat', 'adapt_index')


%% try colored bar plots raterh than lines
load('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/adaptationIndex.mat', 'adapt_index');
runOnce = 0;
if ~runOnce
    adapt_index = adapt_index./max(adapt_index);
    cmap10 = flipud(cmap10(1:7,:));
    figure; imagesc(1:size(cmap10,1)); colormap(cmap10); title('resorted cmap10')
    runOnce = 1;
end
 %plot bars, however, colors are going to be added manually later
 
 
 figure;
bar(adapt_index)
axis off
export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/adaptationIndex_sortedClusters_BARS.eps')


%% double-check sorting again --> OK
% tracesallSorted = allChirpResponsesTypes(clusterOrder,:);
% tracesJONsorted = tracesallSorted(1:7,:);
% tracesCNSsorted = tracesallSorted(8:end,:);
% corr_laggedSorted = corr_lagged(clusterOrder);
% corr_laggedJONSorted = corr_laggedSorted(1:7);
% corr_laggedCNSSorted = corr_laggedSorted(8:end);
% 
% [hfig, hax] = figureTracesI_PP_pixels( 7, 160, 200);
% for k = 1:7
%     axes(hax(k,1)); hold on
%     plot(tracesJONsorted(k,:))
%     ylabel(sprintf('%d: %0.2f', k, corr_laggedJONSorted(k)))
% end
% 
% [hfig, hax] = figureTracesI_PP_pixels( 17, 160, 200);
% for k = 1:17
%     axes(hax(k,1)); hold on
%     plot(tracesCNSsorted(k,:))
%     ylabel(sprintf('%d: %0.2f', k, corr_laggedCNSSorted(k)))
% end

%% figure 6D
load('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/adaptationIndex.mat', 'adapt_index');
adapt_index = adapt_index./max(adapt_index);

% recover amplitude curves
JONamplCurves = load('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/amplitudeBestCDFStim_JON.mat', 'amplCurves'); % unsorted
CNSamplCurves = load('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/amplitudeBestCDFStim_CNS.mat', 'amplCurves');

JON.amplCurves = JONamplCurves.amplCurves;
CNS.amplCurves = CNSamplCurves.amplCurves;
clear JONamplCurves CNSamplCurves

% sort JON curves as in adapt index, and remove clusters
JON.amplCurves = JON.amplCurves(iJON.includeKs,:);
JON.amplCurves = JON.amplCurves(fliplr([4 1 2 3 7 5 6]),:);

CNS.amplCurves = CNS.amplCurves(iCNS.includeKs,:);

allAmplitCurves = cat(1, JON.amplCurves, CNS.amplCurves);

%% hyperb model - fixed exponent
% % x = [-flipud(chunks.pips.amplitudeLevels); chunks.pips.amplitudeLevels]';
% 
% xdata = 15.*(chunks.pips.amplitudeLevels)'; %now in um
% esponente = 1/5;
% fun = @(x,xdata) x(1) .* (xdata.^esponente) ./ (x(2) + (xdata.^esponente) );
% 
% options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
% lb = [];
% ub = [];
% 
% for i = 1:24
%     ydata = allAmplitCurves(i,:);
%     figure
%     plot(xdata, ydata, 'pg')
%     hold on
%     x0 = [(ydata(3)+10), (ydata(3)+10)/2];
%     params2(i,:) = lsqcurvefit(fun,x0,xdata,ydata,lb,ub,options);
%     x = 0:0.005:0.4;
%     plot(x, fun(params2(i,:),x))
%     title(i)
% end

%% hyperb model - constrained parameters! yay... much better

% x = [-flipud(chunks.pips.amplitudeLevels); chunks.pips.amplitudeLevels]';

xdata = 15.*(chunks.pips.amplitudeLevels)'; %now in um
xplot = 0:0.02:3.6;

fun = @(x,xdata) x(1) .* (xdata.^x(2)) ./ (x(3) + (xdata.^x(2)) );

options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
lb = [min(allAmplitCurves(:)), 0.25, 0.05];
ub = [400, 20, 15];

for i = 1:24
    ydata = allAmplitCurves(i,:);
    figure(i)
    plot(xdata, ydata, 'pg')
    hold on
    
    x0 = [(ydata(3)+10), 1.5, (ydata(3)+10)/2];
    params3(i,:) = lsqcurvefit(fun,x0,xdata,ydata,lb,ub,options);
    
    plot(xplot, fun(params3(i,:),xplot))
    title(i)
end


%% now regress against adapt index
halfMaxAmplitudes = params3(:,3) .^ (1 ./ params3(:,2) );
for i = 1:24
    figure(i)
    plot([0,4], [params3(i,1)/2, params3(i,1)/2], '-y')
    plot(halfMaxAmplitudes(i), fun(params3(i,:),halfMaxAmplitudes(i)), 'xr')
end


%%
mkdir('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/fittingCurves')
cd('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/fittingCurves')
for i = 1:24
    figure(i)
    export_fig(sprintf('fc_%02d.eps',i))
end


%% regress
outliers = [13, 18]; % these are the ones for which one parameter takes the limit of the constraining range

keep_ampl = halfMaxAmplitudes(setdiff(1:24, outliers));
keep_adapt_index = adapt_index(setdiff(1:24, outliers));

figure; 
scatter(keep_ampl, keep_adapt_index, 'filled', 'MarkerFaceColor', 'b')
hold on
scatter(halfMaxAmplitudes(outliers), adapt_index(outliers), 'filled', 'MarkerFaceColor', 'r')
legend('data', 'outliers')
xlabel('amplitudes that yields half-max response')
ylabel('adaptation index')
% export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/adaptationIndex_vs_amplitudeHM_fitSigm_outliers.eps')

% regress
figure; 
scatter(keep_ampl, keep_adapt_index, 'filled', 'MarkerFaceColor', 'b')
hold on
x = keep_ampl./max(keep_ampl);
y = keep_adapt_index;
X = [ones(length(x),1) x];
b = X\y
yCalc = X*b;

Rsq = 1 - sum((y - yCalc).^2)/sum((y - mean(y)).^2)

plot(x,yCalc)

xlabel('amplitudes that yields half-max response')
ylabel('adaptation index')
% export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/adaptationIndex_vs_amplitudeHM_fitSigm.eps')


%% recover and save halfMax amplitude vs frequency
% % recover and save  center frequncy for each cluster
% % - CNS:
% [xvar, ia, ic] = unique([chunks.pips.carrierFreqLevels; chunks.tones.carrierFreqLevels(1)]);
% load('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/downsampled/rawDataLinkage/ward_ZSCpix/fTC_data_01_06secAfterOnset.mat')
% for i = 1 : length(fTC)
%     for k = 1 : 19
%         fTC(i).dffPeak(:,:,k) = ( (max(fTC(i).movAvgResponse(:,:,k)) - fTC(i).meanNoise(:,:,k)) ./ fTC(i).meanNoise(:,:,k) ).* 100;
%     end
% end
% %get maximum
% for Ki = 1:19 
%     % amplitude 3: pips
%     yvar3 = [fTC(3).dffPeak(:,:,Ki), NaN(1,8)]; %dff values
%     yvar3 = yvar3(ia);
%     % add tone 4Hz
%     yvar2 = fTC(5).dffPeak(:,:,Ki); %dff values
%     % put together
%     yvar = [yvar2(1), yvar3(2:end)];
%     TCfinal(Ki,:) = yvar;
% end
% 
% 
% onset = 4;
% offset = 600; %cropping before, makes things worse
% cutoff = 0.5;
% x2 = logspace(log10(onset), log10(offset), 2000);
% clear cdfs
% 
% figure; hold on
% plot([onset,offset], [cutoff cutoff], ':', 'Color', [0.6 0.6 0.6])
% for Ki = 1:19
%     yvar2 = interp1(xvar, TCfinal(Ki,:), x2 );
%     cdfs(Ki,:) = cumsum(yvar2)./max(cumsum(yvar2));
%     plot(x2, cdfs(Ki,:), 'Color', cmap(Ki,:));
% %     xvar3 = 1:0.01:19;
% %     yvar3 = interp1(1:19, cdfs(Ki,:), xvar3 );
%     [~, halff] = find(cdfs(Ki,:)>=cutoff,1);
%     halfpoint(Ki) = x2(halff);
%     plot([halfpoint(Ki), halfpoint(Ki)], [0,cutoff], 'Color', cmap(Ki,:));
%     drawnow;
% end
% [a,b] = sort(halfpoint);
% 
% % amplitHalfpoint = halfMaxAmplitudes;
% freqHalfpoint = halfpoint(iCNS.includeKs);
% %% - JON:
% clusterOrder = [4 1 2 3 7 5 6 8 9 10];
% cd(JON.clusterFolder)
% load('fTC_data_01_06secAfterOnset.mat')
% for i = 1 : 5
%     for k = 1 : 10
%         fTC(i).dffPeak(:,:,k) = ( (max(fTC(i).movAvgResponse(:,:,k)) - fTC(i).meanNoise(:,:,k)) ./ fTC(i).meanNoise(:,:,k) ).* 100;
%     end
% end
% for Ki = 1:10
%     % amplitude 3: pips
%     yvar3 = [fTC(3).dffPeak(:,:,Ki), NaN(1,8)]; %dff values
%     yvar3 = yvar3(ia);
%     % add tone 4Hz
%     yvar2 = fTC(5).dffPeak(:,:,Ki); %dff values
%     % put together
%     yvar = [yvar2(1), yvar3(2:end)];
%     TCfinal(Ki,:) = yvar;
% end
% 
% clear cdfs halfpoint
% 
% figure; hold on
% plot([onset,offset], [cutoff cutoff], ':', 'Color', [0.6 0.6 0.6])
% for Ki = iJON.includeKs
%     K2use = clusterOrder(Ki);
%     yvar2 = interp1(xvar, TCfinal(K2use,:), x2 );
%     cdfs(K2use,:) = cumsum(yvar2)./max(cumsum(yvar2));
%     plot(x2, cdfs(K2use,:), 'Color', cmap10(Ki,:));
% %     xvar3 = 1:0.01:19;
% %     yvar3 = interp1(1:19, cdfs(Ki,:), xvar3 );
%     [~, halff] = find(cdfs(K2use,:)>=cutoff,1);
%     halfpoint(Ki) = x2(halff); %NOTE ORDER is FINAL
%     plot([halfpoint(Ki), halfpoint(Ki)], [0,cutoff], 'Color', cmap10(Ki,:));
%     drawnow;
% end
% 
% % put together and save
% freqHalfpoint = cat(2, fliplr(halfpoint), freqHalfpoint); %JONs resorted low to high
% amplitHalfpoint = halfMaxAmplitudes;
% save('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/freqAmplit_24clusters_JON_CNS_halfpoints.mat', 'freqHalfpoint','amplitHalfpoint')
% 

%% plot halfMax amplitude vs frequency


% run block at line 182


load('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/freqAmplit_24clusters_JON_CNS_halfpoints.mat', 'freqHalfpoint','amplitHalfpoint')
outliers = [13, 18]; 
runOnce = 0;
if ~runOnce
    adapt_index = adapt_index./max(adapt_index);
    cmap10 = flipud(cmap10(1:7,:));
%     figure; imagesc(1:size(cmap10,1)); colormap(cmap10); title('resorted cmap10')    
    colorMap = cat(1, cmap10, cmap(3:end,:));
    runOnce = 1;
end

clustNames = {'iv', 'v', 'vi', 'vii' 'viii', 'ix', 'x','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19'};

figure; hold on
for i = setdiff(1:24, outliers)
    scatter(freqHalfpoint(i), amplitHalfpoint(i), 200, 'filled', 'MarkerFaceColor', colorMap(i,:))
    str = sprintf('   %s',clustNames{i});
    text(freqHalfpoint(i), amplitHalfpoint(i), str)
end

% export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/freqVSamplitude_halfMaxs.eps')

ax = gca;
ax.TickDir='out';


ax.YAxis.Scale = 'log';
ax.YAxis.Limits = [0.15, 4.2];
ax.YAxis.TickValues = [0.15 0.225, 0.450, 0.9, 1.800, 3.6];


ax.XAxis.Scale = 'log';
ax.XAxis.Limits = [8, 250];
ax.XAxis.MinorTickValues = [9   10   20    30    40    50    60    70    80    90   100   200];
ax.XAxis.TickValues = [8,16,25,50,100,225];


export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/freqVSamplitude_halfMaxs_logXYAxis.eps')






