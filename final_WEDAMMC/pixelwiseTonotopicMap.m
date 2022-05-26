% plot pixel-wise tonotopy


%% load basic stuff - DO ONCE

%first block of clustering_on_rawData_WEDAMMC.m

cd('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/Downsampled/rawDataLinkage/ward_ZSCpix')
folder2tonotopicMaps = '/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/FigSupplem/pixelwiseTonotopies_toBeRotatedWithGoodAntialiasing';

pip3 = load('FTdata3.mat');     %pips, high amplitude
tones2 = load('FTdata5.mat');   %tones, high amplitude


%% colormap
N = 14;
xvar = xvar(1:N);
extrapoints = 2 - mod(N,2);
cmap = brewermap(N + extrapoints,'RdYlGn');
cmap( (N+extrapoints)/2 : (N+extrapoints)/2 -1+extrapoints, : ) = [];
assert(size(cmap,1) == N);
figure; imagesc(1:size(cmap,1)); colormap(cmap)
ax = gca;
ax.YTick = [];
ax.XTick = 1:N;
ax.XTickLabel = xvar;
ax.Box = 'off';
export_fig(fullfile(folder2tonotopicMaps, 'PixelW_tonotopicMaps_colorbar.tif'))
cmap = cat(1, [1 1 1], cmap);



%% SELECT RUN AND PLOT
cmap = flipud(bone(N+1));
flyNum = 123;
runNum = 2;

z = find(dataTable.fly == flyNum & dataTable.run == runNum) ;

pixel2keep_z = pixel2keep(:,z);

pp.baselines = squeeze(pip3.baselines(:,z,:,:));
pp.responses = squeeze(pip3.responses(:,z,:,:));
tn.baselines = squeeze(tones2.baselines(:,z,:,:));
tn.responses = squeeze(tones2.responses(:,z,:,:));

% get peak tuning for each sig pixel here
peakPrefFreq = pixbypix_freqtuning(pixel2keep_z, pp, tn);
[xvar, ~, ~] = unique([chunks.pips.carrierFreqLevels; chunks.tones.carrierFreqLevels(1)]);



% saturate frequencies bigger than 275Hz
tonotopicMap = zeros(size(pixel2keep_z));
tonotopicMap(pixel2keep_z) = peakPrefFreq.idx;
tonotopicMap = uint8(tonotopicMap);
tonotopicMap = reshape(tonotopicMap, sizeMap);

tonotopicMap(tonotopicMap>N) = N;

figure; imshow(tonotopicMap, cmap);
% figure; image(tonotopicMap); colormap(cmap); %same
%
export_fig(fullfile(folder2tonotopicMaps, sprintf('PixelW_tonotopicMaps_c_%s_bone.tif', basenames{z})))
