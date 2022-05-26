% load up
load_panSpec_downsampled;
load('/Users/galileo/Dropbox (HMS)/Data/chunks_stimfinal.mat', 'chunks')
% load('/Users/galileo/Dropbox (HMS)/Data_Raw_Metadata_BackedUp/BackedUp_Metadata/fly171_run01_metadata_2016-12-20_163728.mat', 'stimuli')
load('/Users/galileo/Dropbox (HMS)/Data/cmap19_spectral.mat', 'cmap10')

iRun = 14;
aZ = matfile(datalist_allJON{iRun}); %174_run01
T = aZ.T;
fastStimMeta = aZ.fastStimMeta;
clear aZ
for z = 1:2% NaZ
    aZ{z} = matfile(datalist_allJON{z});
end

load(fullfile(Folder2Save, 'R_matrix_Downsampled_smallWindow.mat'), 'pixel2keep', 'iKeep', 'pxKeep', 'oldpixel2keep')
    
clusterFolder = '/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix';
treeFileName = 'output_ward_movAVG_ZSCpix_linkage.mat';
cd(clusterFolder)
load('klust_maxClust10.mat')
nK = length(klust);
% load('alignedMatrix_Images_maxClust15.mat')
% cmap10 = Col(2:end,:);
% assert(size(R,1) == sum(pxKeep))


folder2stimuli = '/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/figure2n';
% folder2figure1 = '/Users/galileo/Dropbox (HMS)/figures/figure1';
folder2figure2 = '/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/figure2n';
% folder2figure3 = '/Users/galileo/Dropbox (HMS)/figures/figure3';
cropbaseline = - 1.7; %first point give errors
saveName = fullfile(Folder2Save, 'R_chunks.mat');
ylims = [-50, 400];

clusterOrder = [4 1 2 3 7 5 6 8 9 10]; % this is K2use!! bottom-up
cOrder2 = [1 2 3 4 8 5 7 6 9 10];
[~, cmapOrder] = sort(cOrder2);
cmap10 = cmap10(cmapOrder, :);      %color resorted based on cdf
cmap10 = cmap10([4:end, 3 2 1], :);     %index with Ki
figure; imagesc(1:size(cmap10,1)); colormap(cmap10); title('resorted cmap10')


%% 1. ONLY check peak order for now.
clust = remapIndexedPointsToRun_RoisMaps(Tparent, pixel2keep, sizeMap);
ROIs = clust.ROIs; %these are still unsorted, so do not remove the last one!
trialTypes = {chunks.pips.amplitude_1, chunks.pips.amplitude_2, chunks.pips.amplitude_3, ...
              chunks.tones.amplitude_1_carrPhase_1, chunks.tones.amplitude_2_carrPhase_1 };
[xvar, ia, ic] = unique([chunks.pips.carrierFreqLevels; chunks.tones.carrierFreqLevels(1)]);

load('fTC_data_01_06secAfterOnset.mat')
for i = 1 : length(trialTypes)
    for k = 1 : size(ROIs,2)
        fTC(i).dffPeak(:,:,k) = ( (max(fTC(i).movAvgResponse(:,:,k)) - fTC(i).meanNoise(:,:,k)) ./ fTC(i).meanNoise(:,:,k) ).* 100;
    end
end
for Ki = 1:nK
    % amplitude 3: pips
    yvar3 = [fTC(3).dffPeak(:,:,Ki), NaN(1,8)]; %dff values
    yvar3 = yvar3(ia);
    % add tone 4Hz
    yvar2 = fTC(5).dffPeak(:,:,Ki); %dff values
    % put together
    yvar = [yvar2(1), yvar3(2:end)];
    TCfinal(Ki,:) = yvar;
end
[~, peakFreqIdx] = max(TCfinal, [], 2);
disp('sorted peaks, from push to low freq:')
peakFreqIdx(clusterOrder)


%% try considering cdfs
% first, interpolate tuning curves at linearly-or-log spaced intervals.
onset = 4;
offset = 600; %cropping before, makes things worse
cutoff = 0.5;

x2 = logspace(log10(onset), log10(offset), 2000);
clear cdfs

figure; hold on
plot([onset,offset], [cutoff cutoff], ':', 'Color', [0.6 0.6 0.6])
for Ki = [1:7, 9:10]
    K2use = clusterOrder(Ki)
    yvar2 = interp1(xvar, TCfinal(K2use,:), x2 );
    cdfs(K2use,:) = cumsum(yvar2)./max(cumsum(yvar2));
    plot(x2, cdfs(K2use,:), 'Color', cmap10(Ki,:));
%     xvar3 = 1:0.01:19;
%     yvar3 = interp1(1:19, cdfs(Ki,:), xvar3 );
    [~, halff] = find(cdfs(K2use,:)>=cutoff,1);
    halfpoint(K2use) = x2(halff);
    disp(halfpoint(K2use))
    plot([halfpoint(K2use), halfpoint(K2use)], [0,cutoff], 'Color', cmap10(Ki,:));
    drawnow;
    pause
end
[a,b] = sort(halfpoint);
b
h = gca;
h.XTick = a;
h.XTickLabel = b;

xlim([onset,offset])
ylabel('frequency tuning cdf')
xlabel('frequency log-spaced between 4 and 600 Hz. Ticks represent K2use')
% export_fig(fullfile(folder2figure2, 'cdf_freqTuningCurved.eps'))


h.XAxis.MinorTick = 'off';
h.XTick =  [8, 16, 25, 50, 100, 150, 225, 350, 600];
h.XTickLabel =  [8, 16, 25, 50, 100, 150, 225, 350, 600];
h.XGrid = 'on';
xlabel('');
% export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/cdf_freqTuningCurved_FREQAXISCORRECT.eps')

h.XTick =  [4, 8, 16, 25, 50, 100, 150, 225, 350, 600];
h.XTickLabel =  [4, 8, 16, 25, 50, 100, 150, 225, 350, 600];
h.XGrid = 'on';
xlabel('');

h.XAxis.Scale = 'log';
h.XLim = [4,600];

h.XAxis.MinorTick = 'off';
h.XMinorGrid = 'off';
% export_fig('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/cdf_freqTuningCurved_JON_LOG.eps')


%% make normalized tuning curves distribution (ex fig 6)
TCnorm = bsxfun(@rdivide, TCfinal, max(TCfinal,[],2) );

crop_N_lastFreqs = 0;
[hfig, hax] = figureTracesI_PP_pixels( 1, 400, 100 );
Kgroup = 1;
axes(hax(Kgroup)), hold on
for Ki = [1:7, 9:10]
    K2use = clusterOrder(Ki);
    plot(xvar(1:end-crop_N_lastFreqs), TCnorm(K2use,1:end-crop_N_lastFreqs), 'Color', cmap10(Ki,:))
end
% plot([xvar(1), xvar(end-crop_N_lastFreqs)],[0 0], '-k'); 

hax(Kgroup).XLim = [xvar(1), xvar(end-crop_N_lastFreqs)];
hax(Kgroup).YLim = [min( [0, min(min(TCnorm))] ), 1 ];
hax(Kgroup).XAxis.Visible = 'on';
hax(Kgroup).XTick = xvar(1:end-crop_N_lastFreqs);
hax(Kgroup).YTick = 0:0.5:1;
hax(Kgroup).XAxis.Scale = 'log';
hax(Kgroup).XMinorGrid = 'off';
hax(Kgroup).XMinorTick = 'off';
hax(Kgroup).XGrid = 'off';
% hax(Kgroup).XTickLabel = xvar(1:end-crop_N_lastFreqs);
hax(Kgroup).XTickLabelRotation = 90;
hax(Kgroup).XTickLabel  = [];

export_fig(fullfile(folder2figure2, 'normalizedTuningCurves_JONs.eps'))



%% FIX: RECOLOR ALL MAPS ACTUALLY
for z = 1:length(datalist)
    [runfolders{z}, ~] = fileparts(datalist{z});
end
clust = remapIndexedPointsToRun_RoisMaps(Tparent, pixel2keep, sizeMap);

% make maps (divided old/new because different dimensions)
klust = makeplot_cMaps_superKlusts_singleIteration(clust, klust, dendrOrder);
folder = 'alignedMaps';
mkdir(folder)
cd(folder)

nK = length(klust);

% cmap10 = cmap11([2,3,1,4:11],:); %new cmap 20170914
Col = cat(1,[1,1,1],cmap10);

regGeneralData = matfile('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/anatomyImages/anatomy_alignment_metadata.mat');
Zclust = 1;

%
for f = 1:NaF

    flyfolder = flyfoldersUnique{f};
    ZZ = strfind(datalist, flyfolder);
    zs_fly = zeros(size(ZZ));
    for i = 1:length(zs_fly)
        if ~isempty(ZZ{i})
            zs_fly(i) = 1; % indices into datalist relative to fly f
        end
    end
    zs_fly = logical(zs_fly);

    % load fly's stack2fly150 transformation
    stack2fly150_info = matfile(fullfile(flyfolder, ['stack2fly150sStack_' basenames{find(zs_fly,1)}(1:6) '.mat']));
    if ~isfield(stack2fly150_info, 'Roriginal') && ~strcmp(basenames{find(zs_fly,1)}(1:6), 'fly150')
        stack2fly150_info.Properties.Writable = true;
        stack2fly150_info.Roriginal = imref2d(size(stack2fly150_info.recovered));
    end
    %% sort runs from ventral to dorsal
    zetas = find(zs_fly);
    sliceNums = regGeneralData.alignedSliceNumbers(1,zetas);
    [sliceNums,b] = sort(sliceNums);
    zetas = zetas(b);

    %% load specific runs data and transformations and compute, correcting for downsampling
    for iz = 1:sum(zs_fly)
        z = zetas(iz);
        sliceNum = sliceNums(iz);
        % load specific run's transformation
        regData = matfile(fullfile(runfolders{z}, ['stackRegistration_' basenames{z} '.mat']));
        if ~isfield(regData, 'Roriginal')
            regData.Properties.Writable = true;
            regData.Roriginal = imref2d(size(regData.target));
        end

        al2fly150_mapAll = [];
        for Ki = 1:length(unique(Tparent))
            K2use = clusterOrder(Ki);
            mapClean = klust(K2use).mapsZold(:,:,z);
            mapClean(mapClean == nK+1) = 0;
            mapClean(mapClean>=1) = 1; %binary so far

            %convert from downsampled to fullsize
            mapClean = cat(1, mapClean, zeros(1,size(mapClean, 2))); %61 x 86
            mapClean = cat(2, mapClean, zeros(size(mapClean, 1), 4)); %61 x 90
            mapClean = imresize(mapClean, [86 128]); %should become 86 x 128. The quality sucks, but so be it
            %re-crop
            mapClean(end,:) = [];
            mapClean(:,end-5:end) = [];
            mapClean(mapClean<0.3) = 0;
            mapClean(mapClean>=0.1) = 1;                        %note 1 and not K2use!

            mapClean = imrotate(mapClean, regGeneralData.totAng(1,z)); %preliminary transformation
%             recovered_mapClean = imwarp(mapClean,regData.tform,'OutputView',regData.Roriginal,'FillValues',0); %indexed image
            recovered_mapClean = imwarp(mapClean,regData.tform,'OutputView',regData.Roriginal,'FillValues',0); %indexed image
            
            if strcmp(basenames{find(zs_fly,1)}(1:6), 'fly150')
                al2fly150_mapClean = recovered_mapClean;
            else
                al2fly150_mapClean = imwarp(recovered_mapClean,stack2fly150_info.tform,'OutputView',stack2fly150_info.Roriginal,'FillValues',0);
            end
            al2fly150_mapClean(al2fly150_mapClean<=0.5) = 0;
            al2fly150_mapClean(al2fly150_mapClean>=0.1) = Ki; %dendrogram is colored in order from 1st to nK-esim  in colormap. Clusters are colored accordingy.
            if isempty(al2fly150_mapAll)
                al2fly150_mapAll = al2fly150_mapClean;
            else
                al2fly150_mapAll(al2fly150_mapClean == Ki) = Ki; %dendrogram is colored in order from 1st to nK-esim  in colormap. Clusters are colored accordingy.
            end
        end

        AllClusterMaps(:,:,z) = al2fly150_mapAll;

        AlphaMatrix = zeros(size(al2fly150_mapAll(:,:,1)));
        AlphaMatrix(al2fly150_mapAll(:,:,1) > 0) = 1;

        al2fly150_mapAll = ind2rgb(uint8(al2fly150_mapAll), Col);
        savename = sprintf('alignALLklustsMap_fly%02d_ordZ%02d_%s.png', f, z, basenames{z});
        imwrite(al2fly150_mapAll, savename, 'Alpha', double(AlphaMatrix))
    end
end




%% chunk 1: frequency sweeps
CN = 'chirps_down_up';
% chunks.(CN)
% stimDur = 12;
sortedStimNs = {62, 61}
sortedClassNs = zeros(1, length(sortedStimNs));
for i = 1: length(sortedStimNs)
    sortedClassNs(i) = aZ{1}.stimuli2Class(sortedStimNs{i}(1), 1);
end
% [R_chirps, ~] = makeSaveR_mf_sortedStimuli( aZ, sets, pixel2keep, sortedStimNs, sortedClassNs, cropbaseline );
% save(saveName, 'R_chirps') %will be updated
load(saveName, 'R_chirps', 'timestamps_chirps')

% CHIRPS - plot
[stimTraces, ts_Stim_all, stN, clN, zeros_st_idxs] = makeSortedStimulusSet_protot(T, sortedStimNs, sortedClassNs, cropbaseline);
stimTraces = stimTraces.*3; %um

sortedStimNs = cat(2, sortedStimNs{:});
i = 1;
effectiveFirst  = T(sortedClassNs(i)).ts_dec( find(T(sortedClassNs(i)).ts_dec >= cropbaseline, 1, 'first') );
effectiveLast   = T(sortedClassNs(i)).ts_dec( find(T(sortedClassNs(i)).ts_dec <= T(sortedClassNs(i)).stimDur + 10, 1, 'last') );
iKeepResponse = T(sortedClassNs(i)).ts_dec >= effectiveFirst   & T(sortedClassNs(i)).ts_dec <= effectiveLast;
ts_dec = T(sortedClassNs(i)).ts_dec(iKeepResponse) ;
deltaEndTss = T(sortedClassNs(i)).ts(end) - T(sortedClassNs(i)).ts_dec(end);
deltaStartTss = -(cropbaseline - ts_dec(1));

ts_dec_all = ts_dec - ts_dec(1) + deltaStartTss; %now it is correct!!
zero_tr_idxs = find(ts_dec>=0,1);
stN_tr = sortedStimNs(i) *ones(length(ts_dec), 1);

for i = 2 : length(sortedStimNs) %assuming they ally belong to the same or analogous class
    zero_tr_idxs = cat(1, zero_tr_idxs, find(ts_dec>=0,1) + length(ts_dec_all));
    ts_dec_all = cat(1, ts_dec_all, ts_dec - ts_dec(1) + ts_dec_all(end) + deltaStartTss + deltaEndTss );
    stN_tr = cat(1, stN_tr, sortedStimNs(i)*ones(length(ts_dec), 1) );
end

% start plotting
[hfig, hax] = figureTracesI_PP_pixels( nK+1, 200, 50 );
hax = flipud(hax); % 1 is bottom, as indendrogram
for Ki = 1:nK
    axes(hax(Ki))
    hold off
    K2use = clusterOrder(Ki);
    allR = R_chirps(Tparent==K2use,:);
    avg = mean(allR,1);
    plot([ts_dec_all(1); ts_dec_all; ts_dec_all(end); ts_dec_all(1)], [0, avg, 0, 0], 'Color', cmap10(Ki,:)), hold on, axis tight
    
    hax(Ki).XLim =  [ts_Stim_all(1), ts_Stim_all(end)];
    hax(Ki).XAxis.Visible = 'off';
    hax(Ki).Box = 'off';
    drawnow
end
for Ki = 1:nK
    hax(Ki).YLim = ylims;
    hax(Ki).YTick = [];
    hax(Ki).YAxis.Visible = 'off';
end
YaxReference = floor(diff(ylims)/100)*100;
plot([0,0], [ylims(1), ylims(1) + YaxReference], '-k')
hax(Ki).YAxis.Label.String = YaxReference;
hax(Ki).YAxis.Label.Visible = 'on';
hax(Ki).YAxis.Label.FontSize = 9;
XaxReference = min([floor(hax(Ki).XLim(end)/2)*2, 10]);
plot([0,XaxReference], [ylims(1), ylims(1)], '-k')
hax(Ki).XTick = [];
hax(Ki).XAxis.Label.String = sprintf('%d / %.2f',XaxReference, hax(Ki).XLim(end));
hax(Ki).XAxis.Label.Visible = 'on';
hax(Ki).XAxis.Label.FontSize = 9;
pause

% export and save!
export_fig(fullfile(folder2figure2, 'resp_chirps.eps'))

%% chunk 2: adaptation to long tones
CN = 'tones';
% chunks.(CN)
sortedStimNs = {[85 89 93 75 77 79 81 83], [86 90 94 76 78 80 82 84]};
% sortedStimNs = {76};
sortedClassNs = zeros(1, length(sortedStimNs));
for i = 1: length(sortedStimNs)
    sortedClassNs(i) = aZ{1}.stimuli2Class(sortedStimNs{i}(1), 1);
end
% tic
% [R_tones, ~] = makeSaveR_mf_sortedStimuli( aZ, sets, pixel2keep, sortedStimNs, sortedClassNs, cropbaseline );
% toc
% save(saveName, 'R_tones', '-append')
load(saveName, 'R_tones');
[stimTraces, ts_Stim_all, stN, clN, zeros_st_idxs] = makeSortedStimulusSet_protot(T, sortedStimNs, sortedClassNs, cropbaseline);
stimTraces = stimTraces.*3; %um
sortedStimNs = cat(2, sortedStimNs{:});

% TONE(s) - plot
% timestamps!
i = 1;
effectiveFirst  = T(sortedClassNs(i)).ts_dec( find(T(sortedClassNs(i)).ts_dec >= cropbaseline, 1, 'first') );
effectiveLast   = T(sortedClassNs(i)).ts_dec( find(T(sortedClassNs(i)).ts_dec <= T(sortedClassNs(i)).stimDur + 10, 1, 'last') );
iKeepResponse = T(sortedClassNs(i)).ts_dec >= effectiveFirst   & T(sortedClassNs(i)).ts_dec <= effectiveLast;
ts_dec = T(sortedClassNs(i)).ts_dec(iKeepResponse) ;
deltaEndTss = T(sortedClassNs(i)).ts(end) - T(sortedClassNs(i)).ts_dec(end);
deltaStartTss = -(cropbaseline - ts_dec(1));

ts_dec_all = ts_dec - ts_dec(1) + deltaStartTss; %now it is correct!!
zero_tr_idxs = find(ts_dec>=0,1);
stN_tr = sortedStimNs(i) *ones(length(ts_dec), 1);

for i = 2 : length(sortedStimNs) %assuming they all belong to the same or analogous class
    zero_tr_idxs = cat(1, zero_tr_idxs, find(ts_dec>=0,1) + length(ts_dec_all));
    ts_dec_all = cat(1, ts_dec_all, ts_dec - ts_dec(1) + ts_dec_all(end) + deltaStartTss + deltaEndTss );
    stN_tr = cat(1, stN_tr, sortedStimNs(i)*ones(length(ts_dec), 1) );
end

%so far all the stimuli from the previos section are included (so both R
%and stimTraces are in the same framework. However, only selectStimuli_stN
%will be plotted down here.
% start plotting


% sel_stN = 76;
sel_stN = [86 90 94 80];
i_St = ismember(stN, sel_stN);
i_Tr = ismember(stN_tr, sel_stN);

ts_Stim_all_sel = ts_Stim_all(1:sum(i_St));
% onset_idx_sel = find(ismember( ts_Stim_all_sel, ts_Stim_all(zeros_st_idxs)));
ts_Stim_all_sel = ts_Stim_all_sel - ts_Stim_all_sel(1);
ts_dec_all_sel = ts_dec_all(1:sum(i_Tr));
ts_dec_all_sel = ts_dec_all_sel-ts_dec_all_sel(1);


for Ki = 1:nK
    axes(hax(Ki))
    hold off
    K2use = clusterOrder(Ki);
    allR = R_tones(Tparent==K2use,:);
    avg = mean(allR,1);
    plot([ts_dec_all_sel(1); ts_dec_all_sel; ts_dec_all_sel(end); ts_dec_all_sel(1)], [0, avg(i_Tr), 0, 0], 'Color', cmap10(Ki,:)), hold on, axis tight

    hax(Ki).XLim =  [ts_Stim_all_sel(1), ts_Stim_all_sel(end)];
    hax(Ki).XAxis.Visible = 'off';
    hax(Ki).Box = 'off';
    drawnow
end
for Ki = 1:nK
    hax(Ki).YLim = ylims;
    hax(Ki).YTick = [];
    hax(Ki).YAxis.Visible = 'off';
end
YaxReference = floor(diff(ylims)/100)*100;
plot([0,0], [ylims(1), ylims(1) + YaxReference], '-k')
hax(Ki).YAxis.Label.String = YaxReference;
hax(Ki).YAxis.Label.Visible = 'on';
hax(Ki).YAxis.Label.FontSize = 9;
XaxReference = min([floor(hax(Ki).XLim(end)/2)*2, 10]);
plot([0,XaxReference], [ylims(1), ylims(1)], '-k')
hax(Ki).XTick = [];
hax(Ki).XAxis.Label.String = sprintf('%d / %.2f',XaxReference, hax(Ki).XLim(end));
hax(Ki).XAxis.Label.Visible = 'on';
hax(Ki).XAxis.Label.FontSize = 9;
pause

% export and save!
export_fig(fullfile(folder2figure2, sprintf('resp_tones_%s_st76_125HzHighAmpl.eps',num2str(sel_stN))))


%% chunk 3: songs
CN = 'courtships';
% chunks.(CN)
% stimDur = 12;
sortedStimNs = {74, 73}
sortedClassNs = zeros(1, length(sortedStimNs));
for i = 1: length(sortedStimNs)
    sortedClassNs(i) = aZ{1}.stimuli2Class(sortedStimNs{i}(1), 1);
end
% [R_songs, ~] = makeSaveR_mf_sortedStimuli( aZ, sets, pixel2keep, sortedStimNs, sortedClassNs, cropbaseline );
% save(saveName, 'R_songs', '-append') %will be updated
load(saveName, 'R_songs')
R_songs_remain = R_songs;

[stimTraces, ts_Stim_all, stN, clN, zeros_st_idxs] = makeSortedStimulusSet_protot(T, sortedStimNs, sortedClassNs, cropbaseline);
stimTraces = stimTraces.*3; %um
sortedStimNs = cat(2, sortedStimNs{:});

% SONGS - plot
% timestamps!

for i = 1:2
    sel_stN = sortedStimNs(i);
    stimDur = T(sortedClassNs(i)).stimDur - 2.6;
    effectiveFirst  = T(sortedClassNs(i)).ts_dec( find(T(sortedClassNs(i)).ts_dec >= cropbaseline, 1, 'first') );
    effectiveLast   = T(sortedClassNs(i)).ts_dec( find(T(sortedClassNs(i)).ts_dec <= T(sortedClassNs(i)).stimDur + 10, 1, 'last') );
    iKeepResponse = T(sortedClassNs(i)).ts_dec >= effectiveFirst   & T(sortedClassNs(i)).ts_dec <= effectiveLast;
    ts_dec = T(sortedClassNs(i)).ts_dec(iKeepResponse) ;
    deltaEndTss = T(sortedClassNs(i)).ts(end) - T(sortedClassNs(i)).ts_dec(end);
    deltaStartTss = -(cropbaseline - ts_dec(1));
    
    ts_dec_all_sel = ts_dec - ts_dec(1) + deltaStartTss; %now it is correct!!
    zero_tr_idxs = find(ts_dec>=0,1);
    stN_tr = sortedStimNs(i) *ones(length(ts_dec), 1);
    i_St = ismember(stN, sel_stN);
    ts_Stim_all_sel = ts_Stim_all(i_St);
    % onset_idx_sel = find(ismember( ts_Stim_all_sel, ts_Stim_all(zeros_st_idxs)));
    ts_Stim_all_sel = ts_Stim_all_sel - ts_Stim_all_sel(1);
    R_songs = R_songs_remain(:,1:length(ts_dec));
    R_songs_remain = R_songs_remain(:,length(ts_dec)+1:end);
    
    [hfig, hax] = figureTracesI_PP_pixels( nK+1, 200, 50 );
    hax = flipud(hax); % 1 is bottom, as indendrogram
    
    for Ki = 1:nK
        axes(hax(Ki))
        hold off
        K2use = clusterOrder(Ki);
        allR = R_songs(Tparent==K2use,:);
        avg = mean(allR,1);
        plot([ts_dec_all_sel(1); ts_dec_all_sel; ts_dec_all_sel(end); ts_dec_all_sel(1)], [0, avg, 0, 0], 'Color', cmap10(Ki,:)), hold on, axis tight       
        hax(Ki).XLim =  [ts_Stim_all_sel(1), ts_Stim_all_sel(end)];
        hax(Ki).XAxis.Visible = 'off';
        hax(Ki).Box = 'off';
        drawnow
    end
    
    for Ki = 1:nK
        hax(Ki).YLim = ylims;
        hax(Ki).YTick = [];
        hax(Ki).YAxis.Visible = 'off';
    end
    YaxReference = floor(diff(ylims)/100)*100;
    plot([0,0], [ylims(1), ylims(1) + YaxReference], '-k')
    hax(Ki).YAxis.Label.String = YaxReference;
    hax(Ki).YAxis.Label.Visible = 'on';
    hax(Ki).YAxis.Label.FontSize = 9;
    XaxReference = min([floor(hax(Ki).XLim(end)/2)*2, 10]);
    plot([0,XaxReference], [ylims(1), ylims(1)], '-k')
    hax(Ki).XTick = [];
    hax(Ki).XAxis.Label.String = sprintf('%d / %.2f',XaxReference, hax(Ki).XLim(end));
    hax(Ki).XAxis.Label.Visible = 'on';
    hax(Ki).XAxis.Label.FontSize = 9;
    pause 
    
    % export and save!
    export_fig(fullfile(folder2figure2, sprintf('resp_songs_%s.eps',num2str(sel_stN))))
end


%% chunk 4: steps
CN = 'steps';
chunks.(CN)
% stimDur = 12;
sortedStimNs = {55:60}
sortedClassNs = zeros(1, length(sortedStimNs));
for i = 1: length(sortedStimNs)
    sortedClassNs(i) = aZ{1}.stimuli2Class(sortedStimNs{i}(1), 1);
end
% [R_steps, ~] = makeSaveR_mf_sortedStimuli( aZ, sets, pixel2keep, sortedStimNs, sortedClassNs, cropbaseline );
% save(saveName, 'R_steps', '-append') %will be updated
load(saveName, 'R_steps')

% STEPS - plot
% timestamps!
[stimTraces, ts_Stim_all, stN, clN, zeros_st_idxs] = makeSortedStimulusSet_protot(T, sortedStimNs, sortedClassNs, cropbaseline);
stimTraces = stimTraces.*3; %um

sortedStimNs = cat(2, sortedStimNs{:});
i = 1;
effectiveFirst  = T(sortedClassNs(i)).ts_dec( find(T(sortedClassNs(i)).ts_dec >= cropbaseline, 1, 'first') );
effectiveLast   = T(sortedClassNs(i)).ts_dec( find(T(sortedClassNs(i)).ts_dec <= T(sortedClassNs(i)).stimDur + 10, 1, 'last') );
iKeepResponse = T(sortedClassNs(i)).ts_dec >= effectiveFirst   & T(sortedClassNs(i)).ts_dec <= effectiveLast;
ts_dec = T(sortedClassNs(i)).ts_dec(iKeepResponse) ;
deltaEndTss = T(sortedClassNs(i)).ts(end) - T(sortedClassNs(i)).ts_dec(end);
deltaStartTss = -(cropbaseline - ts_dec(1));

ts_dec_all = ts_dec - ts_dec(1) + deltaStartTss; %now it is correct!!
zero_tr_idxs = find(ts_dec>=0,1);
stN_tr = sortedStimNs(i) *ones(length(ts_dec), 1);

for i = 2 : length(sortedStimNs) %assuming they ally belong to the same or analogous class
    zero_tr_idxs = cat(1, zero_tr_idxs, find(ts_dec>=0,1) + length(ts_dec_all));
    ts_dec_all = cat(1, ts_dec_all, ts_dec - ts_dec(1) + ts_dec_all(end) + deltaStartTss + deltaEndTss );
    stN_tr = cat(1, stN_tr, sortedStimNs(i)*ones(length(ts_dec), 1) );
end


% start plotting
[hfig, hax] = figureTracesI_PP_pixels( nK+1, 200, 50 );
hax = flipud(hax); % 1 is bottom, as indendrogram

for Ki = 1:nK
    axes(hax(Ki))
    hold off
    K2use = clusterOrder(Ki);
    allR = R_steps(Tparent==K2use,:);
    avg = mean(allR,1);
    plot([ts_dec_all(1); ts_dec_all; ts_dec_all(end); ts_dec_all(1)], [0, avg, 0, 0], 'Color', cmap10(Ki,:)), hold on, axis tight
    
    hax(Ki).XLim =  [ts_Stim_all(1), ts_Stim_all(end)];
    hax(Ki).XAxis.Visible = 'off';
    hax(Ki).Box = 'off';
    drawnow
end

for Ki = 1:nK
    hax(Ki).YLim = ylims;
    hax(Ki).YTick = [];
    hax(Ki).YAxis.Visible = 'off';
end
YaxReference = floor(diff(ylims)/100)*100;
plot([0,0], [ylims(1), ylims(1) + YaxReference], '-k')
hax(Ki).YAxis.Label.String = YaxReference;
hax(Ki).YAxis.Label.Visible = 'on';
hax(Ki).YAxis.Label.FontSize = 9;
XaxReference = min([floor(hax(Ki).XLim(end)/2)*2, 10]);
plot([0,XaxReference], [ylims(1), ylims(1)], '-k')
hax(Ki).XTick = [];
hax(Ki).XAxis.Label.String = sprintf('%d / %.2f',XaxReference, hax(Ki).XLim(end));
hax(Ki).XAxis.Label.Visible = 'on';
hax(Ki).XAxis.Label.FontSize = 9;
pause 

% export and save!
export_fig(fullfile(folder2figure2, 'resp_steps.eps'))




%% remake the entire calculation, not based on single runs, but on the whole cluster.
% pipeline:
% for each cluster
% 1. store ext baseline points -> average -> calculate std (and store it)
% 2. store the fucking response points.
% 3. use a fucking mov average to find 5-points peak -> average (and store) 
% 
% 
sizeMap = [60, 86];
clust = remapIndexedPointsToRun_RoisMaps(Tparent, pixel2keep, sizeMap);
ROIs = clust.ROIs; %these are still unsorted, so do not remove the last one!

trialTypes = {chunks.pips.amplitude_1, chunks.pips.amplitude_2, chunks.pips.amplitude_3, ...
              chunks.tones.amplitude_1_carrPhase_1, chunks.tones.amplitude_2_carrPhase_1 }; %keep separate subsets
saveNameFTData = fullfile(clusterFolder, 'FTdata.mat');

% freqTuning_pipsTones_extendedBaselines_extractTraces(trialTypes, aZ, saveNameFTData); % only run if not already.
movingwindowsize = 1;   %for moving average
responseWindow = 1;     %seconds after stim onset. Empty [] for entire default window
if ~isempty(responseWindow)
    responseWindowPoints = sum(ts_dec>=0 & ts_dec <=responseWindow);
else
    responseWindowPoints = [];
end
% fTC = freqTuning_pipsTones_extendedBaselines(trialTypes, ROIs, movingwindowsize, saveNameFTData, responseWindowPoints); % baseline+ext baseline: 102 points
% load('fTC_data_1secResponse.mat') 
% load('fTC_data_02_07secAfterOnset.mat')
load('fTC_data_01_06secAfterOnset.mat')

% calcluate perc dff of peak response
for i = 1 : length(trialTypes)
    for k = 1 : size(ROIs,2)
        fTC(i).dffPeak(:,:,k) = ( (max(fTC(i).movAvgResponse(:,:,k)) - fTC(i).meanNoise(:,:,k)) ./ fTC(i).meanNoise(:,:,k) ).* 100;
    end
end


%% find max dff values for each cluster from figure 2A
maxAllStimuli = [];
negMaxAllStimuli = [];
m = matfile(saveName);
for i = 1 : 3
    switch i
        case 1
            data = m.R_chirps;
        case 2
            data = m.R_songs;
        case 3
            data = m.R_steps;
    end
    for k = 1 : size(ROIs,2)
        maxAllStimuli(k,i) = max(mean(data(klust(k).k,:)));
        negMaxAllStimuli(k,i) = min(mean(data(klust(k).k,:)));
    end
end
% sanity check\ %plot step responses
% [hfig, hax] = figureTracesI_PP_pixels( 11, 200, 30 );
% hax = flipud(hax); % 1 is bottom, as indendrogram
% for Ki = 1:11
%     axes(hax(Ki))
%     plot(mean(data(klust(Ki).k,:)))
% end
%OK
% incorporate peaks from pips and tones
for i = 1:5
    for k = 1 : size(ROIs,2)
        maxAllStimuli(k,i+3) = max(fTC(i).dffPeak(:,:,k));
        negMaxAllStimuli(k,i+3) = min(fTC(i).dffPeak(:,:,k));
    end
end

dffLimits = [min(negMaxAllStimuli,[], 2), max(maxAllStimuli, [],2)];

%% rescale by baseline noise
% recalculate a single baseline noise level for each cluster
for k = 1 : size(ROIs,2)
    basTemp = [];
    for i = 1:5
        basTemp = cat(1, basTemp, reshape(fTC(i).baselines(:,:,k), [],1,1));
    end
    fullBaselPipsTones(:, k) = basTemp;
end
noiseClusters = std(fullBaselPipsTones); %ok.
scaledNoiseClusters = ( (noiseClusters-min(noiseClusters)) ./ max((noiseClusters-min(noiseClusters))) ) / 2 ; %scaled between 0 (min noise) and 0.5 (max noise). [To be subtracted]
% scaledNoiseClusters =
%     0.0434    0.0339    0.0133    0.5000    0.1400    0.1466    0.0201    0.0395    0.2126    0.0610         0




%% single curve with 4Hz from tones and rest from pips, max amplitude only
[xvar, ia, ic] = unique([chunks.pips.carrierFreqLevels; chunks.tones.carrierFreqLevels(1)]);

[hfig, hax] = figureTracesI_PP_pixels( 11, 200, 50 );
hax = flipud(hax); % 1 is bottom, as indendrogram   %check
for Ki = 1:nK   %fix
    axes(hax(Ki)), hold on
    K2use = clusterOrder(Ki);
    % amplitude 3: pips
    yvar3 = [fTC(3).dffPeak(:,:,K2use), NaN(1,8)]; %dff values
    yvar3 = yvar3(ia);
    % add tone 4Hz
    yvar2 = fTC(5).dffPeak(:,:,K2use); %dff values
    % put together
    yvar = [yvar2(1), yvar3(2:end)];
    TCfinal(K2use,:) = yvar;
    plot(xvar, yvar, 'Color', cmap10(Ki,:))
%     plot(xvar(1), yvar(1), 'x', 'Color', cmap11(Ki,:), 'MarkerSize',3)
    plot([xvar(1), xvar(end)],[0 0], '-k');   
    
    hax(Ki).XLim = [xvar(1), xvar(end)];
    hax(Ki).YLim = dffLimits(K2use, :).*0.65;
    hax(Ki).Box = 'off';
    hax(Ki).XAxis.Visible = 'off';
    hax(Ki).YTick = [];
    hax(Ki).XAxis.Scale = 'log';
    hax(Ki).XMinorGrid = 'off';
    set(gca, 'XTick', xvar)
    hax(Ki).XGrid = 'off'; %was on
    hax(Ki).TickDir = 'out';
    hax(Ki).XLim = [0,xvar(end)];
    if Ki == 1
        hax(Ki).XTickLabel = xvar;
        hax(Ki).XTickLabelRotation = 90;
        hax(Ki).XAxis.FontSize = 4;
        hax(Ki).XAxis.Visible = 'on';
    else
        hax(Ki).XTickLabel = [];
        hax(Ki).XAxis.Visible = 'off';
    end
    hax(Ki).XMinorTick = 'off';
    drawnow
end
pause
export_fig(fullfile(folder2figure2, 'tuningCurves_01_06sec_4HzAndPips_MaxAmplit_65%.eps'))

%% find the maximum at highest amplitude
[~, peakFreqIdx] = max(TCfinal, [], 2);
[hfig, hax] = figureTracesI_PP_pixels( 11, 100, 50 );
hax = flipud(hax); 

for Ki = 1:nK 
    K2use = clusterOrder(Ki);
    if peakFreqIdx(K2use) == 1 
        peakFreqIdx(K2use) = 2; % the closest
    end
    % amplitude 1: pips
    yvar1 = [fTC(1).dffPeak(:,:,K2use), NaN(1,8)]; %dff values
    yvar1 = yvar1(ia);
    yvar1 = yvar1(peakFreqIdx(K2use));
    
    % amplitude 2: pips
    yvar2 = [fTC(2).dffPeak(:,:,K2use), NaN(1,8)]; %dff values
    yvar2 = yvar2(ia);
    yvar2 = yvar2(peakFreqIdx(K2use));
    
    % amplitude 3: pips
    yvar3 = [fTC(3).dffPeak(:,:,K2use), NaN(1,8)]; %dff values
    yvar3 = yvar3(ia);
    yvar3 = yvar3(peakFreqIdx(K2use));
    
    amplCurves(K2use,1) = yvar1;
    amplCurves(K2use,2) = yvar2;
    amplCurves(K2use,3) = yvar3;
    
    
    % plot
    axes(hax(Ki)), hold on
    plot(1:3, [yvar1, yvar2, yvar3], 'Color', cmap10(Ki,:))
    
    plot([0.8, 3.2],[0 0], '-k');
    hax(Ki).XLim = [0,4];
    hax(Ki).YLim = dffLimits(K2use, :).*0.65;
    hax(Ki).Box = 'off';
    hax(Ki).XAxis.Visible = 'off';
    hax(Ki).YTick = [];
%     hax(Ki).XAxis.Scale = 'lin';
    hax(Ki).XMinorGrid = 'off';
    set(gca, 'XTick', 1:3)
    hax(Ki).XGrid = 'off';
    hax(Ki).TickDir = 'out';
    hax(Ki).YAxis.Visible = 'on';
    if Ki == 1
        hax(Ki).XTickLabel = {'low', 'mid', 'high'};
        hax(Ki).XTickLabelRotation = 90;
        hax(Ki).XAxis.FontSize = 5;
        hax(Ki).XAxis.Visible = 'on';
    else
        hax(Ki).XTickLabel = [];
        hax(Ki).XAxis.Visible = 'off';
    end
    hax(Ki).XMinorTick = 'off';
    drawnow 
end

save('/Users/galileo/Dropbox (HMS)/Paola_Manuscript/figures/pixelwise_freqTuning/amplitudeBestCDFStim_JON.mat', 'amplCurves')
pause
export_fig(fullfile(folder2figure2, 'tuningAmplitude_01_06sec_65%.eps'))

