% %% spec-JON downsampling
% [~,localhost] = system('hostname');
% if strcmp(localhost(1:end-1), 'Paolas-MacBook-Pro.local') || strcmp(localhost(1:end-1), 'paolas-mbp.private.wireless.med.harvard.edu')
%     ParentFolder = '/Users/galileo/Dropbox (HMS)/Data';
%     DataListFolder = fullfile(ParentFolder,'JON_spec');
% else %must be Archimedes
%     ParentFolder = '/Users/paola/WorkingMemory';
%     DataListFolder = '/Users/paola/WorkingMemory/DataList';
% end
% clear localhost
% 
% load(fullfile(DataListFolder, 'datalist_SPEC_JON.mat'));     %DIR#-required list only
% 
% sets.baseline_threshold = [ 2.2000    2.2000    4.0000    4.0000    8.5000    6.0000    3.9000    3.9000    1.3200    2.2000    1.3200    0.9000    1.3000    1.7000    1.3000    1.3000 ...
%     1.3000    1.7000    1.7000    1.7000    7.0000    4.0000    1.3900    1.2400    1.2400    1.2400    1.2400    0.5000    0.5000    0.6000    0.5000    0.9000 ...
%     1.2500    1.2500    1.5000    2.0000    2.0000    2.0000    1.8400    1.8400    1.8400    1.1500    5.1900    5.0000    1.6700    1.6700    1.6700    1.6700];
% 
% % run sortOut and parametricMaps (unatantum)
% edit sortOut_6F_GREENregistration_specJONdata_spDownsampling.m
% 
% 
% %% panJON signPixels of downsampled data
% ParentFolder = '/Users/galileo/Dropbox (HMS)/Data';
% DataListFolder = fullfile(ParentFolder,'JON_pan');
% load(fullfile(DataListFolder, 'datalist_PAN_JON.mat'));  
% 
% % updtate them with relevant content (steps-only analysis file pointers)
% % old
% for i = 1 : length(datalist_PAN_JON)
%     if isdir(datalist_PAN_JON{i})   % true
%         d = dir([datalist_PAN_JON{i} '/*_analysis_downsampled_stepsOnly.mat']);
%         datalist_PAN_JON{i} = fullfile(datalist_PAN_JON{i}, d.name);
%         sets.analysislist{i} = datalist_PAN_JON{i};
%     else
%         [sets.analysislist{i}, ~] = fileparts(datalist_PAN_JON{i}); % runfolder address
%         d = dir([sets.analysislist{i} '/*_analysis_downsampled_stepsOnly.mat']);
%         datalist_PAN_JON{i} = fullfile(sets.analysislist{i}, d.name);
%     end
% end
% %
% pUSE = 1e-14;     %'ks' 5e-11; %1e-10 used for oldJONs
% nUSE = 1;
% baseline_threshold = 1;
% for z = 12 : length(datalist_PAN_JON)
%     disp(z)
%     fprintf('loading aZ...\n')
%     aZ = load(datalist_PAN_JON{z});
%     fprintf('saving %s...\n', [aZ.basename '_analysis_downsampled.mat'])
%     save(fullfile(sets.analysislist{z}, [aZ.basename '_analysis_downsampled.mat']), '-struct', 'aZ', '-v7.3');
%     disp('computing sign pixels...')
%     aZ = parametricMaps_beta_downsampled(aZ, pUSE, nUSE, 'tt', baseline_threshold);
% end

clear
%% restart all-JON after downsampling is done
ParentFolder = '/Users/galileo/Dropbox (HMS)/Data';
DataListFolder = fullfile(ParentFolder,'JON_PAN&SPEC');
Folder2Save = fullfile(DataListFolder, 'Downsampled');
mkdir(Folder2Save)
cd(Folder2Save)

% get datalist - %provide datalist.mat or otherwise need to DIR##data
load(fullfile(DataListFolder, 'datalist_allJON.mat'));     %DIR#-required list only
NaZ = length(datalist_allJON);
% convert it to downsampled version
for i = 1 : length(datalist_allJON)   % true
    [sets.analysislist{i}, ~] = fileparts(datalist_allJON{i}); % runfolder address
    d = dir([sets.analysislist{i} '/*_analysis_downsampled.mat']);
    datalist_allJON{i} = fullfile(sets.analysislist{i}, d.name);
end

sets.datalist = datalist_allJON; %analysis.mat file
[flyfolders, basenames] = extractRunfoldersBasenames(sets.datalist);
sets.runs_basename = extractRunsBasename_comboFlies( sets, basenames, basenames);
sets.flyfolder = Folder2Save;
clear i d 

%
for z = 1 : NaZ
    aZ{z} = matfile(datalist_allJON{z});
end


%% settings (needed?)
% sets.clustering             = 'GMM';    % 'GMM'; anything else goes for hierarchical
% % sets.matchGMM2Hierarchical  = 1;  % this is by default here.
% sets.loadexistingM          = 1;
% sets.loadexistingclusters   = 1;
% 
% sets.prefix_nK              = 20; % use BIC to automatically determine best # of clusters
% 
% sets.saveM2disk             = 1;
% sets.redefine_pix2use       = 0;  %make sure you specify the combination you want down here and use 1 to do it automatically; otherwise do it manually.
sets.useSignificantPxls     = 2;    % 0: use all the pixels.
                                    % 2: each z gets its own.
                                    % 1: use signif pixels of z==1 only.
                                    % 3: use the set union of all the significant pixels. as in the case of multiple runs at the same z.
                                    % just specify the full set of pixels you want to consider.
% % sets.useSignificantPxls     = ROIs(:,2);
% sets.downsample             = 0;
% sets.export                 = 1;
% sets.chooseClusts           = 0;
% % set all of the following ones in combination (need to have only one 1):
% sets.irespSec               = [0.2,0.2]; %if sets.iresp_center specs are 0, it adds to 0 and stimDur respectively. [0.2, 0.2] standard.
% sets.iresp_center2stOnset   = 1; % always available: both sets.irespSec add to 0.
% sets.iresp_center2stOnsetANDpip = 0; % also centers to pips with
%                                      % standard [0.2,0.2]
% sets.iresp_center2stOffset  = 0; % always available: both sets.irespSec add to stimDur.
% sets.iresp_center2Pips      = 0; % available only if stim is 'PI_DCoffset_PipStimulus'.
% 
% set(0,'DefaultFigureWindowStyle','docked')

%% una tantum - R amtrix
% % make R matrix, and save it
[pixel2keep, iKeep, R] = makeSaveR_mf_zscored( aZ, sets );
% tag = '';
% % % % boxcar out courtship songs.
% % % % guide to choosing stClasses:
% % % 1. 'Chirp_down'
% % % 2. 'Chirp_up'
% % % 3. 'CourtshipSong_A'
% % % 4. 'CourtshipSong_B'
% % % 5. 'PI_DCoffset'
% % % 6. 'PI_DCoffset_PipStimulus'
% % % 7. 'PipStimulus'
% % % 8. 'Tone' (of 'fakeChirp' in older versions)
% % tag = 'noCourts';
% % clear stKeep110
% % T = aZ{1}.T;
% % for stCl = [1, 2, 5:8]
% %     stKeep110(stCl).idx = 1:length(T(1,stCl).trialNums);
% %     stKeep110(stCl).trN = T(1,stCl).trialNums;
% % end
% % boxcar110 = rereferenceBoxcar_82vs110( iKeep(1), stKeep110, {T.trialNums});
% % nameRefSave = fullfile(DataListFolder, sprintf('boxcarring_fly110_%s.mat', tag));
% % save(nameRefSave, 'stKeep110', 'boxcar110', '-v7.3')
% 
% % I will deal with the stimulus later, and save it here as well.

%% import R matrix, and boxcar it, dim reduction
load(fullfile(Folder2Save, 'R_matrix_Downsampled_.mat'));
% 
% addpath('/Users/galileo/GitHub_PP/2P_analysis/drtoolbox/techniques')
% tic
% pca_dims = ceil(intrinsic_dim(R, 'EigValue'));
% toc
% disp(pca_dims) %its 5

rmpath('/Users/galileo/GitHub_PP/2P_analysis/drtoolbox/techniques')
% [coeff,score,latent,~,explained,mu] = pca(R); %using SVD. So far I've
% been using EIG
tic
[coeff,scores,latent,~,explained,mu] = pca(R, 'Algorithm', 'eig', 'NumComponents', 50); %, 'NumComponents',pca_dims);
toc
save('pca50components_eig_output.mat', 'coeff', 'scores', 'latent', 'explained', 'mu')
%% scree plot and explained var
figure; subplot(2,1,1)
plot(latent(1:50))
ylabel('eigenvalues')
subplot(2,1,2)
plot(cumsum(explained(1:50)))
ylabel('cumulative explained variance')

export_fig('PCA_screeplot_explainedVariance.pdf')

%% PC plot
pca_dims = 50;
for ch_i = 1 : 5 : pca_dims
    figure; hold on
    for i = 1:5
        subplot(5, 1, i)
        plot(coeff(:, ch_i + i - 1))
        ax = gca;
        ax.XGrid = 'on';
%         ax.XTick = 0:38:228;
        ax.XLim = [0, 3069];
        ylabel(sprintf('PC %d', ch_i + i - 1))
        
    end
    if ch_i == 1
        export_fig('PCs.pdf')
    else
        export_fig('PCs.pdf', '-append')
    end
end
%%
pca_dims = 19;
coeff = coeff(:,1:pca_dims);
scores = scores(:,1:pca_dims);
explained = sum(explained(1:pca_dims));
 
%% tsne settings
theta           = 0.25;
perplexity      = 50;
tsne_dim        = 3;
reps            = 100000; 
namerTSNESave = fullfile(Folder2Save, sprintf('ready4tSNE_finalPanSpec_downsampled.mat'));
save(namerTSNESave, 'scores', 'tsne_dim', 'pca_dims', 'perplexity', 'theta', 'reps') 
% disp('saved:')
% disp(namerTSNESave)
% 
% 
% %% performing tSNE
% clearvars -except ParentFolder DataListFolder namerTSNESave theta perplexity tsne_dim reps scores pca_dims
% 
% identifier_tSNE = mod(round(datenum(datetime('now'))*1e4), 1e3);
% disp('computing fast tSNE')
% t_tsne = tic;
% M = fast_tsne(scores, tsne_dim, pca_dims, perplexity, theta, [], reps, []); 
% % mappedX = fast_tsne(X, no_dims, initial_dims, perplexity, theta, alg, max_iter, tag)
% toc(t_tsne)
% %
% disp('saving M to disk')
% filenameM = sprintf('tSNE_onRegularPCs_%dD_%d.mat', tsne_dim, identifier_tSNE);
% save(fullfile(DataListFolder, filenameM), 'M', 'perplexity', 'theta', 'pca_dims', 'reps','identifier_tSNE', '-v7.3')
% 

%% do clustering if needed, otherwise skip section....
edit hierarchical_structuredWard_PANSPEC_JONs_ward.m

%% save chunks data for these clusters
edit chunks_response_frequencyTuning_poolOut.m %this also loads all relevant stuff
edit chunks_PLOT_response_frequencyTuning.m

%% if need to load all relevant stuff
Folder2Save = '/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled';
cd(Folder2Save)
clusterFolder = '/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/singleIteration_clustEuclidWard_nodist_noNoisypixels/distance_10000.00';
[datalist_allJON, basenames, dataTable] = loadDatalist_panspec_downsampled;
sets.datalist = datalist_allJON; %analysis.mat file
[flyfolders, ~] = extractRunfoldersBasenames(sets.datalist);
sets.analysislist = flyfolders;
sets.runs_basename = extractRunsBasename_comboFlies( sets, basenames, basenames);
sets.flyfolder = Folder2Save;
load('tSNE_onRegularPCs_3D_191.mat', 'pca_dims', 'identifier_tSNE') %make sure this dataset matches ward_py. There's no control right now.
load('R_matrix_Downsampled_smallWindow.mat', 'R', 'iKeep', 'pixel2keep')
NaZ = length(datalist_allJON);
for z = 1:NaZ
    aZ{z} = matfile(datalist_allJON{z});
end
loadCluster = matfile(fullfile(clusterFolder, 'klustWard_dist10000.00.mat'));
sizeMap = [60, 86];
clust = remapIndexedPointsToRun_RoisMaps(loadCluster.T, pixel2keep, sizeMap);
ROIs = clust.ROIs; %these are still unsorted, so do not remove the last one!
T = loadCluster.T;
load(fullfile(clusterFolder, 'chunks_data_PANSPEC.mat'));
sortFMatrix = [1     2     3     4     5     6    11     8    17     9    12    13    15     7    18    19    16    22    23    24    14   10    20    21];
[~,idx_sortingFlies] = sort(sortFMatrix);
flyfoldersUnique = unique(flyfolders, 'stable');
NaF = length(flyfoldersUnique);
for i=1:NaF, flyNumUnique(i) = str2double(flyfoldersUnique{i}(end-5:end-3)); end    % unsorted by FMatrix, but sortd by flyfoldersUnique!!!
snK = size(ROIs,2);
load(fullfile(clusterFolder,'alignedMatrix_Images.mat'), 'AllClusterMaps', 'flyClusterCount') % flyClusterCount is unsorted by FMatrix, but sortd by flyfoldersUnique!!!

%% plot chunks by k and flies
k = 4;
% find flies that have at least (collectively) 100 pixels in this cluster.
flies = flyClusterCount(:,k) >= 100;
ZS = ismember(dataTable.fly, flyNumUnique(flies));
plot_chunkData('chirps_down_up', loadCluster, clust, chunks, k, find(ZS), dataTable, flyNumUnique(flies));
























%% first GMM pass: calculate or load ROIs - fix this. Seems redundant, but parts may be useful
% t1 = tic;
% 
% importlist = [];
% if sets.loadexistingclusters
%     importlist = dir([sets.flyfolder '/*expROIdata_comboRuns_*' sets.runs_basename '*.mat']);
%     if length(importlist) > 1
%         str = {importlist.name};
%         [selection,validation] = listdlg('PromptString','Select a file:',...
%             'SelectionMode','single','ListSize', [500 100], ...
%             'ListString',str);
%         if validation
%             importlist = importlist(selection);
%         else
%             error('Define the ROIs datafile...')
%         end
%     end
% end
% 
% if ~isempty(importlist)                             % load clusters
%     load(fullfile(sets.flyfolder, importlist.name));
%     nK = size(ROIs,2);
%     if strfind(identifierALL,'.')
%         identifier_tSNE = identifierALL(1:strfind(identifierALL,'.')-1);
%         if isfield(sets, 'identifier_tSNE')
%             if ~strcmp(identifier_tSNE, num2str(sets.identifier_tSNE))
%                 warning('GMM and hierarchical clusters calculated over two different tSNE implementations')
%             end
%         else
%             sets.identifier_tSNE = identifier_tSNE;
%         end
%     end
% else                                                % make clusters
%     % GMM clustering
%     if isempty(sets.prefix_nK)
%         ReplicatesK = 80;
% %         [j, m]   = sortedOut_GMM_selectBIC_JM(        sets.flyfolder, M );
% j = 2;
% m = 2;
%         settings = sortedOut_GMM_selectBIC_K_compute( sets.flyfolder, M, j, m, ReplicatesK );
%         settings.datenow = []; %empty it if you want to use all the available sets in the folder
%         bK       = sortedOut_GMM_selectBIC_K_extract( sets.flyfolder, settings);
%         clusterData = sortedOut_GMM_finalCount( sets, M, j, m, bK );
%     else
%         [j, m]   = sortedOut_GMM_selectBIC_JM(        sets.flyfolder, M, sets.prefix_nK );
%         clusterData = sortedOut_GMM_finalCount( sets, M, j, m, sets.prefix_nK );
%     end
%     identifierALL = clusterData.identifier;
%     sizeMAP = size(aZ{1}.pixels2keep);
%     [ROIs, MAPs, clusterData] = sortedOut_GMM( clusterData, sizeMAP, sets, M, pixel2keep); %using idx now
% %     [ROIs, MAPs, clusterData] = sortedOut_GMM( clusterData, size(aZ(1).pixels2keep), sets, M, pixel2keep); %using idx now
%     
% %     [bootITR, clustSizes] = bootstrapITR(ROIs, aZ, sets.flyfolder, identifierALL, 50, 0, 500);
% %     bootIRR = bootstrapIRR(ROIs, aZ, sets.flyfolder, identifierALL, 50, 0, 500);
% %     clear importlist bK identifier_tSNE j k m
% %     nK = size(ROIs,2);
% %     identifierALL = num2str(identifierALL);
% %     colors = distinguishable_colors(size(MAPs,3));
% %     cMAPs = makeColoredMaps(MAPs, colors);      % convert maps in color
% %     plotClusterMAPs_acrossZs;
% end
% 
% outperm = 1:nK; %no shuffling at this point.
% toc(t1)

%% load existing clustering...

% add from WED_AMMC
% load('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/expROIdata_comboRuns_JONall_7457_12cl_7457.198397.mat', 'MAPs','ROIs', 'identifierALL', 'pixel2keep');
load('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/expROIdata_comboRuns_JONall_7457_15cl_7457.198710.mat', 'MAPs','ROIs', 'identifierALL', 'pixel2keep');
%                 M: [222892x3 double]
%              MAPs: [4-D double]
%              ROIs: [10370x15x79 double]
%       clusterData: [1x1 struct]
%     identifierALL: '7457.165805'
%        pixel2keep: [10370x79 logical]
identifier_tSNE = str2double(identifierALL(1:strfind(identifierALL,'.')-1));
sets.runs_basename = sprintf('JONall_%d',identifier_tSNE);
sets.runs_basename
accessoryfolder = sprintf('clusterMAPS_%s',identifierALL);
accessoryfolder = fullfile(sets.flyfolder, accessoryfolder);

%% plot clustermaps across zs
colors = distinguishable_colors(size(MAPs,3));
cMAPs = makeColoredMaps(MAPs, colors); 

genotype2plot = 0;
tag = sprintf('panJONs_%s', identifierALL);
indices2plot = dataTable.genotype == genotype2plot;
plotClusterMAPs_acrossZs(cMAPs, dataTable, indices2plot, accessoryfolder, tag);


tag = sprintf('JO4_JO32_JO28_JO29_%s', identifierALL);
indices2plot = ( dataTable.genotype == 4 | dataTable.genotype == 32 | dataTable.genotype == 28 | dataTable.genotype == 29 );
plotClusterMAPs_acrossZs(cMAPs, dataTable, indices2plot, accessoryfolder, tag);

tag = sprintf('JO15_JO2_JO22_JO23_JO26_%s', identifierALL);
indices2plot = ( dataTable.genotype == 15 | dataTable.genotype == 2 | dataTable.genotype == 22 | dataTable.genotype == 23 | dataTable.genotype == 26 );
plotClusterMAPs_acrossZs(cMAPs, dataTable, indices2plot, accessoryfolder, tag);


% save figure
%     fname = fullfile(accessoryfolder, sprintf('cluster%02d_MAP_%s.pdf', k, identifierALL));
%     export_fig(fname, '-m10')
%     close


%% una tantum : resave all AVG imgs from registered stack
% 
% for z = 1:length(datalist_allJON)
%     disp(z)
%     FileTif=fullfile(aZ{z}.runfolder,sprintf('registered_%s.tif', basenames{z}));
%     InfoImage=imfinfo(FileTif);
%     mImage=InfoImage(1).Width;
%     nImage=InfoImage(1).Height;
%     NumberImages=length(InfoImage);
%     stack=zeros(nImage,mImage,NumberImages,'uint16');
%     tic
%     TifLink = Tiff(FileTif, 'r');
%     for i=1:NumberImages
%         TifLink.setDirectory(i);
%         stack(:,:,i)=TifLink.read();
%     end
%     TifLink.close();
%     
%     stack = mean(stack, 3); 
%     stack = mat2gray(stack);
%     stack = adapthisteq(stack);
%     stack = imadjust(stack);
% %     figure; imshow(stack)
% 
%     % write tiff
%     writetiff(stack, fullfile(aZ{z}.runfolder, sprintf('AVG_registered_%s.tif', basenames{z})));
% end

%% make and save labeled manual ROIs
t = datetime('now');
t.Format = 'yyMMdd_HHmmss';
accessoryfolder = sprintf('manualROIs_%s',t);
accessoryfolder = fullfile(sets.flyfolder, accessoryfolder);
mkdir(accessoryfolder)

if ~exist('manualROIs', 'var') %declare it
    MIJ.start;
    sizeMap = size(aZ{1}.pixels2keep);
    
    manualROIs = table( repmat({''},100,1), ... % ROI
        repmat({''},100,1), ...       % MAP
        nan(100,1), ...   % z
        repmat({''},100,1), ...       % anatomicalLabel
        'VariableNames', {'poligROI', 'MAP', 'z', 'anatomicalLabel'});
    count_mROI = 0;
    exitManualRoiSelection = 0;
else
    count_mROI = 100 - sum(isnan(manualROIs.z));
    exitManualRoiSelection = 0;
end

% 1. Select which datalist's run to process. 
% 2. Draw ROI.
% 3. Continue the paused program
% 4. Annotate anatomically
% 5. Continue or exit loop
while ~exitManualRoiSelection
    z = input('Choose run from dataTable: ');
    stacked = dataTable.avgF_map{z};
    stacked(:,:,2) = single(aZ{z}.pixels2keep);
    MIJ.createImage(stacked);
    
    disp('after selecting an ROI, hit any key to continue: ')
    pause
        
    count_mROI = count_mROI+1;
    roiIJ = MIJ.getRoi(0);
    MIJ.closeAllWindows();
    pol = roiIJ(:)';
    mapIJ = poly2mask(double(roiIJ(1,:)), double(roiIJ(2,:)), sizeMap(1), sizeMap(2));
    % figure; imagesc(mapIJ); colormap(gray), axis image, axis off
    manualROIs.MAP(count_mROI) = {mapIJ};
    manualROIs.z(count_mROI) = z;
    manualROIs.poligROI(count_mROI) = {pol};
    a = input('anatomical annotation: ');
    manualROIs.anatomicalLabel(count_mROI) = {a};
    exitManualRoiSelection = input('Exit? (1, otherwise 0): ');
end
%% chunk for manual ROIs
if ~exist('chunks','var')
    chunks = chunkStimuli(aZ{1}.fastStimMeta)
end
chunkname = 'steps';

[hfig_chunk, hax_chunk] = plotChunk_mf_manualROIS(chunkname, manualROIs, aZ, sets, ROIs, chunks, dataTable);

figname = fullfile(accessoryfolder, sprintf('%s_%s.pdf', chunkname, t) );
export_fig(figname)
save(fullfile(accessoryfolder, sprintf('manualROIs_%s.mat', t)), 'manualROIs', '-v7.3' )
close

%% chunk stimuli
if ~exist('chunks','var')
    chunks = chunkStimuli(aZ{1}.fastStimMeta)
end
chunkname = 'steps';
clusters2plot = 1:9;
for z = 1:10:length(aZ)
    zplanes2plot = z:min(z+9, length(aZ));   %empty for all
    [hfig_chunk, hax_chunk] = plotChunk_mf(chunkname, clusters2plot, zplanes2plot, aZ, sets, ROIs, chunks);
    figname = fullfile(accessoryfolder, sprintf('%s_%s_fromz%02d.pdf', chunkname, sets.runs_basename, z) );
    export_fig(figname)
    close
end

% hfig_chunk.Position =  [-1    0.05    2    0.85];

%% concatenate stimuli (in their respective relevant window) and calculate psths

z = 1; %assuming now that the stimuli and other conditions are the same across all runs.

n_stimclasses = length(aZ(z).T); %MAJOR TYPES OF STIMULI (e.g PipStimulus, Chirp...)
load(fullfile(aZ(z).runfolder, [aZ(z).basename '_preprocessed.mat']), 'metaimage');
fs_dec = metaimage.acq.frameRate;
fs_acq = 4e4;


% fixing 160919. Goal is to show the whole baseline and poststimwindow too
stimTraces = [];
ts_Stim_all = [];
zeros_st_idxs = [];                                 % changed. Now it tracks the onset of each stimulus
divs_stimuli_idxs = 1;                              % previous zeros_st_idx, onset of a trial chunck (onset of baseline)
divs_classes_idxs = 1;                              % what's this? 
for stCl = 1 : n_stimclasses
    for n = 1 : length(aZ(z).T(stCl).TrTypesThisClass)
        %% stimulus trace 
        sT = aZ(z).T(stCl).sensor(n,:);
        if stCl*n == 1
            DCoffset = mean(sT(20:2e4)); %half a second. Baseline is usually longer than this.
        end
        %         sT = sT(iKeep(z).stClass(stCl).Response_ts); %this was screwing
        %         everything up...
        sT = sT(:);                         % now keeping it all
        stimTraces = cat(1,stimTraces,sT);  % full concatenated vector
        
        
        %% timestamps 1
        %         ts_stim = aZ(z).T(stCl).ts(iKeep(z).stClass(stCl).Response_ts);
        ts_stim = aZ(z).T(stCl).ts;  % keep it all. This starts negative, passes through zero at onset, and it's in seconds.
        if length(ts_stim)~=length(sT)
            croplength = min(length(ts_stim),length(sT));
            ts_stim(croplength+1:end) = [];
            sT(croplength+1:end) = [];
        end
        %% indices 1
        zero_idx = length(ts_Stim_all) + find(ts_stim >= 0, 1);
        zeros_st_idxs = cat(1, zeros_st_idxs, zero_idx);
        
        
        %% timestamps 2
        %make it monotonically increasing (more like absolute timing now)
        ts_stim = ts_stim - ts_stim(1); %starts from zero now
        if isempty(ts_Stim_all)
            ts_Stim_all = ts_stim;
        else
            ts_stim_sum = ts_stim + ts_Stim_all(end) + 1/fs_acq;
            ts_Stim_all = cat(2, ts_Stim_all, ts_stim_sum);
        end
        %% indices 2
        divs_stimuli_idxs = cat(1, divs_stimuli_idxs, length(ts_Stim_all)); %as before, 
                                     % these indices are relative to stimTraces and ts_Stim_all
        
    end
    divs_classes_idxs = cat(1, divs_classes_idxs, n+divs_classes_idxs(end)); %as before
                                     % these indices are relative to divs_stimuli_idxs
end
stimTraces = (stimTraces-DCoffset) * 3; %um conversion %3 for new 30um piezo
zeros_st_ts = ts_Stim_all(zeros_st_idxs);

% % divs_classes_idxs = unique([divs_classes_idxs', 26,42,58,74,90]);     %ad hoc %WEDGE_full
% % divs_classes_idxs = unique([divs_classes_idxs', 9:7:44]);         %ad hoc?
% % divs_classes_idxs = unique([divs_classes_idxs', 11:7:44]);        %ad hoc JONs final
divs_classes_logic = zeros(size(zeros_st_ts));
divs_classes_logic(divs_classes_idxs) = 1;

% divs = divisors(aZ(z).NStimuli);              % what was this>???
% [~,idx] = min(abs(divs - 5.5));
% divs = 1; %divs(idx);




%% plot all



for z = 1 : NaZ
    close all
    %calculate number of repetitions of trials
    N_stimReps = zeros(1,n_stimclasses);
    for stCl = 1 : n_stimclasses   % !!      %  length(aZ)
        N_stimReps(stCl) = max(sum(aZ(z).T(stCl).iTrials));        %SEPARATING stim types
    end
    ncols_traces_DFFs_singK = max(N_stimReps);
    
    % calculate psths
    clear traces_DFFs
    %hybrid version that works for concatenating runs with different Ns.
    ts_dec_all = [];
    for k = 1:size(ROIs,2)
        traces_DFFs_singK = [];
        zeros_stims = 0;
        for stCl = 1 : n_stimclasses   % !!
            % sort stimuli (all, individually)%assuming the same number of
            % repets, otherwise change
            t = squeeze(reshape(aZ(z).T(stCl).tiff, size(aZ(z).T(stCl).tiff,1)*size(aZ(z).T(stCl).tiff,2),1,size(aZ(z).T(stCl).tiff,3), size(aZ(z).T(stCl).tiff,4)));
            t_baseline = mean(squeeze(mean(t( logical(ROIs(:,k,z)),:,aZ(z).T(stCl).indices.ibaseline),1)),2);
            t = squeeze(mean(t( logical(ROIs(:,k,z)),:,:),1));
            t_dFF = bsxfun(@rdivide, (bsxfun(@minus,t,t_baseline)), t_baseline).*100;
            
            traces = [];
            ts_dec = aZ(z).T(stCl).ts_dec - aZ(z).T(stCl).ts(1); %now aligned to the fine ts's zero, as it really is. Does not start exactly from 0.
            deltaEndTss = aZ(z).T(stCl).ts(end) - aZ(z).T(stCl).ts_dec(end);
            for i = 1 : length(aZ(z).T(stCl).TrTypesThisClass)
                tr = t_dFF(aZ(z).T(stCl).iTrials(:,i), :)'; %npointsinResponse x nRepetitions
                traces = cat(1, traces, tr);
                zeros_stims = cat(1, zeros_stims, size(t,2));
                
                if k == 1
                    if isempty(ts_dec_all)
                        ts_dec_all = ts_dec;
                    else
                        ts_dec_sum = ts_dec + ts_dec_all(end) + deltaEndTss + 1/fs_acq; %adding the remaining time at the end plus a small intretrial space, same as I do in the fine ts case
                        ts_dec_all = cat(1, ts_dec_all, ts_dec_sum);
                    end
                end
                
            end
            
            if N_stimReps(stCl) < ncols_traces_DFFs_singK
                addendum = nan(size(traces,1), ncols_traces_DFFs_singK - N_stimReps(stCl));
                traces = cat(2, traces, addendum);
            end
            traces_DFFs_singK = cat(1, traces_DFFs_singK, traces);
        end
        traces_DFFs(:,:,k) = traces_DFFs_singK;
        PS_DFFs_averages = squeeze(nanmean(traces_DFFs,2));
    end
    
    zeros_stims = cumsum(zeros_stims);
    
    clear sT t_dFF t t_baseline tr traces_DFFs_singK i traces_DFFs_sing_z addendum
    
    
    
    %% 0. plot final - setup
    
%     clust_use=[1 5 18 24 38 42 48]; %just the first tot, in order as defined by outperm
    clust_use = 1:nK;
    
    % no stimTable/ALLstim: do here, ignore Error in plotPSTHs_combRuns_showIndicesInUse (line 31)
    
    
    relWidthsColumns = [0.89, 0.11];
    % relWidthsColumns = [0.9, 0.01, 0.09];
    [hfig, hax] = figureTracesI_PP( length(clust_use)+1, relWidthsColumns );
    hfig.Name = 'Selected clusters';
    hfig.Position = [0    0.0044    1, 0.95];
    tbox = annotation('textbox');
    tbox.Position = [0.7667    0.9462    0.1600    0.0500];
    tbox.EdgeColor = [1,1,1];
    tbox.String = sprintf('%s\n%s\nID: %s',aZ(z).basename, sets.runs_basename, identifierALL);
    tbox.Interpreter = 'none';
    tbox.HorizontalAlignment = 'center';
    tbox.FontSize = 13.5;
    
    
    stCl = 1;
    % sets.irespSec               = [0.2, 0.2]; %if sets.iresp_center specs are 0, it adds to 0 and stimDur respectively
    % sets.iresp_center2stOnset   = 0; % always available: both sets.irespSec add to 0.
    % sets.iresp_center2stOffset  = 0; % always available: both sets.irespSec add to stimDur.
    % sets.iresp_center2Pips      = 0; % available only if stim is 'PI_DCoffset_PipStimulus'.
    assert(sets.iresp_center2stOnset+sets.iresp_center2stOffset+sets.iresp_center2Pips <=1, ...
        'Only one single sets.iresp_center2x can be true.')
    if sets.iresp_center2stOnset
        iresp = aZ(z).T(stCl).ts_dec >= sets.irespSec(1) & aZ(z).T(stCl).ts_dec <= sets.irespSec(2);
    elseif sets.iresp_center2stOffset
        iresp = aZ(z).T(stCl).ts_dec >= aZ(z).T(stCl).stimDur + sets.irespSec(1) & aZ(z).T(stCl).ts_dec <= aZ(z).T(stCl).stimDur + sets.irespSec(2);
    else
        iresp = aZ(z).T(stCl).ts_dec >= sets.irespSec(1) & aZ(z).T(stCl).ts_dec <= aZ(z).T(stCl).stimDur + sets.irespSec(2);
        if strcmp(fastStimMeta.className{1}, 'PI_DCoffset_PipStimulus') && sets.iresp_center2Pips
            assert(length(uni.singlePipDur)==1,'Center2Pips with more than one pipDurs is currently not supported.')
            assert(length(uni.pipLatency)==1,'Center2Pips with more than one pipLatencies is currently not supported.')
            iresp = aZ(z).T(stCl).ts_dec >= uni.pipLatency + sets.irespSec(1) & ...
                aZ(z).T(stCl).ts_dec <= uni.pipLatency + uni.singlePipDur + sets.irespSec(2);
        end
    end
    ibasl = aZ(z).T(stCl).ts_dec <=0; %GCamp6S
    
    % plotPSTHs_combRuns_showIndicesInUse( aZ(z), iKeep(z), iresp, ibasl, hax(1,3), DCoffset  );
    
    
    %% 1. concatenated responses
    
    % no stimTable/ALLstim: do this, and then evaluate block:
    % divs = divisors(aZ(z).NStimuli);
    
    % % xticks_st = 0;
    % % for z = 1:n_stimclasses
    % xticks_st = zeros_st; %(:,2);
    % % end % there's some imperfection here....
    % xticks_stMAJ = xticks_st(1:divs:end);
    %
    % xticks_trMAJ = zeros_stims(1:divs:end);
    
    
    for k = 1:length(clust_use)
        ku = clust_use(k);
        axes(hax(k+1,1))
        if outperm(ku) <= size(ROIs,2)
            plot(ts_dec_all,traces_DFFs(:,:,outperm(ku)), 'Color', [0.5 0.5 0.5], 'LineWidth',0.5); hold on
%         plot(ts_dec_all,PS_DFFs_averages(:,clust_use(k)), 'Color', colors(clust_use(k), :), 'LineWidth',2);
            plot(ts_dec_all,PS_DFFs_averages(:,outperm(ku)), 'Color', colors(outperm(ku), :), 'LineWidth',2);
            xlim([ts_dec_all(1), ts_dec_all(end)]);
            %     ylim([yLMin, yLMax]);
            axis tight
            hax(k+1,1).YLabel.String = num2str(outperm(ku));
        end
        hax(k+1,1).YLabel.Units = 'normalized';
%         hax(k+1,1).YLabel.Position = [-0.07 0.02 0];
        hax(k+1,1).YLabel.HorizontalAlignment = 'right';
        hax(k+1,1).YLabel.VerticalAlignment = 'bottom';
        hax(k+1,1).XTick = zeros_st_ts;
        %     hax(k+1,1).GridColor = [0,0,0];
        hax(k+1,1).GridAlpha = 0.3;
    end
    
    
    
    
    % plot spectrogram
    hExtraS = axes('Units', 'normalized','Color','none');
    hExtraS.Position = hax(1,1).Position;
    spectrogram(stimTraces,10000,500,10000,4e4,'yaxis');    % used to be 20000 - 500 - 20000. Trying to distort less the temporal dimension.
    colorbar off
    hExtraS.YLim = [0,0.65];
    hExtraS.YTick = 0:0.2:0.6;
    hExtraS.FontSize = 8;
    hExtraS.XAxisLocation = 'top';
    if isempty(regexpi(hExtraS.XLabel.String, 'secs'))
        hExtraS.XTick = 0:1/6:hExtraS.XTick(end);
        hExtraS.XTickLabel = round(hExtraS.XTick*60);
    end
    hExtraS.XColor = [0.65,0.65,0.65];
    hExtraS.XLabel.String = 'Time (seconds)';
    hExtraS.TickDir = 'out';
    
    
    
    % plot stimulus
    axes(hax(1,1))
    plot(ts_Stim_all,stimTraces, 'k')
    xlim([ts_Stim_all(1), ts_Stim_all(end)]);
    hax(1,1).XTick = zeros_st_ts;
    %     hax(1,1).GridColor = [0,0,0];
    hax(1,1).GridAlpha = 0.3;
    hax(1,1).YAxisLocation = 'right';
    hax(1,1).YLim = [-15, 15];
    ylabel('\mum')
    % hax(1,1).YLabel.Position = 1.0e+04 *[-22.5985    0.0005   -0.0001];  %make the wjole thing normalized in um

    
    hExtraS.Box = 'off';
    
    %% 3. add cluster map and reliability index
    %1: implement inter-trial reliability measure as the mean of all pair-wise
    %correlations between trials
    nTrialsMin = min(N_stimReps); %check, across z
    idxs = ~triu(ones(nTrialsMin));
    for k = 1:size(ROIs,2)
        pwcs = corr(traces_DFFs(:,1:nTrialsMin,k));
%         intertrialRel(k) = mean(mean(pwcs(idxs))); %are you sure you want to take a mean like that? corr coefficients are not additive.
        z_corrs = fisherz(pwcs(idxs));
        intertrialRel(k)= ifisherz(mean(z_corrs));
    end
    % new inter-trial reliability
    
    
    
    
    
    
    for k = 1:length(clust_use)
        ku = clust_use(k);
        axes(hax(k+1,length(relWidthsColumns)))
        hax(k+1,length(relWidthsColumns)).YAxisLocation = 'right';
        hax(k+1,length(relWidthsColumns)).YLabel.HorizontalAlignment = 'left';
        if outperm(ku) <= size(ROIs,2)
            imagesc(MAPs(:,:,outperm(ku),z));
            ylabel(sprintf('%1.2f',intertrialRel(outperm(ku))));
            axis image
        end
        axis off
        hax(k+1,length(relWidthsColumns)).YLabel.Visible = 'on';
        hax(k+1,length(relWidthsColumns)).YLabel.Rotation = 0;
        hax(k+1,length(relWidthsColumns)).YDir = 'reverse';
    end
    k = 1;
    hax(k+1,length(relWidthsColumns)).YLabel.HorizontalAlignment = 'left';
    hax(k+1,length(relWidthsColumns)).YLabel.String = sprintf('inter-trial\nreliability\n%1.2f',intertrialRel(outperm(ku)));
    hax(k+1,length(relWidthsColumns)).YLabel.VerticalAlignment = 'bottom';
    
    delete(hax(1,2))
    
    %% save fig
    figure(hfig)
    name = fullfile(sets.flyfolder, sprintf('analysis_%s_%s',sets.runs_basename, identifierALL));
    % do not overwrite
    d = dir([name,'*']);
    name = sprintf('%s_%d.png', name, length(d) );
    % save
    export_fig(name, '-m3.5')
    % export_fig(name)
    
    %% save fig eps
    
    % name = fullfile(sets.flyfolder, sprintf('analysis_%s_%d',sets.runs_basename, identifier));
    % % do not overwrite
    % d = dir([name,'*']);
    % name = sprintf('%s_%d.eps', name, length(d) );
    % % save
    % export_fig(name)
    %
    % %% save fig eps spectrogram
    %
    % name = fullfile(sets.flyfolder, sprintf('analysis_spectrogram_%s_%d',sets.runs_basename, identifier));
    % % do not overwrite
    % d = dir([name,'*']);
    % name = sprintf('%s_%d.eps', name, length(d) );
    % % save
    % export_fig(name)
    
    %% colored map
    % draw map with colored ROIs
% % old    
%     mapClean = zeros(size(MAPs, 1), size(MAPs,2));
%     for i = 1:nK %index of colors
%         j = outperm(i);     %index of clusters
%         mapClean(MAPs(:,:,j,z)==10) = i; %i.e. when i==1, takes the cluster in first position and gives it color # 1.
%     end
%     Col = [1,1,1];
%     Col = cat(1, Col, colors(1:nK,:)); 


    mapClean = zeros(size(MAPs, 1), size(MAPs,2));
    for k = 1:length(clust_use) %index of cluster and colors within outperm
        ku = clust_use(k);
        if outperm(ku) <= size(ROIs,2)
            mapClean(MAPs(:,:,outperm(ku),z)==10) = ku; %i.e. when i==1, takes the cluster in first position and gives it color # 1.
        end
    end
    Col = [1,1,1];
    Col = cat(1, Col, colors(outperm(1:nK),:)); 
    
    figure('Color', [1,1,1]);
    imagesc(mapClean)
    axis off
    axis image
    colormap(Col(1:max(unique(mapClean)),:))
    title(sprintf('%s\n%s\n%s',aZ(z).basename, sets.runs_basename, identifierALL), 'Interpreter', 'none')
    
    % save fig
    name = fullfile(sets.flyfolder, sprintf('ColoredMap_%s_%s',sets.runs_basename, identifierALL));
    % do not overwrite
    d = dir([name,'*']);
    name = sprintf('%s_%d.png', name, length(d) );
    % save
    export_fig(name, '-m3.5')
    % export_fig(name)
    

end

