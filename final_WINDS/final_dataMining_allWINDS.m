rootname = 'datalist_allWINDS';
rootfolder = 'WINDS_100';

ParentFolder = '/Users/galileo/Dropbox (HMS)/Data';
DataListFolder = fullfile(ParentFolder,rootfolder);
Folder2Save = fullfile(DataListFolder, 'Downsampled');
if exist(Folder2Save, 'dir') ~= 7, mkdir(Folder2Save); end
cd(Folder2Save)

load(fullfile(DataListFolder, [rootname '.mat']));     %DIR#-required list only
datalist = eval(rootname);
clear(rootname)
NaZ = length(datalist);
NaZ_100 = 23;%re-referenced
NaZ_86  = 44; 
% ReferenceFolder = '/Users/galileo/Dropbox (HMS)/Data/fly170_PP/fly170_run01';
% ReferenceData = 'fly170_run01_analysis.mat';

%% run sortOut and parametricMaps (unatantum)
edit sortOut_6F_GREENregistration_resampleWINDS_spDownsampling.m 

%%
for z = 1 : NaZ
    [sets.analysislist{z}, ~] = fileparts(datalist{z}); % runfolder address
    d = dir([sets.analysislist{z} '/*_analysis_downsampled.mat']);
    datalist{z} = fullfile(sets.analysislist{z}, d.name);
end
sets.datalist = datalist;

for z = 1 : NaZ
    aZ{z} = matfile(datalist{z});
end
[flyfolders, basenames] = extractRunfoldersBasenames(sets.datalist);
sets.runs_basename = extractRunsBasename_comboFlies( sets, basenames, basenames);
sets.flyfolder = Folder2Save;
sets.runs_basename = 'allWINDS';
sets.useSignificantPxls     = 2;    % 0: use all the pixels.
                                    % 2: each z gets its own.
set(0,'DefaultFigureWindowStyle','docked')





%% make datalists if not made already... and then retrieve them
% load('/Users/galileo/Dropbox (HMS)/Data/WINDS_100/datalist_WINDS_100.mat');
% load('/Users/galileo/Dropbox (HMS)/Data/WINDS_86/datalist_WINDS_86.mat');
% 
% % update and save datalists for once (HMS)
% pattern_old = '/Users/galileo/Dropbox/';
% pattern_new = '/Users/galileo/Dropbox (HMS)/';
% for i = 1 : length(datalist_WINDS_100)
%     startIdx = strfind(datalist_WINDS_100{i}, pattern_old); %this is usually 1, unless already resaved.
%     if startIdx == 1
%         datalist_WINDS_100{i} = [pattern_new, datalist_WINDS_100{i}(length(pattern_old)+1 : end)];
%     end
% end
% save('/Users/galileo/Dropbox (HMS)/Data/WINDS_100/datalist_WINDS_100.mat', 'datalist_WINDS_100');
%   
% for i = 1 : length(datalist_WINDS_86)
%     startIdx = strfind(datalist_WINDS_86{i}, pattern_old); %this is usually 1, unless already resaved.
%     if startIdx == 1
%         datalist_WINDS_86{i} = [pattern_new, datalist_WINDS_86{i}(length(pattern_old)+1 : end)];
%     end
% end
% save('/Users/galileo/Dropbox (HMS)/Data/WINDS_86/datalist_WINDS_86.mat', 'datalist_WINDS_86');
%%
% datalist2beReref = load('/Users/galileo/Dropbox (HMS)/Data/WINDS_100/datalist_WINDS_100.mat');
% datalistRef = load('/Users/galileo/Dropbox (HMS)/Data/WINDS_86/datalist_WINDS_86.mat');
% datalist2beReref = datalist2beReref.(cell2mat(fieldnames(datalist2beReref)));
% datalistRef = datalistRef.(cell2mat(fieldnames(datalistRef)));
% % updtate them with relevant content (steps-only analysis file pointers)
% % 
% for i = 1 : length(datalist2beReref)
%     if isdir(datalist2beReref{i})   % ok
%         d = dir([datalist2beReref{i} '/*analysis.mat']);
%         datalist2beReref{i} = fullfile(datalist2beReref{i}, d.name);
%     end
% end
% 
% % new: need to swap with steps-only file
% for i = 1 : length(datalistRef)
%     if ~isdir(datalistRef{i})   % ok
%         [folderpath, ~] = fileparts(datalistRef{i});
%         d = dir([folderpath '/*analysis.mat']);
%         datalistRef{i} = fullfile(folderpath, d.name);
%     end
% end
% 
% datalist = cat(2, datalist2beReref, datalistRef); %100 and then 86
% NaZ = length(datalist);
% 
% datalist_allWINDS = datalist;
% save(fullfile(Folder2Save, 'datalist_allWINDS.mat'), 'datalist_allWINDS')
%% save tsRef (unatantum)
% aT = load(fullfile(ReferenceFolder, ReferenceData), 'T');
% tsRef = aT.T(1).ts_dec;
% mkdir(fullfile(Folder2Save, 'NewResamplingData'))
% save(fullfile(Folder2Save, 'NewResamplingData', 'tsRef_fly170.mat'), 'tsRef', '-v7.3');

%% run sortOut and parametricMaps (unatantum)
% edit sortOut_6F_GREENregistration_resampleWINDS_spDownsampling.m 

%% patch my mistake, ma avrei fatto meglio a usare matfile
% for z = 1 : NaZ
%     [sets.analysislist{z}, ~] = fileparts(datalist{z}); % runfolder address (legacy)
%     d1 = dir(fullfile(sets.analysislist{z}, '/*analysis_EXTBSL.mat'));
%     generalAnalysis = load(fullfile(sets.analysislist{z}, d1.name));
%     d = dir(fullfile(sets.analysislist{z}, '/*analysis.mat'));
%     parMapsAnalysis = load(fullfile(sets.analysislist{z}, d.name));
%     
%     generalAnalysis.PixbyPix    = parMapsAnalysis.PixbyPix;
%     generalAnalysis.pixels2keep = parMapsAnalysis.pixels2keep;
%     generalAnalysis.runfolder   = parMapsAnalysis.runfolder;
%     
%     save(fullfile(sets.analysislist{z}, d.name),'-struct', 'generalAnalysis');
%     delete(fullfile(sets.analysislist{z}, d1.name));
%     
% end

%% just redo parametric maps
% pUSE = 0.00005; %86?
% pUSEks = 0.001;
% nUSE = 1;
% baseline_threshold = 1;
% for z = 1 :length(datalist2beReref)
%     aZ = load(datalist2beReref{z});
%     disp(z)
%     disp(aZ.basename)
%     aZ = parametricMaps_beta(aZ, pUSE, nUSE, 'tt', baseline_threshold, pUSEks);
%     close all
% end
% clear aZ
% clear pUSE pUSEks baseline_threshold nUSE

%% make and save R matrix: need to make them separately and then concatenate everything (unatantum)
% 
% [pixel2keep_100, iKeep_100, R_100] = makeSaveR_mf( aZ(1:NaZ_100), sets );
% [pixel2keep_86, iKeep_86, R_86] = makeSaveR_mf( aZ(NaZ_100+1:end), sets );
% 
% % % if need to load
% % R86load = load(fullfile(sets.flyfolder, 'R_matrix_Downsampled_100.mat'));
% % R_86            = R86load.R;
% % iKeep_86        = R86load.iKeep;
% % pixel2keep_86   = R86load.pixel2keep;
% % clear R86load
% 
% 
% % build up
% R = cat(1, R_100, R_86);
% clear R_100 R_86
% 
% z = 1;
% R_100.sizeMAP = size(aZ{z}.pixels2keep);
% R_100.iKeep = iKeep_100;
% R_100.basenames = basenames(1:length(datalist2beReref));
% R_100.pixel2keep = pixel2keep_100;
% R_100.R = R(1:sum(R_100.pixel2keep(:)),:);
% 
% z = NaZ_100+1;
% R_86.sizeMAP = size(aZ{z}.pixels2keep);
% R_86.iKeep = iKeep_86;
% R_86.basenames = basenames(NaZ_100+1:end);
% R_86.pixel2keep = pixel2keep_86;
% R_86.R = R(sum(R_100.pixel2keep(:))+1:end,:);
% 
% % save
% filenameR = fullfile(sets.flyfolder, 'R_matrix_downsampled_winds_100&86.mat');
% save(filenameR, 'R_100', 'R_86', 'R', '-v7.3')
% 
% clear iKeep_86 iKeep_100 pixel2keep_86 pixel2keep_100

% Save zs2 Rmatrix for linkage calculation
load(fullfile(Folder2Save, 'R_matrix_downsampled_winds_100&86.mat')); 

%only consider sustained response phase:
R = reshape(R, size(R,1), [], 4);
R = R(:,12:44,:);
R = reshape(R, size(R,1), [], 1);
Rzs = zscore(R, 0, 2);
% want to clean it up a bit? %maybe not
save('Rzs_winds_downsampled_noprepr.mat', 'Rzs');


%% load R matrix (R, R_100, R_86) - NO NEED TO BOXCAR IT, but should boxcar/produce/save stimtraces
% load(fullfile(DataListFolder, 'R_matrix_downsampled_winds_100&86.mat'));  %'R', 'pixel2keep', 'iKeep'       %not so easy
% pxls100 = sum(R_100.pixel2keep(:));
% pxls86 = sum(R_86.pixel2keep(:));
% sizeMAP = size(aZ{1}.pixels2keep);
% 
% 
% identifier_tSNE = 123;
% sets.identifier_tSNE = identifier_tSNE;
% sets.runs_basename = sprintf('allWINDSl_%d',identifier_tSNE);
% 
% identifierALL = '007';
% 
% 
% load(fullfile(DataListFolder, 'winds8_ManualMapping.mat'))
% % redundBrush = sum(brushedIdx)>1; %not present here
% idx = zeros(1, length(brushedIdx));
% for i = 1 : size(brushedIdx, 1)
%     idx(brushedIdx(i,:)) = i;
% end
% nK = max(idx);

%% PCA step 1
% tic
% pca_dims = ceil(intrinsic_dim(R, 'EigValue'));
% toc
% disp(pca_dims)
% 
% % pca_dims = 2; 
% 
% %% PCA step 2
% pca_dims = 3;
% tic
% [scores, PCs, lambda] = initialPCA_fast_tsne(R, pca_dims); %PCs and lambda are now both sorted in descending order
% toc % <1 sec
% 
% % % plotting is only coarse for now
% % figure;
% % for i = 1 : pca_dims
% %     subplot(pca_dims,1,i)
% %     plot(PCs(:,i));
% % end
% % % clear R PCs lambda i
% 
% %%
% figure; plot(scores(:,1), scores(:,2), '.', 'MarkerSize', 10)
% 
% 
% %% %% tsne settings
% theta           = 0.25;
% perplexity      = 50;
% tsne_dim        = 2;
% reps            = 50000; 
% 
% tag = '';
% namerTSNESave = fullfile(Folder2Save, sprintf('ready4tSNE_allWINDS%s.mat', tag));
% save(namerTSNESave, 'scores', 'tsne_dim', 'pca_dims', 'perplexity', 'theta', 'reps') 
% disp('saved:')
% disp(namerTSNESave)
% 
% 
% %% performing tSNE
% clearvars -except ParentFolder DataListFolder namerTSNESave theta perplexity tsne_dim reps scores pca_dims
% perplexity      = 50;
% clc
% identifier_tSNE = mod(round(datenum(datetime('now'))*1e4), 1e3);
% diary on
% disp('computing fast tSNE')
% disp(identifier_tSNE)
% t_tsne = tic;
% M = fast_tsne(scores(:,1:2), tsne_dim, pca_dims, perplexity, theta, [], reps, []); 
% % mappedX = fast_tsne(X, no_dims, initial_dims, perplexity, theta, alg, max_iter, tag)
% toc(t_tsne)
% %
% disp('saving M to disk')
% filenameM = sprintf('tSNE_onRegularPCs_%dD_%d.mat', tsne_dim, identifier_tSNE);
% save(fullfile(Folder2Save, filenameM), 'M', 'perplexity', 'theta', 'pca_dims', 'reps','identifier_tSNE', '-v7.3')
% diary off


%% converged from final_plotPSTHs_f_winds_noCAT.m, as from 
% edit batchCluster_differentKs_pan_spec.m
% ROIs and MAPs are your main outputs. Same info, two different formats.

%pan
clust(1).ROIs = zeros(size(R_100.pixel2keep, 1), nK, size(R_100.pixel2keep,2));           %2D matrix: rows is your tiff's (x*y), columns is # of z planes considered
clust(1).MAPs = zeros([R_100.sizeMAP, nK, size(R_100.pixel2keep,2)]);     %3D matrix: rows is your tiff's x, columns is your tiff's y, third dimension is # of z planes considered
%spec
clust(2).ROIs = zeros(size(R_86.pixel2keep, 1), nK, size(R_86.pixel2keep,2));           %2D matrix: rows is your tiff's (x*y), columns is # of z planes considered
clust(2).MAPs = zeros([R_86.sizeMAP, nK, size(R_86.pixel2keep,2)]);     %3D matrix: rows is your tiff's x, columns is your tiff's y, third dimension is # of z planes considered


for k = 1:nK        %go through clusters
    idxTot = idx;       % array of clustering indices for all pixels in M
%     for z = 1:size(R_100.pixel2keep,2)     %loop through z planes
%         MAP = zeros(R_100.sizeMAP);   % just your tiff's x by y
%         MAP(R_100.pixel2keep(:,z)) = 1;               % x by y matrix
%         pxp = R_100.pixel2keep(:,z);          % logical (x*y) array of significant pixels for this z plane
%         
%         idxChunk = idxTot(1:sum(pxp(:)));   % sum(pxp(:)) is the number of significant pixels for this z plane. Let's call this number SP.  Take those pixels out.
%         idxTot(1:sum(pxp(:))) = [];     % idxTot stores the remaining clustering indices for the significant pixels of the remaining z planes.
%         Lia = idxChunk == k;            % Lia si a logical, SP-dimensional array. Within the chunk of this z plane, take only the ones that belong to this specific cluster.
%         numLia = find(Lia);  %numbers here still refer to SP and not to the tiff-dimensions.
%         
%         % next you go from a logical (x*y) array [pxp] to an (x*y) array
%         % [pxpCS],
%         % where significant pixels are numbered in ascending order and
%         % everything else is a zero. These numbers will correspond to
%         % numbers in numLia (which is in SP-format).
%         pxpCS = cumsum(pxp(:));     
%         pxpCSDif = diff(pxpCS);
%         on = pxpCSDif == 1;
%         on = [0;on];
%         pxpCS(~on) = 0; %this has the same code as numLia now
% 
%         %find intersection with specific cluster, and store it in ROIs and
%         %MAPs:
%         finalLia = ismember(pxpCS, numLia);  
%         clust(1).ROIs(:,k,z) = finalLia;
%         if sum(finalLia) ~= 0
%             MAP(finalLia)=10;
%         else
%             MAP(1)=10;
%         end
%         clust(1).MAPs(:,:,k,z) = MAP;
%     end
    
    for z = 1:size(R_86.pixel2keep,2)     %loop through z planes
        MAP = zeros(R_86.sizeMAP);   % just your tiff's x by y
        MAP(logical(R_86.pixel2keep(:,z))) = 1;               % x by y matrix
        pxp = R_86.pixel2keep(:,z);          % logical (x*y) array of significant pixels for this z plane
        
        idxChunk = idxTot(1:sum(pxp(:)));   % sum(pxp(:)) is the number of significant pixels for this z plane. Let's call this number SP.  Take those pixels out.
        idxTot(1:sum(pxp(:))) = [];     % idxTot stores the remaining clustering indices for the significant pixels of the remaining z planes.
        Lia = idxChunk == k;            % Lia si a logical, SP-dimensional array. Within the chunk of this z plane, take only the ones that belong to this specific cluster.
        numLia = find(Lia);  %numbers here still refer to SP and not to the tiff-dimensions.
        
        % next you go from a logical (x*y) array [pxp] to an (x*y) array
        % [pxpCS],
        % where significant pixels are numbered in ascending order and
        % everything else is a zero. These numbers will correspond to
        % numbers in numLia (which is in SP-format).
        pxpCS = cumsum(pxp(:));     
        pxpCSDif = diff(pxpCS);
        on = pxpCSDif == 1;
        on = [0;on];
        pxpCS(~on) = 0; %this has the same code as numLia now

        %find intersection with specific cluster, and store it in ROIs and
        %MAPs:
        finalLia = ismember(pxpCS, numLia);  
        clust(2).ROIs(:,k,z) = finalLia;
        if sum(finalLia) ~= 0
            MAP(finalLia)=10;
        else
            MAP(1)=10;
        end
        clust(2).MAPs(:,:,k,z) = MAP;
    end
end

%
    for k = 1 : nK
        figure('Color', [1,1,1], 'WindowStyle', 'normal', 'Position', [1  1  1440  820]); hold on
%         for z = 1 : size(R_100.pixel2keep,2) 
%             hax = subplot(ceil(NaZ/6),6,z); hold on
%             axis image
%             axis off
%             axis ij
%             imagesc(clust(1).MAPs(:,:,k,z)); % no colors
%         end
        for z = 1 : size(R_86.pixel2keep,2) 
%             hax = subplot(ceil(NaZ/6),6, size(R_100.pixel2keep,2)+z);
hax = subplot(ceil(44/6),6, z);
            hold on
            axis image
            axis off
            axis ij
            imagesc(clust(2).MAPs(:,:,k,z)); % no colors
        end
        
        fname = fullfile(sets.flyfolder, sprintf('manualclusts_cluster%02d_MAP_%s.pdf', k, identifierALL));
        export_fig(fname, '-m10')
        close 

    end

%% calculate or load ROIs
sets.flyfolder = fileparts(sets.analysislist{1});
sets.flyfolder = fileparts(sets.flyfolder);
sets.flyfolder = fileparts(sets.flyfolder);
sets.runs_basename = extractRunsBasename_comboFlies( sets, {aZ.runfolder},  {aZ.basename});
if z > 1
    sets.flyfolder = fullfile(sets.flyfolder, sets.runs_basename);
    if exist(sets.flyfolder,'dir') ~= 7
        mkdir(sets.flyfolder);
    end
else
    sets.flyfolder = aZ.runfolder;
end

importlist = [];

if sets.loadexistingclusters
    importlist = dir([sets.flyfolder '/*expROIdata_comboRuns_*' sets.runs_basename '*.mat']);
    if length(importlist) > 1
        str = {importlist.name};
        [selection,validation] = listdlg('PromptString','Select a file:',...
            'SelectionMode','single','ListSize', [500 100], ...
            'ListString',str);
        if validation
            importlist = importlist(selection);
        else
            error('Define the ROIs datafile...')
        end
    end
end

if sets.loadexistingclusters  && ~isempty(importlist)
    load(fullfile(sets.flyfolder, importlist.name));
else
    if sets.redefine_pix2use
%         indices = redefine_azIndices(aZ(1), sets);
        try
%         aZ(1) = parametricMapsF(aZ(1), indices);
        aZ = parametricMapsF_diffTrials(aZ, pUSE, nUSE, test2use);
        catch err
            if strcmp(err.message, 'Subscripted assignment between dissimilar structures.')
                clear aZ
                load_aZ_structure;
            end
        end
    end
    % new
    [M, pixel2keep, iKeep, sets.identifier_tSNE] = sortedOut_tSNE( aZ, sets );
%     if isempty(sets.prefix_nK)
%         [ROIs, MAPs, identifierALL, clusterData] = sortedOut_GMM( aZ, sets, M, pixel2keep, iKeep );
%     else
%         [ROIs, MAPs, identifierALL, clusterData] = sortedOut_GMM( aZ, sets, M, pixel2keep, iKeep, sets.prefix_nK );
%     end
    
    
    
% [pixel2keep, iKeep] = sortedOut_accessory( aZ, sets);   % just temp here
    if isempty(sets.prefix_nK)
        ReplicatesK = 80;
%         [j, m]   = sortedOut_GMM_selectBIC_JM(        sets.flyfolder, M );
j = 2;
m = 2;
        settings = sortedOut_GMM_selectBIC_K_compute( sets.flyfolder, M, j, m, ReplicatesK, 20 );
        settings.datenow = []; %empty it if you want to use all the available sets in the folder
        bK       = sortedOut_GMM_selectBIC_K_extract( sets.flyfolder, settings);
        clusterData = sortedOut_GMM_finalCount( sets, M, j, m, bK );
    else
        [j, m]   = sortedOut_GMM_selectBIC_JM(        sets.flyfolder, M, sets.prefix_nK );
        clusterData = sortedOut_GMM_finalCount( sets, M, j, m, sets.prefix_nK );
    end
    identifierALL = clusterData.identifier;
    sizeMAP = size(aZ(1).pixels2keep);
    [ROIs, MAPs, clusterData] = sortedOut_GMM( clusterData, sizeMAP, sets, M, pixel2keep); %using idx now
%     [ROIs, MAPs, clusterData] = sortedOut_GMM( clusterData, size(aZ(1).pixels2keep), sets, M, pixel2keep); %using idx now
    
%     [bootITR, clustSizes] = bootstrapITR(ROIs, aZ, sets.flyfolder, identifierALL, 50, 0, 500);
%     bootIRR = bootstrapIRR(ROIs, aZ, sets.flyfolder, identifierALL, 50, 0, 500);
%     clear importlist bK identifier_tSNE j k m
%     nK = size(ROIs,2);
%     identifierALL = num2str(identifierALL);
%     colors = distinguishable_colors(size(MAPs,3));
%     cMAPs = makeColoredMaps(MAPs, colors);      % convert maps in color
%     plotClusterMAPs_acrossZs;

    %%
    nK = size(ROIs,2);
    
%     nK = 10:20;
%     for inK = 1 : length(nK);
%         [gen(inK).ROIs, gen(inK).MAPs, gen(inK).identifier, gen(inK).clusterData] = sortedOut_GMM( aZ, sets, M, pixel2keep, iKeep, nK(inK) );
%         close all
%     end

end


%% concatenate stimuli (in their respective relevant window) and calculate psths

z = 1; %assuming now that the stimuli and other conditions are the same across all runs.

n_stimclasses = length(aZ(z).T);
load(fullfile(aZ(z).runfolder, [aZ(z).basename '_preprocessed.mat']), 'metaimage');
fs_dec = metaimage.acq.frameRate;
fs_acq = 4e4;

zeros_st_idxs = 1; %zeros(aZ(z).NStimuli+1,2);
stimTraces = [];
ts_Stim_all = [];
divs_classes_idxs = 1;
for stCl = 1 : n_stimclasses
    for n = 1 : length(aZ(z).T(stCl).TrTypesThisClass)
%         sT = aZ(z).T(stCl).sensor(n,:);
%         if stCl*n == 1
%             DCoffset = mean(sT(20:2*4e4));
%         end
        idx = find(aZ(z).T(stCl).trials == n, 1);
        sT = aZ(z).ALLstimuli(idx).stim.makeDigital;
        
        sT = sT(iKeep(z).stClass(stCl).Response_ts);
        sT = sT(:);
        stimTraces = cat(1,stimTraces,sT);
        
        ts_stim = aZ(z).T(stCl).ts(iKeep(z).stClass(stCl).Response_ts);
        if isempty(ts_Stim_all)
            ts_Stim_all = ts_stim;
        else
            ts_stim = ts_stim + ts_Stim_all(end) + 1/fs_acq;
            ts_Stim_all = cat(2, ts_Stim_all, ts_stim);
        end
        zeros_st_idxs = cat(1, zeros_st_idxs, length(ts_Stim_all));
    end
    divs_classes_idxs = cat(1, divs_classes_idxs, n);
end
% stimTraces = (stimTraces-DCoffset) * 3; %um conversion %3 for new 30um piezo

zeros_st = ts_Stim_all(zeros_st_idxs);
divs_classes_idxs = cumsum(divs_classes_idxs);

% divs_classes_idxs = unique([divs_classes_idxs', 26,42,58,74,90]);     %ad hoc %WEDGE_full
% divs_classes_idxs = unique([divs_classes_idxs', 9:7:44]);         %ad hoc?
% divs_classes_idxs = unique([divs_classes_idxs', 11:7:44]);        %ad hoc JONs final
divs_classes_log = zeros(size(zeros_st));
divs_classes_log(divs_classes_idxs) = 1;

divs = divisors(aZ(z).NStimuli);
[~,idx] = min(abs(divs - 5.5));
divs = 1; %divs(idx);


%%
% define the cluster order for once.
clust_use=1:nK;
outperm = clust_use;
% 
% if nK >1
%     % %either:
%     [ T, outperm ] = sortClusters( MAPs(:,:,clust_use,z));
%     close
%     % %or, manually, but unshuffle first:
%     % outperm = [2 13 15 18 17 14 19 11 8 4 9 6 1 10 16 7 5 12 3 20];
%     % outperm = [11 12 14 6 8 4 15 20 16 19 10 5 7 2 9 18 13 1 3 17];
%     % [~, outperm] = sort(outperm);
% 
%     clust_use = clust_use(outperm);
% else
%     outperm = 1;
% end

% colors = brewermap(nK,'RdYlGn');



%% START HERE
% colors = distinguishable_colors(4);
% colors = [166,206,227; ...
%           31,120,180; ...
%           178,223,138; ...
%           51,160,44]./255;
colors = [208,28,139; ...   %fucsia
          31,120,180; ...   %blu
          230,97,1; ...     % orange
          51,160,44]./255;  % green
      

for z = 1 : length(sets.analysislist)
    close all
    
    ts_dec = aZ(z).T(stCl).ts_dec;
    for k = 1:nK
        traces_DFFs_singK = [];
        zeros_stims = 0;
        for stCl = 1 : n_stimclasses   % !!
            % sort stimuli (all, individually)%assuming the same number of
            % repets, otherwise change
            t = squeeze(reshape(aZ(z).T(stCl).tiff, size(aZ(z).T(stCl).tiff,1)*size(aZ(z).T(stCl).tiff,2),1,size(aZ(z).T(stCl).tiff,3), size(aZ(z).T(stCl).tiff,4)));
            t_baseline = nanmean(squeeze(nanmean(t( logical(ROIs(:,k,z)),:,aZ(z).T(stCl).indices.ibaseline),1)),2);
            t = squeeze(nanmean(t( logical(ROIs(:,k,z)),:,:),1));
            t_dFF = bsxfun(@rdivide, (bsxfun(@minus,t,t_baseline)), t_baseline).*100;

            for i = 1 : length(aZ(z).T(stCl).TrTypesThisClass)
                clust(k).traces(i).traces = t_dFF(aZ(z).T(stCl).iTrials(:,i), :)'; %npointsinResponse x nRepetitions
                clust(k).traces(i).mean = nanmean(clust(k).traces(i).traces,2);
                clust(k).traces(i).std = nanstd(clust(k).traces(i).traces,0,2);
            end
        end
    end
    
    %% plot
    relWidthsColumns = [0.89, 0.11];
    [hfig, hax] = figureTracesI_PP( nK+1, relWidthsColumns );
    hfig.Name = 'Selected clusters';
    hfig.Position = [1.0000   -0.0578    0.4243    1.0311];
    tbox = annotation('textbox');
    tbox.Position = [0.7667    0.9462    0.1600    0.0500];
    tbox.EdgeColor = [1,1,1];
    tbox.String = sprintf('%s\n%s\nID: %d',aZ(z).basename, sets.runs_basename, identifierALL);
    tbox.Interpreter = 'none';
    tbox.HorizontalAlignment = 'center';
    tbox.FontSize = 13.5;
    
    axes(hax(1,1))
    plot([0,unique(aZ(z).fastStimMeta.stimDur)],[0.5,0.5],'-k', 'LineWidth',4)
    xlim([ts_dec(1),ts_dec(end)])
    hax(1,1).XRuler.TickLabel = ceil(ts_dec(1)):1:floor(ts_dec(end));
    hax(1,1).XTick = [ceil(ts_dec(1)):1:floor(ts_dec(end))];
    hax(1,1).XGrid = 'off';
    hax(1,1).YAxis.Visible = 'off';
    hax(1,1).XAxisLocation = 'top';
    hax(1,1).XLabel.String = 'time (s)';
    
    for k = 1:nK
        axes(hax(k+1,1))
        hold on
        for i = 1 : 4 %length(aZ(z).T(stCl).TrTypesThisClass)
%             errorbar(traces(i).mean, traces(i).std);
            
            plot(ts_dec,  clust(k).traces(i).traces, 'LineWidth',0.5, 'Color', colors(i,:).*0.75)
            plot(ts_dec, clust(k).traces(i).mean,'LineWidth',3, 'Color', colors(i,:))
        end
        axis tight
    end
    
    
    
    %%
    
    
    %% 3. add cluster map and reliability index
    %1: implement inter-trial reliability measure as the mean of all pair-wise
    %correlations between trials
%     warning('check this line')
    nTrialsMin = size(clust(1).traces(1).traces, 2); %check, across z
    idxs = ~triu(ones(nTrialsMin));
    for k = 1:nK
        traces_DFFs= [];
        for ist = 1:4
            traces_DFFs = cat(1, traces_DFFs, clust(k).traces(ist).traces);
        end
        pwcs = corr(traces_DFFs);
        intertrialRel(k) = mean(mean(pwcs(idxs)));
    end
%     
    for k = 1:nK
        axes(hax(k+1,length(relWidthsColumns)))
        hax(k+1,length(relWidthsColumns)).YAxisLocation = 'right';
        hax(k+1,length(relWidthsColumns)).YLabel.HorizontalAlignment = 'left';
        imagesc(MAPs(:,:,clust_use(k),z));
        axis image
        axis off
        ylabel(sprintf('%1.2f',intertrialRel(clust_use(k))));
        hax(k+1,length(relWidthsColumns)).YLabel.Visible = 'on';
        hax(k+1,length(relWidthsColumns)).YLabel.Rotation = 0;
        hax(k+1,length(relWidthsColumns)).YDir = 'reverse';
    end
    k = 1;
    hax(k+1,length(relWidthsColumns)).YLabel.HorizontalAlignment = 'left';
    hax(k+1,length(relWidthsColumns)).YLabel.String = sprintf('inter-trial\nreliability\n%1.2f',intertrialRel(clust_use(k)));
    hax(k+1,length(relWidthsColumns)).YLabel.VerticalAlignment = 'bottom';
    
    delete(hax(1,2))
    
    %% save fig
    figure(hfig)
    name = fullfile(sets.flyfolder, sprintf('analysis_%s_%d',sets.runs_basename, identifierALL));
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
    
%     %% colored map
%     % draw map with colored ROIs
%     
%     mapClean = zeros(size(MAPs, 1), size(MAPs,2));
%     for i = [1:nK] %index of colors
%         j = outperm(i);     %index of clusters
%         mapClean(MAPs(:,:,j,z)==10) = i; %i.e. when i==1, takes the cluster in first position and gives it color # 1.
%     end
%     Col = [1,1,1];
%     Col = cat(1, Col, colors(1:nK,:)); 
%    
%     
%     figure('Color', [1,1,1]);
%     imagesc(mapClean)
%     axis off
%     axis image
%     colormap(Col)
%     title(sprintf('%s\n%s\n%d',aZ(z).basename, sets.runs_basename, identifierALL), 'Interpreter', 'none')
%     
%     % save fig
%     name = fullfile(sets.flyfolder, sprintf('ColoredMap_%s_%d',sets.runs_basename, identifierALL));
%     % do not overwrite
%     d = dir([name,'*']);
%     name = sprintf('%s_%d.png', name, length(d) );
%     % save
%     export_fig(name, '-m3.5')
%     % export_fig(name)
%     
    
    

end


