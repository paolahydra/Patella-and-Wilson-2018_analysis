%% flip stimuli 1-3 for L pixels only
% Save zs2 Rmatrix for linkage calculation
load(fullfile(Folder2Save, 'R_matrix_downsampled_unfiltered_winds_100&86.mat'));
% R_100.pixel2keep = load('/Users/galileo/Dropbox (HMS)/Data/WINDS_100/Downsampled/R_matrix_Downsampled_100.mat', 'pixel2keep');
% R_86.pixel2keep = load('/Users/galileo/Dropbox (HMS)/Data/WINDS_100/Downsampled/R_matrix_Downsampled_86.mat', 'pixel2keep');
% R_100.pixel2keep = R_100.pixel2keep.pixel2keep;
% R_86.pixel2keep = R_86.pixel2keep.pixel2keep;
% save(fullfile(Folder2Save, 'R_matrix_downsampled_winds_100&86.mat'), 'R_100', 'R_86', '-append');

R2 = reshape(R, size(R,1), [], 4);
pixels2keep_allZs = cat(2, sum(R_100.pixel2keep), sum(R_86.pixel2keep)); %cool
pix2k_allZs_end = cumsum(pixels2keep_allZs);
pix2k_allZs_start = pix2k_allZs_end+1;
pix2k_allZs_start = [1, pix2k_allZs_start(1:end-1)];


for z = 1:NaZ %easy approach, one z at the time
    if ~dataTable.ipsilateral(z) %if 0, flip
        R2(pix2k_allZs_start(z):pix2k_allZs_end(z),:,:) = R2(pix2k_allZs_start(z):pix2k_allZs_end(z),:,[3,2,1,4]);
    end
end
Rflipped = R2;
Rflipped = reshape(Rflipped, size(Rflipped,1), [], 1);
save(fullfile(Folder2Save, 'R_matrix_downsampled_unfiltered_winds_100&86.mat'), 'Rflipped', '-append');

%% only consider sustained response phase:
% R = reshape(R, size(R,1), [], 4);
m = matfile(fullfile(Folder2Save, 'R_matrix_SubsetFlies_smallWindow.mat'));
T32 = aZ{32}.T;
ts_dec = T32.ts_dec;
ts_dec(ts_dec< m.windowOnset_fromStOnset ) = [];
ts_dec(ts_dec> T32.stimDur + m.windowOffset_fromStOffset ) = [];
clear m T32
save(fullfile(Folder2Save, 'R_matrix_downsampled_winds_100&86.mat'), 'ts_dec', '-append');

zeroidx = find(ts_dec>=0, 1);
tsKeep_i = false(size(ts_dec));
tsKeep_i(1:zeroidx-1) = true;
tsKeep_i(zeroidx+19-1:zeroidx+41-1) = true;

R2 = R2(:,tsKeep_i,:);     %19:41  %12:44
R2 = reshape(R2, size(R2,1), [], 1);
Rzs = zscore(R2, 0, 2);

% now get rid of the baseline for linkage purpose
R2 = reshape(Rzs, size(Rzs,1), [], 4);
R2(:,1:zeroidx-1,:) = [];
R2 = reshape(R2, size(R2,1), [], 1);
Rzs = R2;


clusterFolder = fullfile(Folder2Save, '/rawDataLinkage'); 
mkdir(clusterFolder)
% want to clean it up a bit? %maybe not
save(fullfile(Folder2Save,'Rzs_winds_downsampled_13FLIPPED.mat'), 'Rzs');

%%
% figure; hold on
% for i = 1:4
%     plot(mean(Rflipped(:,:,i)));
% end
% 
% Rzs = reshape(Rzs, size(Rzs,1), [], 1);
% figure; hold on
% for i = 1:4
%     plot(mean(Rzs(:,:,i)));
% end