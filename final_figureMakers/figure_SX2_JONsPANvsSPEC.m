%% load
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


load('/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/anatomyImages/JON_cropMask.mat')

%% reference final clustering
cmap10H = flipud(cmap10);
figure; imagesc(1:size(cmap10H,1)); colormap(cmap10H); title('resorted cmap10 - horizontal')

K1 = 10;
T1 = load(fullfile(clusterFolder, sprintf('klust_maxClust%d.mat', K1)));
load(sprintf('singleKmaps_maxClust%d.mat', K1), 'singleKmaps')
T1.singleKmaps = singleKmaps;
clear singleKmaps
% sort these by color/finalAXIS
cOrderHoriz = fliplr(clusterOrder); %black left, A right
TR.Tparent = T1.Tparent;

for Ki = 1:K1
    K2use = cOrderHoriz(Ki);
    TR.Tparent(T1.Tparent==K2use) = Ki;
    TR.klust(Ki) = T1.klust(K2use);
end
TR.singleKmaps = T1.singleKmaps(:,:,:,cOrderHoriz);
% for Ki = 1:10, figure; imagesc(TR.singleKmaps(:,:,13,Ki)), axis image, title(Ki), end

% OK!!
clear T1
nK = 10;

%% tabulate genotypes
genotypes = unique(dataTable.genotype);
for g = 1 : length(genotypes)
    tabulateGenotypes(:,g) = dataTable.genotype== genotypes(g);
end

%% do average maps by genotype, clusters and z levels as in figure 3 (figure2n) (in cropped region)
maps = zeros(size(cropMask, 1), size(cropMask, 2), nK, length(genotypes), 4);  %   1,2:map;    3:K;    4:genotype;    5:depth
for g = 1 : length(genotypes)
    for Z = 1:4
        zetas = tabulateGenotypes(:,g) & dataTable.zLevel==Z;
        for Ki = 1:nK
            maps(:,:, Ki, g, Z)  = nansum(TR.singleKmaps(:,:,zetas,Ki) , 3)  .* cropMask;
%             maps(:,:, Ki, g, Z)  = (nansum(TR.singleKmaps(:,:,zetas,Ki) , 3) ./ sum(zetas)) .* cropMask;
        end
    end
end
maps(~isfinite(maps)) = 0;
maps(maps>1) = 1;

%% focus on level 3 only for now

%% correlate images (within and between clusters as a baseline) and plot each pmax. May also calculate and plot confidence intervals
Z = 3;
r3 = zeros(nK, nK, length(genotypes));

for Ki = 1:nK
    for Ki2 = 1 : nK
        for g2 = 1 : length(genotypes)
            r3(Ki, Ki2, g2) = corr2( maps(:,:,Ki2,1,Z), maps(:,:,Ki,g2,Z) );
        end
    end
end


for g2 = 1 : length(genotypes)
    figure;
    imshow(r3(:,:,g2))
    axis on
    box off
    ylabel(sprintf('genotype %d\nclusters:',genotypes(g2)))
    xlabel('corr relative to PAN clusters:')
    ax = gca;
    ax.XAxisLocation = 'top';
end
    
%     
% [hfig, hax] = figureTracesI_PP(10, 0.1*ones(1,10), 0.01);
% for Ki = 1:nK
%     for g1 = 1 : length(genotypes)
%         axes(hax(Ki, g1))
%     end
% end

%% Xcorrelate images and take peak
Z = 3;
rX3 = zeros(nK, nK, length(genotypes));

for Ki = 1:nK
    for Ki2 = 1 : nK
        for g2 = 1 : length(genotypes)
            C = normxcorr2( maps(:,:,Ki2,1,Z), maps(:,:,Ki,g2,Z) );
            [ypeak, xpeak] = find(C==max(C(:)));
            yoffSet = ypeak-size(maps,1);
            xoffSet = xpeak-size(maps,2);
            rX3(Ki, Ki2, g2) = max(C(:));
        end
    end
end


for g2 = 1 : length(genotypes)
    figure;
    imshow(rX3(:,:,g2))
    axis on
    box off
    ylabel(sprintf('genotype %d\nclusters:',genotypes(g2)))
    xlabel('XCORR max - relative to PAN clusters:')
    ax = gca;
    ax.XAxisLocation = 'top';
end




