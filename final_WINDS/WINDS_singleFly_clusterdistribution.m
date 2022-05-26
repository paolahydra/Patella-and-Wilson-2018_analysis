% load
% based on: clustering_on_rawData_selectWINDS_FLIPPED
% ---> in folder: analysis/final_WINDS
Folder2Save = '/Users/galileo/Dropbox (HMS)/Data/WINDS_100/Downsampled/SubsetFlies';
load(fullfile(Folder2Save, 'datalist_selectWindFlies.mat'));
load(fullfile(Folder2Save, 'R_matrix_downsampled_winds_100&86_withBaseline.mat'));
R = Rflipped;
clear Rflipped
clear runfolders
for z = 1:length(datalist)
    [runfolders{z}, ~] = fileparts(datalist{z});
end

clusterFolder = fullfile(Folder2Save, '/rawDataLinkage'); % linkage was already calculated on orchestra
load(fullfile(clusterFolder, 'output_ward_Rzsp_linkage_winds_FLIPPED.mat'), 'Zc');
cd(clusterFolder)
pxKeep = true(size(R, 1),1);   % no (easily-detected at least) uncorrelated pixels within this dataset

load('/Users/galileo/Dropbox (HMS)/Data/cmap19_spectral.mat', 'cmapwinds')

%% select a K and plot subclusters (traces)
K = 8;
nK = K;
load(fullfile(clusterFolder, sprintf('klust_maxClust%d.mat', K)));   %'klust', 'dendrOrder', 'Tparent'); % these are sorted as in dedndrogram

cmap = cmapwinds(1:8,:);

orDT = load('/Users/galileo/Dropbox (HMS)/Data/WINDS_100/datalist_allWINDS.mat', 'dataTable');
orDT = orDT.dataTable;
clust = remapIndexedPointsToRun_RoisMaps_oldNew(Tparent, R_100, R_86);      %sukee
% make maps (divided old/new because different dimensions)
klust = makeplot_cMaps_superKlusts_singleIteration_oldNew(clust, klust, dendrOrder, cmap);
regGeneralData = matfile('/Users/galileo/Dropbox (HMS)/Data/WINDS_100/anatomyImages/anatomy_alignment_metadata.mat');
allStacks = regGeneralData.allStacks; %allStacks
allStacks_flyNums = cat(1,allStacks.flyNum);

dWA = matfile('/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/datalist_WEDAMMC_piezo.mat');
alignedAMMCmasks = dWA.alignedAMMCmasks;


colors = [205, 149 43; ...
    148, 111, 176; ...
    237, 28, 48; ...
    55, 160, 72]./255;


%% for each fly, for each cluster, plot  zs(mean sustained phase) for all 4 stimuli
load('Rzs_winds_downsampled_13FLIPPED.mat')
Rzs_sc = reshape(permute(Rzs,[2,3,1]),[], 4, size(Rzs,1) );
Rzs_sc = squeeze(mean(Rzs_sc));

% make corresponding version of R with mean(sustained phase only(dff)) --> R2
load(fullfile(Folder2Save, 'R_matrix_downsampled_winds_100&86.mat'), 'Rflipped'); %no baseline here, to use same indices
R2 = reshape(Rflipped, size(Rflipped,1), [], 4);
R2 = R2(:,19:41,:); 
% R2 = reshape(R2, size(R2,1), [], 1);
R2 = (squeeze(mean(R2,2)))';


flyNums = unique(dataTable.fly, 'stable');



px2k = R_100.pixel2keep(:);
%split Zsc and Tparent
Rzs_sc_100 = Rzs_sc(:, 1:sum(px2k));
R2_100 = R2(:, 1:sum(px2k));

Rzs_sc_86 = Rzs_sc;
Rzs_sc_86(:, 1:sum(px2k))=[];
R2_86 = R2;
R2_86(:, 1:sum(px2k))=[];

Tparent_100 = Tparent(1:sum(px2k));
Tparent_86 = Tparent;
Tparent_86(1:sum(px2k))=[];


%% zscored bar plots
[hfig, hax] = figureTracesI_PP_pixels( K, repmat(80,1,length(flyNums)), 160 );  
hax = flipud(hax);
f = 1;    
for Ki = 1:K
    allRzs_sc = Rzs_sc_100(:,Tparent_100==Ki);
    
    meanZS_perStim = mean(allRzs_sc');
    axes(hax(Ki,f)), hold on
    bar(meanZS_perStim); %no specific colors in 2016b
    hax(Ki,f).YLim = [-2, 2];
    ylabel(Ki)
end
hax(Ki,f).Title.String = flyNums(f);


Rzs_sc_86_remainder = Rzs_sc_86;
Tparent_86_remainder = Tparent_86;
for f = 2:4;  
    zetas_thisfly = dataTable.fly == flyNums(f);
    zetas_thisfly(1:NaZ_100) = [];
    px2kThisFly = sum(sum(R_86.pixel2keep(:,zetas_thisfly) ));
    
    
    
    Rzs_sc_86_thisfly = Rzs_sc_86_remainder(:,1:px2kThisFly);
    Rzs_sc_86_remainder(:,1:px2kThisFly) = [];
    
    Tparent_86_thisfly = Tparent_86_remainder(1:px2kThisFly);
    Tparent_86_remainder(1:px2kThisFly) = [];
    
    for Ki = 1:K
        allRzs_sc = Rzs_sc_86_thisfly(:,Tparent_86_thisfly==Ki);
        meanZS_perStim = mean(allRzs_sc');
        axes(hax(Ki,f)), hold on
        bar(meanZS_perStim); %no specific colors in 2016b
        hax(Ki,f).YLim = [-2, 2];
        axis off
    end
    hax(Ki,f).Title.String = flyNums(f);
    
end

export_fig(sprintf('bars_zs_sust_means_byFly.eps'))


%% dff bar plots
[hfig, hax] = figureTracesI_PP_pixels( K, repmat(80,1,length(flyNums)), 160 );  
hax = flipud(hax);
f = 1;    
for Ki = 1:K
    allRzs_sc = R2_100(:,Tparent_100==Ki);
    
    meanZS_perStim = mean(allRzs_sc');
    axes(hax(Ki,f)), hold on
    bar(meanZS_perStim); %no specific colors in 2016b
%     hax(Ki,f).YLim = [-2, 2];
    ylabel(Ki)
end
hax(Ki,f).Title.String = flyNums(f);


R2_86_remainder = R2_86;
Tparent_86_remainder = Tparent_86;
for f = 2:4;  
    zetas_thisfly = dataTable.fly == flyNums(f);
    zetas_thisfly(1:NaZ_100) = [];
    px2kThisFly = sum(sum(R_86.pixel2keep(:,zetas_thisfly) ));
    
    
    
    R2_86_thisfly = R2_86_remainder(:,1:px2kThisFly);
    R2_86_remainder(:,1:px2kThisFly) = [];
    
    Tparent_86_thisfly = Tparent_86_remainder(1:px2kThisFly);
    Tparent_86_remainder(1:px2kThisFly) = [];
    
    for Ki = 1:K
        allRzs_sc = R2_86_thisfly(:,Tparent_86_thisfly==Ki);
        meanZS_perStim = mean(allRzs_sc');
        axes(hax(Ki,f)), hold on
        bar(meanZS_perStim); %no specific colors in 2016b
        axis off
    end
    hax(Ki,f).Title.String = flyNums(f);
    
end

yaxLims = [min(min(cat(1,hax.YLim))), max(max(cat(1,hax.YLim)))];
hax(1,1).YLim = yaxLims;
linkaxes(hax,'y')

export_fig(sprintf('bars_dff_sust_means_byFly.eps'))

%% boxcar - dff
[hfig, hax] = figureTracesI_PP_pixels( K, repmat(80,1,length(flyNums)), 160 );  
hax = flipud(hax);
f = 1;    
for Ki = 1:K
    axes(hax(Ki,f)), hold on
    allRzs_sc = R2_100(:,Tparent_100==Ki);
%     meanZS_perStim = mean(allRzs_sc');
%     bar(meanZS_perStim); %no specific colors in 2016b
    boxplot(allRzs_sc', 'Colors', colors, 'Widths',0.6, 'OutlierSize', 2)
    plot(1:4, mean(allRzs_sc'), 'xk')
    ylabel(Ki)
    hax(Ki,f).XGrid = 'off';
    hax(Ki,f).XAxis.Visible = 'off';
    hax(Ki,f).Box = 'off';
end
hax(Ki,f).Title.String = flyNums(f);


R2_86_remainder = R2_86;
Tparent_86_remainder = Tparent_86;
for f = 2:4;  
    zetas_thisfly = dataTable.fly == flyNums(f);
    zetas_thisfly(1:NaZ_100) = [];
    px2kThisFly = sum(sum(R_86.pixel2keep(:,zetas_thisfly) ));
    
    
    R2_86_thisfly = R2_86_remainder(:,1:px2kThisFly);
    R2_86_remainder(:,1:px2kThisFly) = [];
    
    Tparent_86_thisfly = Tparent_86_remainder(1:px2kThisFly);
    Tparent_86_remainder(1:px2kThisFly) = [];
    
    for Ki = 1:K
        axes(hax(Ki,f)), hold on
        allRzs_sc = R2_86_thisfly(:,Tparent_86_thisfly==Ki);
        boxplot(allRzs_sc', 'Colors', colors, 'Widths',0.6, 'OutlierSize', 2)
        plot(1:4, mean(allRzs_sc'), 'xk')
        axis off
    end
    hax(Ki,f).Title.String = flyNums(f);
    
end

yaxLims = [min(min(cat(1,hax.YLim))), max(max(cat(1,hax.YLim)))];
hax(1,1).YLim = yaxLims;
linkaxes(hax,'y')

export_fig(sprintf('boxplot_dff_sust_means_byFly.eps'))

%% singleTrials lines - dff
[hfig, hax] = figureTracesI_PP_pixels( K, repmat(80,1,length(flyNums)), 160 );  
hax = flipud(hax);
f = 1;    
for Ki = 6%1:K
    axes(hax(Ki,f)), hold on
    allRzs_sc = R2_100(:,Tparent_100==Ki);
%     meanZS_perStim = mean(allRzs_sc');
%     bar(meanZS_perStim); %no specific colors in 2016b
%     boxplot(allRzs_sc', 'Colors', colors, 'Widths',0.6, 'OutlierSize', 2)
%     plot(1:4, mean(allRzs_sc'), 'xk')
    plot(allRzs_sc, '-k')
    ylabel(Ki)
    hax(Ki,f).XGrid = 'off';
    hax(Ki,f).XAxis.Visible = 'off';
    hax(Ki,f).Box = 'off';
end
hax(Ki,f).Title.String = flyNums(f);


R2_86_remainder = R2_86;
Tparent_86_remainder = Tparent_86;
for f = 2:4;  
    zetas_thisfly = dataTable.fly == flyNums(f);
    zetas_thisfly(1:NaZ_100) = [];
    px2kThisFly = sum(sum(R_86.pixel2keep(:,zetas_thisfly) ));
    
    
    R2_86_thisfly = R2_86_remainder(:,1:px2kThisFly);
    R2_86_remainder(:,1:px2kThisFly) = [];
    
    Tparent_86_thisfly = Tparent_86_remainder(1:px2kThisFly);
    Tparent_86_remainder(1:px2kThisFly) = [];
    
    for Ki = 6% 1:K
        axes(hax(Ki,f)), hold on
        allRzs_sc = R2_86_thisfly(:,Tparent_86_thisfly==Ki);
%         boxplot(allRzs_sc', 'Colors', colors, 'Widths',0.6, 'OutlierSize', 2)
%         plot(1:4, mean(allRzs_sc'), 'xk')
        plot(allRzs_sc, '-k')
        axis off
    end
    hax(Ki,f).Title.String = flyNums(f);
    
end

yaxLims = [min(min(cat(1,hax.YLim))), max(max(cat(1,hax.YLim)))];
hax(1,1).YLim = yaxLims;
linkaxes(hax,'y')

export_fig(sprintf('lines_dff_cl8>6_sust_means_byFly.eps'))

%% singleTrials lines - zsc
[hfig, hax] = figureTracesI_PP_pixels( K, repmat(80,1,length(flyNums)), 160 );  
hax = flipud(hax);
f = 1;    
for Ki = 6%[3,4,6,8]
    axes(hax(Ki,f)), hold on
    allRzs_sc = Rzs_sc_100(:,Tparent_100==Ki);
%     meanZS_perStim = mean(allRzs_sc');
%     bar(meanZS_perStim); %no specific colors in 2016b
%     boxplot(allRzs_sc', 'Colors', colors, 'Widths',0.6, 'OutlierSize', 2)
%     plot(1:4, mean(allRzs_sc'), 'xk')
    plot(allRzs_sc, '-k')
    ylabel(Ki)
    hax(Ki,f).XGrid = 'off';
    hax(Ki,f).XAxis.Visible = 'off';
    hax(Ki,f).Box = 'off';
end
hax(Ki,f).Title.String = flyNums(f);


Rzs_sc_86_remainder = Rzs_sc_86;
Tparent_86_remainder = Tparent_86;
for f = 2:4;  
    zetas_thisfly = dataTable.fly == flyNums(f);
    zetas_thisfly(1:NaZ_100) = [];
    px2kThisFly = sum(sum(R_86.pixel2keep(:,zetas_thisfly) ));
    
    
    Rzs_sc_86_thisfly = Rzs_sc_86_remainder(:,1:px2kThisFly);
    Rzs_sc_86_remainder(:,1:px2kThisFly) = [];
    
    Tparent_86_thisfly = Tparent_86_remainder(1:px2kThisFly);
    Tparent_86_remainder(1:px2kThisFly) = [];
    
    for Ki = 6%1:K
        axes(hax(Ki,f)), hold on
        allRzs_sc = Rzs_sc_86_thisfly(:,Tparent_86_thisfly==Ki);
%         boxplot(allRzs_sc', 'Colors', colors, 'Widths',0.6, 'OutlierSize', 2)
%         plot(1:4, mean(allRzs_sc'), 'xk')
        plot(allRzs_sc, '-k')
        axis off
    end
    hax(Ki,f).Title.String = flyNums(f);
    
end

yaxLims = [min(min(cat(1,hax.YLim))), max(max(cat(1,hax.YLim)))];
hax(1,1).YLim = yaxLims;
linkaxes(hax,'y')

export_fig(sprintf('lines_zsc_cl8>6_sust_means_byFly.eps'))

