% need to revert some changes made by error in the conversion to PANSPEC
% mode

% clusterFolder = '/Users/galileo/Dropbox (HMS)/Data/JON_PAN&SPEC/Downsampled/rawDataLinkage/ward_movAVG_ZSCpix';
% treeFileName = 'output_ward_movAVG_ZSCpix_linkage.mat';
load(fullfile(clusterFolder, treeFileName), 'Zc');
sizeMap = [60, 86];

for K = 9:15
    disp(K)
    clearvars -except R Rzs K pixel2keep iKeep pxKeep total_within Zc clusterFolder sizeMap datalist flyfoldersUnique flyNumUnique
    NaF = length(flyfoldersUnique);
    maxClust = K;
    maxLeavesDend = maxClust;
    distanceMetrics = 'euclidean';
    linkageMethod = 'ward';
    % temporary dendr
    figure; [~,Tdend,outpermLeaves] = dendrogram(Zc, maxLeavesDend, 'Orientation', 'left');
    axDend = gca;
    axDend.YTickLabel = 1:maxLeavesDend;
    export_fig(sprintf('dendrogram_maxLeaves%d.pdf',maxLeavesDend))

    %% cluster (main), sort klusts based on dendrogram -- AND THEN RE-SORT again [DO]
    Tparent = cluster(Zc,'maxclust',maxClust);
    for Ki = 1:length(unique(Tparent))
        tempKs = find(Tparent==Ki);
        klust(Ki).k = (tempKs(:))';
    end
    for Ki = 1:length(unique(Tparent)), allElements(Ki, 1:length(klust(Ki).k)) = klust(Ki).k; end
    
    countLeaves = 1;
    for Ki = 1:length(unique(Tparent))   %while ~isempty(leavesorder)
        lengthLeave = 0;
        leavesorder = find(Tdend == outpermLeaves(countLeaves));
        countLeaves = countLeaves+1;
        lengthLeave = lengthLeave + length(leavesorder);
        firstLeaf(Ki) = leavesorder(1);
        [dendrOrder(Ki), ~] = find(allElements==firstLeaf(Ki));
        while sum(Tparent == dendrOrder(Ki)) > lengthLeave
            leavesorder = find(Tdend == outpermLeaves(countLeaves));
            countLeaves = countLeaves+1;
            lengthLeave = lengthLeave + length(leavesorder);
        end
    end
    clear leavesorder countLeaves leavesorder
    
    % Also rename T once and for all according to dendrogram
    T = zeros(size(Tparent));
    for Ki = 1:length(unique(Tparent))
        T(Tparent==dendrOrder(Ki)) = Ki;
        tempKs = find(T==Ki);
        klustSorted(Ki).k = (tempKs(:))';
    end
    Tparent = T;
    klust = klustSorted;
    dendrOrder = 1:length(unique(Tparent));
    clear klustSorted T
    clear outpermLeaves %to avoid errors
    nK = length(unique(Tparent));
    colors = distinguishable_colors(nK);
    save(fullfile(clusterFolder, sprintf('klust_maxClust%d.mat', maxClust)), 'Tparent', 'dendrOrder', 'klust')


%     load(fullfile(clusterFolder, sprintf('klust_maxClust%d.mat', maxClust)))
    nK = length(unique(Tparent));
    
    %% plot all average traces with errors, sorted from bottom to up as in dendrogram! [may skip if already done]
    [hfig, hax] = figureTracesI_PP( length(unique(Tparent)), 1 );
    hax = flipud(hax); % index 1 is bottom axis
    % hfig.Position   =   [1,5,1800,700];
    for Ki = 1:length(unique(Tparent))
        axes(hax(Ki)), hold on
        %     hax(Ki).XTick = 0:38:228;
        hax(Ki).XGrid = 'on';
        hax(Ki).FontSize = 3;
        
        allR = R(Tparent==Ki,:);
        avg = mean(allR,1);
        upperbound = prctile(allR, 75, 1);
        lowerbound = prctile(allR, 25, 1);
        plot(avg, '-r'), hold on, axis tight
        plot(upperbound, '-k')
        plot(lowerbound, '-k')
        axis tight
        ylabel(Ki)
        ylabel(sprintf('clust %d\n(px %d)',Ki, sum(Tparent==Ki)))
        hax(Ki).YLabel.FontSize = 6.5;
    end
%     hax(nK).Title.String = sprintf('total within distance: %.2e',total_within(nK));
%     hax(nK).Title.FontSize = 10;
    hfig.Position = [ 0    0.0044    0.3257    0.8900];
    export_fig(fullfile(clusterFolder, sprintf('clustMain_%d_mtraces.pdf', nK)))
    close
    
    %%
    maxClustComponents = 150;
    folder = fullfile(clusterFolder, sprintf('%dclusters', nK));
    if exist(folder) ~= 7
        mkdir(folder);
    end
    cd(folder)
    
    Tcomponents = cluster(Zc,'maxclust',maxClustComponents);
    [table_T10030] = crosstab(Tcomponents, Tparent);
    clust30 = remapIndexedPointsToRun_RoisMaps(Tparent, pixel2keep, sizeMap, pxKeep);
    clust100 = remapIndexedPointsToRun_RoisMaps(Tcomponents, pixel2keep, sizeMap, pxKeep);
    sortFliesGenoType = 0;
    [~, flyClusterCountDivZs100] = barcluster_panSpec(Tcomponents, datalist, flyfoldersUnique, clust100.ROIs, sortFliesGenoType);
    [~, flyClusterCountDivZs30] = barcluster_panSpec(Tparent, datalist, flyfoldersUnique, clust30.ROIs, sortFliesGenoType);
    
    figure; imagesc(table_T10030)
    axis square
    ax = gca;
    ax.XTick = 1:nK;
    ax.XTickLabel = sum(table_T10030);
    ax.XTickLabelRotation = 90;
    export_fig('table_main_components.png')
    %% display component traces
    % sortedflyMeta = {'pan' 'pan' 'pan' 'pan' 'pan' 'pan' 'A22' 'A23' 'A23' 'A26' 'AB15' 'AB15' 'AB15' 'BA28' 'B2' 'B2' 'ACE4' 'ACE4' 'ACE4' 'CE32' 'CE32' 'AD29'  'AD29' 'AD29'};
    for k = 1:nK
        % [hfig, hax] = figureTracesI_PP( sum(table_T10030(:,k)>0)+1, [0.66 0.34] );
        [hfig, hax] = figureTracesI_PP_pixels( sum(table_T10030(:,k)>0)+1, [600 300], 30 );
        axes(hax(1,1)), hold on
        allR = R(Tparent==k,:);
        avg = mean(allR,1);
        upperbound = prctile(allR, 75, 1);
        lowerbound = prctile(allR, 25, 1);
        plot(avg, '-m', 'LineWidth', 2), hold on, axis tight
        plot(upperbound, '-k')
        plot(lowerbound, '-k')
        axis tight
        ylabel(sprintf('T30 k %d\n(%d)',k, sum(Tparent==k)))
        
        axes(hax(1,2))
        bar(1:length(flyfoldersUnique), flyClusterCountDivZs30(:,k))
        hax(1,2).XTick = 1:length(flyfoldersUnique);
        hax(1,2).FontSize = 7;
        hax(1,2).XTickLabel = flyNumUnique;
        hax(1,2).XTickLabelRotation = 90;
        
        
        compNums = find(table_T10030(:,k));
        for i = 1 : sum(table_T10030(:,k)>0) %number of component clusters
            axes(hax(i+1,1)), hold on
            allR = R(Tcomponents==compNums(i),:);
            avg = mean(allR,1);
            upperbound = prctile(allR, 75, 1);
            lowerbound = prctile(allR, 25, 1);
            plot(avg, '-r', 'LineWidth', 2), hold on, axis tight
            plot(upperbound, '-k')
            plot(lowerbound, '-k')
            axis tight
            ylabel(sprintf('T60 k %d\n(%d)',compNums(i), sum(Tcomponents==compNums(i))))
            drawnow
            
            axes(hax(i+1,2))
            bar(1:NaF, flyClusterCountDivZs100(:,compNums(i)))
            hax(i+1,2).XTick = 1:NaF;
            hax(i+1,2).FontSize = 7;
            hax(i+1,2).XTickLabel = [];
            hax(i+1,2).XTickLabelRotation = 90;
        end
        % hax(i+1,2).XTickLabel = sortedflyMeta;
        
        export_fig(sprintf('clust30_%02d_R_ResponseAndDistribution.pdf',k))
        close
    end
    
    
end
