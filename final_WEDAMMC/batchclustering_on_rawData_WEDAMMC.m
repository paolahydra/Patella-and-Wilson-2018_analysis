clusterFolder = '/Users/galileo/Dropbox (HMS)/Data/WEDAMMC_piezo/Downsampled/rawDataLinkage/ward_ZSCpix';
treeFileName = 'output_ward_Rzsp_linkage_wedammc.mat';
load(fullfile(clusterFolder, treeFileName), 'Zc');

for K = 15:30
    disp(K)
    clearvars -except R Rzs K pixel2keep iKeep pxKeep total_within Zc
    maxClust = K;
    % maxClustComponents = 150;
    maxLeavesDend = maxClust;
    distanceMetrics = 'euclidean';
    linkageMethod = 'ward';
    % temporary dendr
    figure; [~,Tdend,outpermLeaves] = dendrogram(Zc, maxLeavesDend, 'Orientation', 'left');
    axDend = gca;
    axDend.YTickLabel = 1:maxLeavesDend;
    export_fig(sprintf('dendrogram_maxLeaves%d.pdf',maxLeavesDend))
    close
    %% cluster (main), sort klusts based on dendrogram -- AND THEN RE-SORT again [DO]
    Tparent = cluster(Zc,'maxclust',maxClust);
    
    % GOAL: Ki order of superclusts is fixed. Starting from the left-most
    % cluster-branch in dendrogram, we need to figure out which klust(Ki) it
    % corresponds to.
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
    save(sprintf('klust_maxClust%d.mat', maxClust), 'klust', 'dendrOrder', 'Tparent');
    nK = length(unique(Tparent));
    
    %% total within distance
%     for k = 1:nK
%         disp(k)
%         withinDistances = bsxfun(@minus, Rzs(Tparent==k,:),  mean(Rzs(Tparent==k,:)) );
%         withinDistances = withinDistances.^2;
%         withinDistances = sum(withinDistances, 2);
%         withinDistances = withinDistances.^(1/2);
%         wDs(k) = sum(withinDistances);
%         %     pw = pdist(Rzs(Tparent==k,:));
%         %     wr_total(k) = sum(pw);
%     end
%     total_within(nK) = sum(wDs)
%     save(sprintf('klust_maxClust%d.mat', maxClust), 'wDs', '-append');
%     
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
    export_fig(sprintf('clustMain_%d_mtraces.pdf', nK))
    close
end
