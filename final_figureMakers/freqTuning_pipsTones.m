function fTC = freqTuning_pipsTones(trialTypes, aZ, ROIs, windowsize, savename) 
% trialTypes is a cell with different subsets of stimuli
sortedClassNs = zeros(1, length(trialTypes));
for i = 1 : length(trialTypes)
    sortedClassNs(i) = aZ{1}.stimuli2Class(trialTypes{i}(1), 1);
end

for z = 1 : length(aZ)
    fprintf('%d ',z)
    classNum = sortedClassNs(1);
    T = aZ{z}.T(1,classNum);
    t1 = T.ts_dec;
    if z == 1, t1size = length(t1); end
    iresp = t1 >= 0  &  t1 <= T.stimDur + 0.5;
    ibasl = t1 >= -1.7145 & t1 <= 0.12;
    
    for i = 1 : length(trialTypes)
        if sortedClassNs(i) ~= classNum %update T
            classNum = sortedClassNs(i);
            T = aZ{z}.T(1,classNum);
            t1 = T.ts_dec;
            if z == 1, t1size = length(t1); end
            iresp = t1 >= 0  &  t1 <= T.stimDur + 0.5;
            ibasl = t1 >= -1.7145 & t1 <= 0.12;
        end
        stimNums = trialTypes{i};
        % let's not use extended_baselines for now
        
        t = squeeze(reshape(T.tiffSM, size(T.tiffSM,1)*size(T.tiffSM,2),1,size(T.tiffSM,3), size(T.tiffSM,4)));
        % t is now (pixels, allStimsThisClass*reps, timepoints)
%         if size(t,3) < t1size %sometimes baseline has one less point
%         %HERE WAS THE NAN PROBLEM
%             addbefore = t1size - size(t,3); 
%             t = cat(3, nan(size(t,1), size(t,2), addbefore), t);
%         end
        
        t_baseline = t( :,:,ibasl); %all trials this class
        t_response = t( :,:,iresp);
        % sort trials out
        for ist = 1 : length(stimNums)
            tT(i).baselines(:,z,:,ist) = squeeze(nanmean(t_baseline(:, T.trials==stimNums(ist), :), 2)); % (all5160pixels, timepoints), trial averaged per st type
            tT(i).responses(:,z,:,ist) = squeeze(nanmean(t_response(:, T.trials==stimNums(ist), :), 2)); % (all5160pixels, timepoints), trial averaged per st type
        end
    end
end

% %% debug nans
% z = 5;
% classNum = sortedClassNs(1);
% T = aZ{z}.T(1,classNum);
% for ist = 1 : 18*9
%     tiff1 = squeeze(T.tiffSM(:,:,ist,:));
%     countNans(ist) = sum(isnan(tiff1(:)));
% end
% sum(countNans)
% % % MIJ.createImage(tiff1)
% % RESULT: no NANs here

%% perhaps save tT now, as it is long to compute
save(savename, 'tT')

%% divide by cluster-ROIs and average pixels within, calculate and store statistics
clear fTC
for i = 1 : length(trialTypes)
    baselines = squeeze(reshape(tT(i).baselines, [], 1, size(tT(i).baselines,3), size(tT(i).baselines,4) ) ) ;
    responses = squeeze(reshape(tT(i).responses, [], 1, size(tT(i).responses,3), size(tT(i).responses,4) ) ) ;
    for k = 1 : size(ROIs,2)
        ROI = logical(reshape(squeeze(ROIs(:,k,:)), [], 1));
        fTC(i).baselines(:,:,k) = squeeze(nanmean(baselines(ROI,:,:)));
        fTC(i).responses(:,:,k) = squeeze(nanmean(responses(ROI,:,:)));
%         fTC(i).stdNoise(:,:,k) = std(fTC(i).baselines(:,:,k));
%         fTC(i).meanNoise(:,:,k) = nanmean((fTC(i).baselines(:,:,k))); % NB it was just f, not df
% pool across stimtypes
        fTC(i).stdNoise(:,:,k) = std(reshape(fTC(i).baselines(:,:,k), [],1,1) );
        fTC(i).meanNoise(:,:,k) = nanmean(reshape(fTC(i).baselines(:,:,k), [],1,1 ) ); % NB it was just f, not df
        
        fTC(i).movAvgResponse(:,:,k) = movmean(fTC(i).responses(:,:,k), windowsize);
        fTC(i).normPeak(:,:,k) = (max(fTC(i).movAvgResponse(:,:,k)) - fTC(i).meanNoise(:,:,k)) ./ fTC(i).stdNoise(:,:,k);
    end
end

%% check output
[hfig, hax] = figureTracesI_PP_pixels( 11, 200, 30 );
hax = flipud(hax); % 1 is bottom, as indendrogram
for Ki = 1:11
    axes(hax(Ki)), hold on
    for i = 1:3
        plot(fTC(i).normPeak(:,:,Ki))
    end
end

% %% check cluster 4
% MIJ.start;
% 
% i = 3;
% baselines =tT(i).baselines;
% responses = tT(i).responses;
% size(baselines)
% 
% k = 4;
% ROI = logical(ROIs(:,k,:));
% 
% z = 1; 
% BasZ = reshape(baselines(:,z,:,1), 60, 86, []); size(BasZ)
% ResZ = reshape(responses(:,z,:,1), 60, 86, []); size(ResZ)
% MIJ.createImage(cat(3, BasZ, ResZ))
% 
% % only show clusterin image
% avgBasZ = nanmean(baselines(:,z,:,:),4); size(avgBasZ)
% avgBasZ(~ROI(:,:,z),:,:) = 0;
% avgBasZ=reshape(avgBasZ, 60, 86, [], 1);
% MIJ.createImage(avgBasZ)
% 
% % show full image
% avgBasZ = nanmean(baselines(:,z,:,:),4); size(avgBasZ)
% avgBasZ=reshape(avgBasZ, 60, 86, [], 1);
% MIJ.createImage(avgBasZ)
% 
% 
% %% debug cluster 4
% k = 4;
% basel_allStims = mean(tT(i).baselines, 4);
% size(basel_allStims)    %   5160          79          22
% size(ROIs)              %   5160          11          79
% ROI4byz = logical(squeeze(ROIs(:,k,:)));
% size(ROI4byz)           %   5160          79
% 
% for z = 1:79
%     bas(:,z) = squeeze(mean(basel_allStims(ROI4byz(:,z),z,:),1)); %only timepoints by z plane, averaged across cluster 4 pixels and all stimuli
% end
% % RESULT: NaNs all over in almost all zplanes, sometimes in all timepoints,
% % more often only in ~fisrt half of timepoints.



