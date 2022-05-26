function fTC = freqTuning_pipsTones_extendedBaselines(trialTypes, ROIs, movingwindowsize, savename, responseWindowPoints) 

savenameGen = savename(1:end-4);
for i = 1 : length(trialTypes)
    mClass(i).matfile = matfile(sprintf('%s%d.mat',savenameGen, i));
end

%% divide by cluster-ROIs and average pixels within, calculate and store statistics
clear fTC
for i = 1 : length(trialTypes)
    disp(i)
    baselines = mClass(i).matfile.baselines;
    responses = mClass(i).matfile.responses;
%     if ~isempty(responseWindowPoints)
%         responses = responses(:,:,1:responseWindowPoints,:);
    responses = responses(:,:,3:7,:); %4:9 for 02 07 sec
%     end
    baselines = squeeze(reshape(baselines, [], 1, size(baselines,3), size(baselines,4) ) ) ;
    responses = squeeze(reshape(responses, [], 1, size(responses,3), size(responses,4) ) ) ;
    size(baselines)
    size(responses)
    for k = 1 : size(ROIs,2)
        ROI = logical(reshape(squeeze(ROIs(:,k,:)), [], 1));
        fTC(i).baselines(:,:,k) = squeeze(nanmean(baselines(ROI,:,:)));
        fTC(i).responses(:,:,k) = squeeze(nanmean(responses(ROI,:,:)));
%         fTC(i).stdNoise(:,:,k) = std(fTC(i).baselines(:,:,k));
%         fTC(i).meanNoise(:,:,k) = nanmean((fTC(i).baselines(:,:,k))); % NB it was just f, not df
% pool across stimtypes
        fTC(i).stdNoise(:,:,k) = std(reshape(fTC(i).baselines(:,:,k), [],1,1) );
        fTC(i).meanNoise(:,:,k) = nanmean(reshape(fTC(i).baselines(:,:,k), [],1,1 ) ); % NB it was just f, not df
        
        fTC(i).movAvgResponse(:,:,k) = movmean(fTC(i).responses(:,:,k), movingwindowsize);
        fTC(i).normPeak(:,:,k) = (max(fTC(i).movAvgResponse(:,:,k)) - fTC(i).meanNoise(:,:,k)) ./ fTC(i).stdNoise(:,:,k);
    end
end

disp('not saving fTC to disk')
% save('fTC_data_01_06secAfterOnset.mat', 'fTC');


