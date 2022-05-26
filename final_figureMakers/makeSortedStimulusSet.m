function [stimTraces, ts_Stim_all, stN, clN, zeros_st_idxs] = makeSortedStimulusSet(T, sortedStimNs, sortedClassNs)
fs_acq = 4e4;
% ts_dec = T.ts_dec;
% % if classnum == 7 || classnum == 8
% %     iresp = ts_dec >= 0.25  &  ts_dec <= T.stimDur + 0.25;
% if classnum == 6
%     pipLatency = 1.25;
%     pipDur = 1;
%     iresp = ts_dec >= pipLatency+0.25  &  ts_dec <= pipLatency + pipDur + 0.25;
% else
%     iresp = ts_dec >= 0.25  &  ts_dec <= T.stimDur + 0.25;
% end
% ibasl = ts_dec <=0.15;
%% from trace_ChunkSubstim
stimTraces = [];
ts_Stim_all = [];           % x for f(x) =  stimTraces;
zeros_st_idxs = [];                         % changed. Now it tracks the onset of each stimulus
stN = [];
clN = [];
% divs_stimuli_idxs = 1;                      % previous zeros_st_idx, onset of a trial chunck (onset of baseline)
% divs_classes_idxs = 1;                      % what's this? 
for i = 1 : length(sortedClassNs)
    stCl = sortedClassNs(i);
    trialTypes = sortedStimNs{i};
    for n = trialTypes
        st_n = find(T(stCl).trialNums == n);
        sT = T(stCl).sensor(st_n,:);
        sT = sT(:);                         % now keeping it all
        if i == 1  &&  n == trialTypes(1)
            DCoffset = mean(sT(100:4e4)); % Baseline is usually longer than this.
        end
        stimTraces = cat(1,stimTraces,sT);  % full concatenated vector
        %% timestamps 1
        %         ts_stim = T(stCl).ts(iKeep(z).stClass(stCl).Response_ts);
        ts_stim = T(stCl).ts;  % keep it all. This starts negative, passes through zero at onset, and it's in seconds.
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
        stN_thisOne = n*ones(size(sT));
        stN = cat(1, stN, stN_thisOne);
        clN_thisOne = stCl*ones(size(sT));
        clN = cat(1, clN, clN_thisOne);
    end
end
stimTraces = (stimTraces-DCoffset) * 3; %um conversion %3 for new 30um piezo
end
