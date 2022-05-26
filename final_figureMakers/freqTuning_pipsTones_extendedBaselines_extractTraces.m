function freqTuning_pipsTones_extendedBaselines_extractTraces(trialTypes, aZ, savename) 
% trialTypes is a cell with different subsets of stimuli
sortedClassNs = zeros(1, length(trialTypes));
for i = 1 : length(trialTypes)
    sortedClassNs(i) = aZ{1}.stimuli2Class(trialTypes{i}(1), 1);
end


savenameGen = savename(1:end-4);
for i = 1 : length(trialTypes)
    mClass(i).matfile = matfile(sprintf('%s%d.mat',savenameGen, i), 'Writable',true);
end

for z = 1 : length(aZ)
    fprintf('%d ',z)
    classNum = sortedClassNs(1);
    BX = aZ{z}.baselinesExtended;
    BX = squeeze(reshape(BX, size(BX,1) * size(BX,2), 1,size(BX,3), size(BX,4)));
    T = aZ{z}.T(1,classNum);
    t1 = T.ts_dec;
    if z == 1, t1size = length(t1); end
    iresp = t1 >= 0  &  t1 <= T.stimDur + 0.5;
    ibasl = t1 >= -1.7145 & t1 <= 0.12;
    
    for i = 1 : length(trialTypes)
        m = mClass(i).matfile;
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
        baseline = [];
        response = [];
        for ist = 1 : length(stimNums)
            trialAvg_baseline = squeeze(nanmean(t_baseline(:, T.trials==stimNums(ist), :), 2)); % (all5160pixels, timepoints), trial averaged per st type
            basel_this = cat(2, trialAvg_baseline, BX(:,1:80,stimNums(ist)) );
            basel_this = permute(basel_this, [1,3,2]);
            baseline = cat(4, baseline, basel_this);
            
            resp_this = squeeze(nanmean(t_response(:, T.trials==stimNums(ist), :), 2));
            resp_this = permute(resp_this, [1,3,2]);
            response = cat(4, response, resp_this);
        end
        %write
        if z == 1
            m.baselines = baseline; %remove 4 extra points to level every run up without nans
            m.responses = response; % (all5160pixels, timepoints), trial averaged per st type
        else
            size(response)
            size(baseline)
            if i <= 3
                if size(response, 3) > 18 %this could be either 18 or 29
                    response(:,:,19:end,:) = [];
                elseif size(response,3) < 18
                    response = cat(3, response(:, :, 1:18-size(response,3), :), response);
                end
            else
                if size(response, 3) > 29 %this could be either 18 or 29
                    response(:,:,30:end,:) = [];
                elseif size(response,3) < 29
                    response = cat(3, response(:, :, 1:29-size(response,3), :), response);
                end
            end
            if size(baseline, 3) > 102
                baseline(:,:,102:end,:) = [];
            elseif size(baseline,3) < 102
                baseline = cat(3, baseline(:, :, 1:102-size(baseline,3), :), baseline);
            end
            
            m.baselines(:,z,:,:) = baseline; %remove 4 extra points to level every run up without nans
            m.responses(:,z,:,:) = response; % (all5160pixels, timepoints), trial averaged per st type
        end
    end
end

end


