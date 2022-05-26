function [registeredTarget, varargout] = stackRegistration( referenceStack, targetStack, folder, varargin )
%% new simple translation applied to target as of Sept 22, 2016 - PP
if nargin == 4
    ws = varargin{1};
else
    ws = 2;
end

if nargin == 5
    saveShifts = varargin{2};
else
    saveShifts = 1;
end

%%
nframesTot = size(referenceStack,3);
for j = 0:ws-1
    sprintf('\nIn progress: %3.1f percent', (j+1)/ws *100)
    rem = mod (nframesTot-j, ws); %this is your remainder.
    stack_zproj = group_z_project(referenceStack(:, :, j+1 : end-rem), ws); 
    stack_zproj = cast(stack_zproj, 'uint16');
    
    Gmag = zeros(size(stack_zproj));
    parfor z = 1 : size(stack_zproj,3)
        [Gmag(:,:,z),~]  = imgradient(imadjust(adapthisteq(stack_zproj(:,:,z))), 'sobel');
    end
    
    if j == 0
        ztg = floor(size(Gmag,3)/2)-4 : floor(size(Gmag,3)/2)+4;
        target = mean(Gmag(:,:,ztg),3);
    end
    
    [s(j+1).outs, s(j+1).stack] = stackRegisterMA_PP(Gmag,target);
end


%% interpolate x-y shifts, 
errors     = nan(ws,nframesTot);
diffPhases = nan(ws,nframesTot);
rowsShifts = nan(ws,nframesTot);
colsShifts = nan(ws,nframesTot);

for j = 0:ws-1
    rem = mod (nframesTot-j, ws); %this is your remainder.
    
    rS = repmat(s(j+1).outs(:,3)',ws,1);
    cS = repmat(s(j+1).outs(:,4)',ws,1);
    dF = repmat(s(j+1).outs(:,2)',ws,1);
    eR = repmat(s(j+1).outs(:,1)',ws,1);
    
    rowsShifts(j+1, j+1 : end-rem) = rS(:);
    colsShifts(j+1, j+1 : end-rem) = cS(:);
    diffPhases(j+1, j+1 : end-rem) = dF(:);
    errors(j+1, j+1 : end-rem) = eR(:);
end
if ws>1
    shifts(:,1) = nanmean(errors)';
%     shifts(:,2) = nanmean(diffPhases)'; %not sure here
    shifts(:,2) = circ_mean(diffPhases)'; 
    shifts(:,3) = nanmean(rowsShifts)';
    shifts(:,4) = nanmean(colsShifts)';
else
    shifts(:,1) = (errors)';
    shifts(:,2) = (diffPhases)'; 
    shifts(:,3) = (rowsShifts)';
    shifts(:,4) = (colsShifts)';
end

if saveShifts
    save(fullfile(folder,'registrations_shifts.mat'), 'shifts')
end

%% apply shift to original target stack
% [~, registeredTarget] = stackRegisterMA_PP(targetStack, [], [], shifts);
disp('Translating target stack...')
registeredTarget = zeros(size(targetStack));
colsShifts = shifts(:,4);
rowsShifts = shifts(:,3);
colsShifts( abs(colsShifts) < 0.001 ) = 0;
rowsShifts( abs(rowsShifts) < 0.001 ) = 0;
parfor pl = 1:size(shifts, 1)
    registeredTarget(:,:,pl) = imtranslate(targetStack(:,:,pl), [colsShifts(pl), rowsShifts(pl)]);
end
disp('Done.')
if nargout > 1
    varargout{1} = shifts;
end
