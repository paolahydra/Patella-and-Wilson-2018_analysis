[imglist{1},folder] = uigetfile('/Users/galileo/Dropbox/Data/*RED*.tif','Select the RED stack .tiff file:');
[imglistGREEN{1},folder] = uigetfile([folder, '*.tif'],'Select the GREEN stack .tiff file:');

ws = 15;

%% load RED stack to be registered
FileTif=fullfile(folder,imglist{1});%'ImageStack.tif';
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
stack=zeros(nImage,mImage,NumberImages,'uint16');
tic
TifLink = Tiff(FileTif, 'r');
for i=1:NumberImages
    TifLink.setDirectory(i);
    stack(:,:,i)=TifLink.read();
end
TifLink.close();

stack = squeeze(stack);
nframesTot = size(stack,3);


%% load GREEN stack to be registered
FileTif=fullfile(folder,imglistGREEN{1});%'ImageStack.tif';
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
stackGREEN=zeros(nImage,mImage,NumberImages,'uint16');
tic
TifLink = Tiff(FileTif, 'r');
for i=1:NumberImages
    TifLink.setDirectory(i);
    stackGREEN(:,:,i)=TifLink.read();
end
TifLink.close();

stackGREEN = squeeze(stackGREEN);

%% window size == 15

for j = 0:ws-1;
    rem = mod (nframesTot-j, ws); %this is your remainder.
    stack_zproj = group_z_project(stack(:, :, j+1 : end-rem), ws); 
    if j == 0;
        target = stack_zproj(:,:,1);
    end
    [s(j+1).outs, s(j+1).stack] = stackRegisterMA(stack_zproj,target);
    writetiff(s(j+1).stack, fullfile(folder,sprintf('registered_%2d_%s', j, imglist{1})));
end
% output =  [error, diffphase, net_row_shift, net_col_shift]
% error     Translation invariant normalized RMS error between f and g
% diffphase     Global phase difference between the two images (should be
%               zero if images are non-negative).
% net_row_shift net_col_shift   Pixel shifts between images

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
shifts(:,1) = nanmean(errors)';
shifts(:,2) = nanmean(diffPhases)'; %not sure here
shifts(:,3) = nanmean(rowsShifts)';
shifts(:,4) = nanmean(colsShifts);
mean(shifts)
std(shifts)
save(fullfile(folder,'registrations_shifts.mat'), 'shifts')

% % Too much:
% mean:   0.8078    0.0000    0.3972    0.1744
% std:    0.0161    0.0004    0.4709    0.3855

%% apply shift to original RED and GREEN stacks
[~, stack2] = stackRegisterMA(stack,target, [], shifts);
writetiff(stack2, fullfile(folder,sprintf('registered_%s', imglist{1})));

[~, stackGR2] = stackRegisterMA(stackGREEN, [], [], shifts);
writetiff(stackGR2, fullfile(folder,sprintf('registered_%s', imglistGREEN{1})));



