%this version calculates shifts over smoothed tiff and applies also onto original unsmoothed tiff
[imglistTarget{1},folder] = uigetfile('/Users/galileo/Dropbox/Data/*.tif','Select the unsmoothed stack .tiff file:');
if strcmp(folder(end),'/')
    folder(end) = [];
end
[rootfolder, ~] = fileparts(folder);

imglistGREEN = dir([folder '/smoothed_' imglistTarget{1}]); %reference
imglistGREEN = imglistGREEN.name;


ws = 12;

%% load GREEN stacks to be registered
%target
FileTif=fullfile(folder,imglistTarget{1});%'ImageStack.tif';
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


%reference
FileTif=fullfile(folder,imglistGREEN);%'ImageStack.tif';
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


%%

nframesTot = size(stackGREEN,3);
for j = 0:ws-1;
    disp(j+1)
    rem = mod (nframesTot-j, ws); %this is your remainder.
    stack_zproj = group_z_project(stackGREEN(:, :, j+1 : end-rem), ws); 
    
    Gmag = zeros(size(stack_zproj));
    parfor jig = 1 : size(stack_zproj,3)        
%         [~, threshold] = edge(stack_zproj(:,:,jig), 'sobel');
%         fudgeFactor = .25;
%         BW = edge(stack_zproj(:,:,jig), 'sobel', threshold*fudgeFactor);
%         figure; imshow(BW)
        [Gmag(:,:,jig),~] = imgradient(stack_zproj(:,:,jig), 'sobel');
%         figure; imshowpair(Gmag, Gdir, 'montage');
    end

    if j == 0;
        ztg = floor(size(Gmag,3)/2)-4 : floor(size(Gmag,3)/2)+4;
        target = mean(Gmag(:,:,ztg),3);
%         target2 = Gmag(:,:,1);
%         figure; imshowpair(target2, target, 'montage');
    end
    
    [s(j+1).outs, s(j+1).stack] = stackRegisterMA(Gmag,target);
end


% for j = 0:ws-1;
%     writetiff(s(j+1).stack, fullfile(folder,sprintf('registered_%2d_%s', j, imglistGREEN)));
% end


% % output =  [error, diffphase, net_row_shift, net_col_shift]
% % error     Translation invariant normalized RMS error between f and g
% % diffphase     Global phase difference between the two images (should be
% %               zero if images are non-negative).
% % net_row_shift net_col_shift   Pixel shifts between images

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
shifts(:,4) = nanmean(colsShifts)';

save(fullfile(folder,'registrations_shifts.mat'), 'shifts')
    colsShifts(j+1, j+1 : end-rem) = cS(:);
    diffPhases(j+1, j+1 : end-rem) = dF(:);

%% apply shift to original RED and GREEN stacks
[~, stackUnSmoo] = stackRegisterMA(stack, [], [], shifts);
writetiff(stackUnSmoo, fullfile(folder,sprintf('registered_%s', imglistTarget{1})));


[~, stackGR2] = stackRegisterMA(stackGREEN, [], [], shifts);
writetiff(stackGR2, fullfile(folder,sprintf('registered_%s', imglistGREEN)));

% save zProject too
stack_zproj_stackGR2 = group_z_project(stackGR2); 
stack_zproj_stackGR2 = stack_zproj_stackGR2/max(stack_zproj_stackGR2(:));
stack_zproj_stackGR2 = adapthisteq(stack_zproj_stackGR2);
% figure; imagesc(stack_zproj_stackGR2); colormap(gray);
imwrite(uint16(stack_zproj_stackGR2*(2^16-1)),fullfile(rootfolder,sprintf('AVG_registered_%s', imglistGREEN)),'tif')

%% RED?
[~,nb] = fileparts(imglistGREEN);
d = dir(fullfile(folder, [nb 'RED.tif']));
if length(d) == 1
    % load RED stack to be registered
    FileTif=fullfile(folder,d.name);%'ImageStack.tif';
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
    % register and save
    [~, stackUnSmoo] = stackRegisterMA(stack,target, [], shifts);
    writetiff(stackUnSmoo, fullfile(folder,sprintf('registered_G2R_%s', d.name)));
    
    % save zProject too
    stack_zproj_stackR2 = group_z_project(stackUnSmoo);
    stack_zproj_stackR2 = stack_zproj_stackR2/max(stack_zproj_stackR2(:));
    stack_zproj_stackR2 = adapthisteq(stack_zproj_stackR2);
%     figure; imagesc(stack_zproj_stackR2); colormap(gray);
    imwrite(uint16(stack_zproj_stackR2*(2^16-1)),fullfile(rootfolder,sprintf('AVG_registered_G2R_%s', d.name)),'tif')
    
    % save composite
    composite(:,:,1) = uint8(stack_zproj_stackR2*255);
    composite(:,:,2) = uint8(stack_zproj_stackGR2*255);
    composite(:,:,3) = uint8(zeros(size(stack_zproj_stackGR2)));
%     figure; imagesc(composite); axis image; axis off; title(aZ.basename, 'Interpreter', 'none')
    imwrite(composite,fullfile(rootfolder,sprintf('AVG_composite_registered_G2R_%s', nb)),'tif')
end


%%
mean(shifts)
std(shifts)

figure('WindowStyle', 'normal', 'Position', [  50,  30, 1350, 720]); hold on
subplot(3,1,1); plot(shifts) 
title(imglistGREEN, 'Interpreter', 'none')
subplot(3,1,2); plot(shifts) 
title('distance')
plot(distances)
subplot(3,1,3); plot(shifts) 
title('frame to frame absolute distances')
plot(abs(diff(distances)))

name = [folder '/shifts_' imglistTarget{1}];
% save
export_fig(name, '-m3.5')

% % Too much:
% mean:   0.8078    0.0000    0.3972    0.1744
% std:    0.0161    0.0004    0.4709    0.3855

% % ?? looks good, given actual big shift. 
% mean:   0.5729   -0.0000    0.0362   -0.2232
% std:    0.0441    0.0001    0.0693    0.2658