% reconstruct totAng as used during registration for flies 100, when I had
% a bug in anglesStck_single. And the bug should be that it had not been
% defined as:
% anglesStck_single = cat(2,allStacks.manualRotation); % bug fix
% therefore it should have only been zeros.

totAng_saved = regGeneralData.totAng % this I think was fixed later (before registering flies 86), and saved in its correct version
% totAng_saved =
%          0         0         0         0         0         0         0         0         0         0  -15.1200  -15.1200  -15.1200  -15.1200
%   -15.1200  -15.1200  -15.1200  -15.1200  -15.1200  -15.1200  -15.1200  -15.1200  -15.1200    1.0200    1.0200    1.0200    1.0200    1.0200
%     1.0200    1.0200    1.0200    1.0200         0         0         0         0         0         0         0         0         0         0
%          0         0         0         0         0         0         0         0         0         0         0         0         0         0
%          0         0         0         0         0         0         0         0         0         0         0

% totAng = %bugged
%       0         0         0         0         0         0         0         0         0         0  -15.1200  -15.1200  -15.1200  -15.1200
%   -15.1200  -15.1200  -15.1200  -15.1200  -15.1200  -15.1200  -15.1200  -15.1200  -15.1200  -10.0800  -10.0800  -10.0800  -10.0800  -10.0800
%   -10.0800  -10.0800  -10.0800  -10.0800         0         0         0         0         0         0         0         0         0         0
%          0         0         0         0         0         0         0         0         0         0         0         0         0         0
%          0         0         0         0         0         0         0         0         0         0         0

% totAng = %correct
%          0         0         0         0         0         0    2.0000    2.0000    2.0000    2.0000    1.8800    1.8800    1.8800    1.8800
%     1.8800    1.8800    1.8800    1.8800    1.8800    1.8800    1.8800    1.8800    1.8800    1.0200    1.0200    1.0200    1.0200    1.0200
%     1.0200    1.0200    1.0200    1.0200         0         0         0         0         0         0         0         0         0         0
%          0         0         0         0         0         0         0         0         0         0         0         0         0         0
%          0         0         0         0         0         0         0         0         0         0         0


% result: for flies 100, the reconstructed bugged version is actually the same as
% the saved one



%% reconstruct bugged version:
anglesStck_single = zeros(1, NaF);
acqAnglesStacks = zeros(1, NaF);

for f =1:NaF
    flyfolder = flyfoldersUnique{f};
    [~, flyname] = fileparts(flyfolder);
    flyname = flyname(1:6);
    ZZ = strfind(datalist, flyfolder);
    zs_fly = zeros(size(ZZ));
    for i = 1:length(zs_fly)
        if ~isempty(ZZ{i})
            zs_fly(i) = 1; % indices into datalist relative to fly f
        end
    end
    zs_fly = logical(zs_fly);
    
    %%  get each funcional run's acquisition angle
    for z = find(zs_fly)
        metaimage = aZ{z}.metaim;
        anglesRuns(z) = metaimage.scanRotation;
    end
    
    
    acqAnglesStacks(f) = allStacks(f).metaimage.scanRotation;
    
    anglesStck(zs_fly) = acqAnglesStacks(f) + anglesStck_single(f);  %this has been used incorrectly for flies 100.
    totAng = anglesRuns-anglesStck; %this will be updated for each fly
end



%% reconstruct unbugged version:
anglesStck_single = cat(2,allStacks.manualRotation); % bug fix
acqAnglesStacks = zeros(1, NaF);

for f =1:NaF
    flyfolder = flyfoldersUnique{f};
    [~, flyname] = fileparts(flyfolder);
    flyname = flyname(1:6);
    ZZ = strfind(datalist, flyfolder);
    zs_fly = zeros(size(ZZ));
    for i = 1:length(zs_fly)
        if ~isempty(ZZ{i})
            zs_fly(i) = 1; % indices into datalist relative to fly f
        end
    end
    zs_fly = logical(zs_fly);
    
    %%  get each funcional run's acquisition angle
    for z = find(zs_fly)
        metaimage = aZ{z}.metaim;
        anglesRuns(z) = metaimage.scanRotation;
    end
    
    
    acqAnglesStacks(f) = allStacks(f).metaimage.scanRotation;
    
    anglesStck(zs_fly) = acqAnglesStacks(f) + anglesStck_single(f);  %this has been used incorrectly for flies 100.
    totAng = anglesRuns-anglesStck; %this will be updated for each fly
end




%% check by aligning AVG functional images to all stacks
cropRows = [0,1];   % cropRows(1) is #rows to crop from start, cropRows(2) is #rows to crop from end.
cropCols = [0,6];
        
theseFlyNums = [166   169   175 ]; %  157   1701  1700   179   180]; 

for f = 1:length(theseFlyNums)
    flyNum = theseFlyNums(f);
    disp('--- fly: -------------------------------')
    disp(flyNum)
    disp('-')
    i_allSt = find(allStacks_flyNums == flyNum); %ef fix for fly 118
    if flyNum > 1000
        hemisph = mod(flyNum,10);
        flyNum = floor(flyNum/10);
        zs_fly = dataTable.fly==flyNum & dataTable.ipsilateral == hemisph; %ordinal numberwithin datalist
        runNames = dataTable.run(dataTable.fly==flyNum & dataTable.ipsilateral == hemisph); %actual run number, to reconstruct flyname
    else
        zs_fly = dataTable.fly==flyNum; %ordinal numberwithin datalist
        runNames = dataTable.run(dataTable.fly==flyNum);
    end  
    
    %% sort runs from ventral to dorsal
    zetas = find(zs_fly);
    sliceNums = regGeneralData.alignedSliceNumbers(1,zetas);
    [sliceNums,b] = sort(sliceNums);
    zetas = zetas(b);
   
    %% reg_avgGREEN (need to make it first)
    for iz = 1:sum(zs_fly)
        z = zetas(iz);          % change for flies 86
        sliceNum = sliceNums(iz);
        % load specific run's transformation
        regData = matfile(fullfile(runfolders{z}, ['stackRegistration_' basenames{z} '.mat']));
        disp(z)
        disp(aZ{z}.runfolder)
        
        FileTif=fullfile(aZ{z}.runfolder, sprintf('registered_%s.tif', basenames{z}));
        InfoImage=imfinfo(FileTif);
        mImage=InfoImage(1).Width;
        nImage=InfoImage(1).Height;
        NumberImages=length(InfoImage);
        tiffMaxi = zeros(nImage,mImage,1,NumberImages,'uint16');
        TifLink = Tiff(FileTif, 'r');
        for im=1:NumberImages
            TifLink.setDirectory(im);
            tiffMaxi(:,:,1,im)=TifLink.read();
        end
        TifLink.close();
        tiffMaxi = squeeze(tiffMaxi);
        img = mean(tiffMaxi, 3);
        % crop!!
        img(1:cropRows(1),:,:,:) = [];
        img(end-cropRows(2)+1:end,:,:,:) = [];
        img(:,1:cropCols(1),:,:) = [];
        img(:,end-cropCols(2)+1:end,:,:) = [];
        % improve
        img = mat2gray(img);
        img = adapthisteq(img);
        img = imadjust(img);
        % apply transformations
        rotated_map = imrotate(img, regGeneralData.totAng(1,z)); %preliminary transformation
        recovered_map = imwarp(rotated_map,regData.tform,'OutputView',regData.Roriginal); %,'FillValues',[1,1,1]);
        al2fly119_map = imwarp(recovered_map, allStacks(i_allSt).tform,'OutputView',allStacks(i_allSt).Roriginal); %,'FillValues',[1,1,1]);
        
        dataTable.reg_avgGREEN(z) = {al2fly119_map};
        savename = sprintf('/Users/galileo/Dropbox (HMS)/Data/WINDS_100/Downsampled/LinkedImages/alignedFunctionalMaps/alignfunc_fly%02d_ordZ%02d_%s.png', f, z, basenames{z});
        imwrite(al2fly119_map, savename)
    end
end
            

