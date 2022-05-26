load_winds_downsampled;
rootfolder = 'WINDS_100/anatomyImages';
ParentFolder = '/Users/galileo/Dropbox (HMS)/Data';
DataListFolder = fullfile(ParentFolder,rootfolder);
runfolders = sets.analysislist;
saveName = fullfile(DataListFolder, 'anatomy_alignment_metadata.mat');
regGeneralData = matfile(saveName);
datalist_address = '/Users/galileo/Dropbox (HMS)/Data/WINDS_100/datalist_allWINDS.mat';

folderTD = '/Users/galileo/Dropbox (HMS)/Data/TDTstacks/';
load('/Users/galileo/Dropbox (HMS)/Data/TDTstacks/anatomy_alignment_metadata.mat') %allStacks
allStacks_flyNums = cat(1,allStacks.flyNum);
warning off

%% matrix zetas
% by 5 um - OK
dataTable.zOffsetAlStack(dataTable.fly == 157) = 0;
dataTable.zOffsetAlStack(dataTable.fly == 179) = 4;
dataTable.zOffsetAlStack(dataTable.fly == 180) = 4;
for f = [157, 179, 180]
    dataTable.matrixsZ(dataTable.fly == f) = dataTable.stacksZ(dataTable.fly == f) + dataTable.zOffsetAlStack(dataTable.fly == f);
end


% by 10 um - OK
dataTable.zOffsetAlStack(dataTable.fly == 169) = 6;
dataTable.zOffsetAlStack(dataTable.fly == 170) = 1;
% separate offset for fly 170 L
dataTable.zOffsetAlStack(dataTable.fly == 170 & dataTable.ipsilateral == 0) = 5;
for f = [169, 170] %this is actually 170R only 
    % dataTable.stacksZ remains unchanged
    dataTable.matrixsZ(dataTable.fly == f) = dataTable.stacksZ(dataTable.fly == f)*2-1 + dataTable.zOffsetAlStack(dataTable.fly == f);
end


% flies 166 and 175 will get their own alignmen  (e.g. (runs 3-1-2:
% mZ:18-20-21) ) - OK
f = 166;
dataTable.matrixsZ(dataTable.fly == f) = [15 17 19 20 22 23 ];

f = 175;
dataTable.matrixsZ(dataTable.fly == f) = [20 19 18 17 16 15 14 20 19 21 22 17 15]; %some of these from the same nominal z have been distributed arbitrary to different matrixZ. Some runs could actually be better displayed with a flipped order.

save(datalist_address, 'dataTable', '-append')


%%
greenfolder = '/Users/galileo/Dropbox (HMS)/Data/WINDS_100/Downsampled/linkedIMAGES/alignedAnatomicalMaps_GREEN';
redfolder   = '/Users/galileo/Dropbox (HMS)/Data/WINDS_100/Downsampled/linkedIMAGES/alignedAnatomicalMaps_RED';
masksFolder = '/Users/galileo/Dropbox (HMS)/Data/WINDS_100/Downsampled/linkedIMAGES/alignedTDTMASKS';
tdtBorderFd = '/Users/galileo/Dropbox (HMS)/Data/WINDS_100/Downsampled/linkedIMAGES/alignedBorders';
% mkdir(greenfolder)
% mkdir(redfolder)
% mkdir(masksFolder)
% mkdir(tdtBorderFd)

%% new boundary making, plotting and saving into dataTable

%% make cluster maps - flies 100 only
theseFlyNums = [166   169   175 ]; %  157   1701  1700   179   180]; 

for f = 1:length(theseFlyNums)
    flyNum = theseFlyNums(f);
    i_allSt = find(allStacks_flyNums == flyNum); %ef fix for fly 118
    if flyNum > 1000
        hemisph = mod(flyNum,10);
        flyNum = floor(flyNum/10);
        zs_fly = dataTable.fly==flyNum & dataTable.ipsilateral == hemisph; %ordinal numberwithin datalist
%         runNames = dataTable.run(dataTable.fly==flyNum & dataTable.ipsilateral == hemisph); %actual run number, to reconstruct flyname
    else
        zs_fly = dataTable.fly==flyNum; %ordinal numberwithin datalist
%         runNames = dataTable.run(dataTable.fly==flyNum); %shouldn't these be sorted too?
    end  
    disp('--- fly: -------------------------------')
    disp(flyNum)
    disp('-')
    %% sort runs from ventral to dorsal
    zetas = find(zs_fly);
    sliceNums = regGeneralData.alignedSliceNumbers(1,zetas);
    [sliceNums,b] = sort(sliceNums);
    zetas = zetas(b);
        
    %% load transformation data and make maps
    for iz = 1:sum(zs_fly)
        z = zetas(iz);          % change for flies 86
        sliceNum = sliceNums(iz);
        Zaligned = dataTable.matrixsZ(z);
        Zwithin  = dataTable.stacksZ(z);
        % load specific run's transformation
        regData = matfile(fullfile(runfolders{z}, ['stackRegistration_' basenames{z} '.mat']));
        tagName = basenames{z};
        
        mapGreen = allStacks(i_allSt).alignedstackG(:,:, Zwithin);
        mapWhite = ones(size(mapGreen));
        figure; imshow(mapWhite); hold on
        
        mapmask = ones([99, 122]); %non-downsampled, cropped image size
        mapmask = imrotate(mapmask, regGeneralData.totAng(1,z)); %preliminary transformation
        mapmask = imwarp(mapmask,regData.tform,'OutputView',regData.Roriginal,'FillValues',0); %indexed image
        mapmask = imwarp(mapmask,allStacks(i_allSt).tform,'OutputView',allStacks(i_allSt).Roriginal,'FillValues',0); %ok
        mapmask(mapmask<0.7) = 0;
        mapmask = logical(mapmask);
        [B,L] = bwboundaries(mapmask,'noholes');
        boundary = B{1};
        dataTable.boundaries(z) = {boundary};
        plot(boundary(:,2), boundary(:,1), '-c')
        
        t = annotation('textbox');
        t.FontSize = 15;
        t.HorizontalAlignment = 'center';
        t.String = tagName;
        t.Interpreter = 'none';
        t.LineStyle = 'none';
        t.Position =  [0.4852    0.1389    0.0857    0.0320];
        t.Color = 'c';
        mapName = fullfile(tdtBorderFd, sprintf('alignALL_fly%03d_absZ%02d_%s_OwnRedBound&rectangleRun.eps',flyNum, Zaligned, tagName));
        export_fig(mapName)
        close
    end
end


%% make cluster maps - flies 86 only
theseFlyNums = [ 157   1701  1700   179   180]; 

for f = 1:length(theseFlyNums)
    flyNum = theseFlyNums(f);
    F = flyNum;
    if flyNum > 1000
        hemisph = mod(flyNum,10);
        flyNum = floor(flyNum/10);
        zs_fly = dataTable.fly==flyNum & dataTable.ipsilateral == hemisph;      %logical indexing within datalist
%         runNames = dataTable.run(dataTable.fly==flyNum & dataTable.ipsilateral == hemisph);     %actual run number, to reconstruct flyname
    else
        zs_fly = dataTable.fly==flyNum;     %logical indexing within datalist
%         runNames = dataTable.run(dataTable.fly==flyNum);
    end 
    i_allSt = find(allStacks_flyNums == flyNum); %ef fix for fly 118
    disp('--- fly: -------------------------------')
    disp(flyNum)
    disp('-')
    %% sort runs from ventral to dorsal
    zetas = find(zs_fly); 
    sliceNums = dataTable.stacksZ(zs_fly);   
    [sliceNums,b] = sort(sliceNums);        % stack slices paired to any functional run, (sorted ventral to dorsal)
    zetas = zetas(b);                       % ordinal numbers relative to datalist, (sorted ventral to dorsal)
    
    %% load transformation data and make maps
    for iz = 1:sum(zs_fly)          %iz is local counter within fly
        z = zetas(iz);              % relative to datalist
        z_86 = z - NaZ_100;         % now this is relative to the R_86 subset
        sliceNum = sliceNums(iz);   % still relative to its own stack
        Zaligned = dataTable.matrixsZ(z);
        Zwithin  = dataTable.stacksZ(z);
        % load specific run's transformation
        regData = matfile(fullfile(runfolders{z}, ['stackRegistration_' basenames{z} '.mat']));
        tagName = basenames{z};
        
        mapGreen = allStacks(i_allSt).alignedstackG(:,:, Zwithin);
        mapWhite = ones(size(mapGreen));
        figure; imshow(mapWhite); hold on
        
        mapmask = ones([85, 122]); %non-downsampled, cropped image size
        mapmask = imrotate(mapmask, regGeneralData.totAng(1,z)); %preliminary transformation
        mapmask = imwarp(mapmask,regData.tform,'OutputView',regData.Roriginal,'FillValues',0); %indexed image
        mapmask = imwarp(mapmask,allStacks(i_allSt).tform,'OutputView',allStacks(i_allSt).Roriginal,'FillValues',0); %ok
        mapmask(mapmask<0.7) = 0;
        mapmask = logical(mapmask);
        [B,L] = bwboundaries(mapmask,'noholes');
        boundary = B{1};
        dataTable.boundaries(z) = {boundary};
        plot(boundary(:,2), boundary(:,1), '-c')
        
        t = annotation('textbox');
        t.FontSize = 15;
        t.HorizontalAlignment = 'center';
        t.String = tagName;
        t.Interpreter = 'none';
        t.LineStyle = 'none';
        t.Position =  [0.4852    0.1389    0.0857    0.0320];
        t.Color = 'c';
        mapName = fullfile(tdtBorderFd, sprintf('alignALL_fly%03d_absZ%02d_%s_OwnRedBound&rectangleRun.eps',F, Zaligned, tagName));
        export_fig(mapName)
        close
    end
end








