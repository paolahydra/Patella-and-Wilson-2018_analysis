
%% 
rootname = 'datalist_allWINDS';
rootfolder = 'WINDS_100';

ParentFolder = '/Users/galileo/Dropbox (HMS)/Data';
DataListFolder = fullfile(ParentFolder,rootfolder);
Folder2Save = fullfile(DataListFolder, 'Downsampled');
cd(Folder2Save)

load(fullfile(DataListFolder, [rootname '.mat']));     %DIR#-required list only
datalist = eval(rootname);

NaZ = length(datalist);
NaZ_100 = 23;%re-referenced
NaZ_86  = 44; 

for z = 1 : NaZ
    [sets.analysislist{z}, ~] = fileparts(datalist{z}); % runfolder address
    d = dir([sets.analysislist{z} '/*_analysis_downsampled.mat']);
    datalist{z} = fullfile(sets.analysislist{z}, d.name);
end
sets.datalist = datalist;

for z = 1 : NaZ
    aZ{z} = matfile(datalist{z});
end
[flyfolders, basenames] = extractRunfoldersBasenames(sets.datalist);
sets.runs_basename = extractRunsBasename_comboFlies( sets, basenames, basenames);
sets.flyfolder = Folder2Save;
sets.runs_basename = 'allWINDS';
sets.useSignificantPxls     = 2;    % 0: use all the pixels.
                                    % 2: each z gets its own.
 
%% datatable
dataTable = table( nan(length(datalist),1), ... % fly
    nan(length(datalist),1), ...                % run
    nan(length(datalist),1), ...                % set (100 or 86)
    nan(length(datalist),1), ...                % datalist #
    nan(length(datalist),1), ...                % alphanum #
    nan(length(datalist),1), ...                % stackZ
    nan(length(datalist),1), ...                % matrix z
    nan(length(datalist),1), ...                % R imaging side (ex ipsi)
    repmat({''}, length(datalist), 1), ...      % regist pix2keep
    repmat({''}, length(datalist), 1), ...      % boundaries for regist FOV mask
    repmat({''}, length(datalist), 1), ...      % regist avgF_map GREEN
    repmat({''}, length(datalist), 1), ...      % regist avgF_map RED, if any
    repmat({''}, length(datalist), 1), ...      % regist thresholded RED map, if any
    nan(length(datalist),1), ...                % zOffsetAlStack
    'VariableNames', {'fly', 'run', 'set', 'datalistRank', 'alphanumeralRank', 'stacksZ', 'matrixsZ', 'ipsilateral' ...
                       'reg_pixels2keep', 'boundaries', 'reg_avgGREEN', 'reg_avgRED', 'reg_maskRED', 'zOffsetAlStack' });
                  
% build up data table
for z = 1: length(datalist)
    dataTable.fly(z) = str2double(basenames{z}(4:6));
    dataTable.run(z) = str2double(basenames{z}(11:12));
end
%%
for z = 1:NaZ_100
    dataTable.set(z) = 100;
end
for z = NaZ_100+1:NaZ
    dataTable.set(z) = 86;
end
%%
for z = 1:NaZ
    dataTable.datalistRank(z) = z;
end

%%
save(fullfile(DataListFolder, [rootname '.mat']), 'dataTable', '-append')

%% registration 
targetParentFolder = fullfile(DataListFolder, 'anatomyImages');
saveName = fullfile(targetParentFolder, 'anatomy_alignment_metadata.mat');
regGeneralData = matfile(saveName); %single z planes to their own stack, general infos

% % update and save regGeneralData for once (HMS)
% regGeneralData.Properties.Writable = true;
% pattern_old = '/Users/galileo/Dropbox/';
% pattern_new = '/Users/galileo/Dropbox (HMS)/';
% flyfolders = regGeneralData.flyfolders;
% for i = 1 : length(flyfolders)
%     startIdx = strfind(flyfolders{i}, pattern_old); %this is usually 1, unless already resaved.
%     if startIdx == 1
%         flyfolders{i} = [pattern_new, flyfolders{i}(length(pattern_old)+1 : end)];
%     end
% end
% regGeneralData.flyfolders = flyfolders;
% flyfoldersUnique = regGeneralData.flyfoldersUnique;
% for i = 1 : length(flyfoldersUnique)
%     startIdx = strfind(flyfoldersUnique{i}, pattern_old); %this is usually 1, unless already resaved.
%     if startIdx == 1
%         flyfoldersUnique{i} = [pattern_new, flyfoldersUnique{i}(length(pattern_old)+1 : end)];
%     end
% end
% regGeneralData.flyfoldersUnique = flyfoldersUnique;
% regGeneralData.Properties.Writable = false;


flyfoldersUnique = regGeneralData.flyfoldersUnique;
% % %sort them numerically
% % flyfoldersUnique = sort(flyfoldersUnique);  % I don't want to do this.
% Prone to errors.

% flies2process = 1 : length(flyfoldersUnique);
% [flies2process,~] = listdlg('PromptString','Select flies to process:',...
%              'ListString',flyfoldersUnique, 'InitialValue', flies2process, 'ListSize', [300 340] );

%% set up dataTable (sorted by datalist)
dataTable = table( nan(length(datalist),1), ... % fly
    nan(length(datalist),1), ...                % run
    nan(length(datalist),1), ...                % datalist #
    nan(length(datalist),1), ...                % alphanum #
    nan(length(datalist),1), ...                % stackZ
    nan(length(datalist),1), ...                % matrix z
    nan(length(datalist),1), ...                % tdTom positive?
    nan(length(datalist),1), ...                % ipsi?
    nan(length(datalist),1), ...                % AMMC
    nan(length(datalist),1), ...                % tonotopy                      %TODO
    nan(length(datalist),1), ...                % steps                         %TODO
    nan(length(datalist),1), ...                % other kind of responses?      %TODO
    repmat({''}, length(datalist), 1), ...      % regist pix2keep
    repmat({''}, length(datalist), 1), ...      % boundaries for regist FOV mask
    repmat({''}, length(datalist), 1), ...      % regist avgF_map GREEN
    repmat({''}, length(datalist), 1), ...      % regist avgF_map RED, if any
    repmat({''}, length(datalist), 1), ...      % regist thresholded RED map, if any
    'VariableNames', {'fly', 'run', 'datalistRank', 'alphanumeralRank', 'stacksZ', 'matrixsZ', 'tdTomato', 'ipsilateral', ...
                      'AMMC', 'tonotopy', 'steps', 'otherResponses', ...
                       'reg_pixels2keep', 'boundaries', 'reg_avgGREEN', 'reg_avgRED', 'reg_maskRED' });
%                   
% %% build up data table
% for z = 1 : length(datalist)
%     dataTable.fly(z) = str2double(basenames{z}(4:6));
%     dataTable.run(z) = str2double(basenames{z}(11:12));
% end
% dataTable.datalistRank = (1:length(datalist))';

save(fullfile(DataListFolder, [rootname '.mat']), 'dataTable', '-append')    %ad interim safety save

%% split datalist by flies
cropRows = [0,1];   % cropRows(1) is #rows to crop from start, cropRows(2) is #rows to crop from end.
cropCols = [0,6];
warning off

for f = 1 : length(flyfoldersUnique)       % note: f goes by fly; 
                            %       z refers to DataListFolder (goes by runs)
    flyfolder = flyfoldersUnique{f};
    ZZ = strfind(datalist, flyfolder);
    zs_fly = zeros(size(ZZ));
    for i = 1:length(zs_fly)
        if ~isempty(ZZ{i})
            zs_fly(i) = 1; % indices into datalist relative to fly f
        end
    end
    zs_fly = logical(zs_fly);
    zetas = find(zs_fly);
    %%
    for iz = 1:sum(zs_fly)
        z = zetas(iz);
        disp(basenames{z})
        regData = matfile(fullfile(sets.analysislist{z}, ['stackRegistration_' basenames{z} '.mat']));
        stack2fly150_info = matfile(fullfile(flyfolder, ['stack2fly150sStack_' basenames{find(zs_fly,1)}(1:6) '.mat']));  
        
        %% reg_pixels2keep
        rotated_map = imrotate(aZ{z}.pixels2keep, regGeneralData.totAng(1,z)); %preliminary transformation
        rotated_map = double(rotated_map);
        recovered_map = imwarp(rotated_map,regData.tform,'OutputView',regData.Roriginal); %,'FillValues',[1,1,1]);
        al2fly150_map = imwarp(recovered_map,stack2fly150_info.tform,'OutputView',stack2fly150_info.Roriginal); %,'FillValues',[1,1,1]);
        dataTable.reg_pixels2keep(z) = {al2fly150_map};
        
        %% reg_avgGREEN (need to make it first)
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
        al2fly150_map = imwarp(recovered_map,stack2fly150_info.tform,'OutputView',stack2fly150_info.Roriginal); %,'FillValues',[1,1,1]);
        
        dataTable.reg_avgGREEN(z) = {al2fly150_map};
        
        if dataTable.tdTomato(z)
            %% reg_avgRED
            FileTif=fullfile(aZ{z}.runfolder, sprintf('AVG_regist_%sRED.tif', basenames{z}));
            TifLink = Tiff(FileTif, 'r');
            img=TifLink.read();
            TifLink.close();
            % improve
            img = mat2gray(img);
            img = adapthisteq(img);
            img = imadjust(img);
            % apply transformations
            rotated_map = imrotate(img, regGeneralData.totAng(1,z)); %preliminary transformation
            recovered_map = imwarp(rotated_map,regData.tform,'OutputView',regData.Roriginal); %,'FillValues',[1,1,1]);
            al2fly150_map = imwarp(recovered_map,stack2fly150_info.tform,'OutputView',stack2fly150_info.Roriginal); %,'FillValues',[1,1,1]);
            
            dataTable.reg_avgRED(z) = {al2fly150_map};
            %% reg_maskRED
            d = dir(fullfile(aZ{z}.runfolder, sprintf('*AVG_regist_thresholded_%sRED.tif', basenames{z})));
            FileTif=fullfile(aZ{z}.runfolder, d.name);
            TifLink = Tiff(FileTif, 'r');
            img=TifLink.read();
            TifLink.close();
            % apply transformations
            rotated_map = imrotate(img, regGeneralData.totAng(1,z)); %preliminary transformation
            recovered_map = imwarp(rotated_map,regData.tform,'OutputView',regData.Roriginal); %,'FillValues',[1,1,1]);
            al2fly150_map = imwarp(recovered_map,stack2fly150_info.tform,'OutputView',stack2fly150_info.Roriginal); %,'FillValues',[1,1,1]);
            
            dataTable.reg_maskRED(z) = {al2fly150_map};
        end
        
        %% boundaries
        mapmask = ones(size(aZ{z}.pixels2keep));
        % load specific run's transformation
        mapmask = imrotate(mapmask, regGeneralData.totAng(1,z));
        mapmask = imwarp(mapmask,regData.tform,'OutputView',regData.Roriginal);
        mapmask = imwarp(mapmask,stack2fly150_info.tform,'OutputView',stack2fly150_info.Roriginal);
        mapmask = logical(mapmask);
        [B,L] = bwboundaries(mapmask,'noholes');
        dataTable.boundaries(z) = {B{1}};
        
        %%
        if 0    % plot routine
            figure;
            imagesc(al2fly150_map);
            colormap(gray);
            axis image
            ax = gca;
            ax.Visible = 'on';
            ax.XTick = [];
            ax.YTick = [];
        end
    
    end
    
end
save(fullfile(DataListFolder, [rootname '.mat']), 'dataTable', '-append')       

%% look at reg_pixels2keep to select tonotopy+ runs 9both ipsi and contra)
close all
zetas = find(~dataTable.ipsilateral & dataTable.matrixsZ>9);
for iz = 1:length(zetas)
    z = zetas(iz);
    figure;
    imagesc(dataTable.reg_pixels2keep{z});
    colormap(gray);
    axis image
    ax = gca;
    ax.Visible = 'on';
    ax.XTick = [];
    ax.YTick = [];
    title(z)
end

%manually changed dataTable.tonotopy and resaved dataTable.
% included all possible runs that could show freq responses lateral wed and
% more dorsal (as much as we get), both ipsi and contra.
% eventually restrict to the runs I want to really use to show tonotopies

