rootname = 'datalist_allWINDS';
rootfolder = 'WINDS_100';

ParentFolder = '/Users/galileo/Dropbox (HMS)/Data';
DataListFolder = fullfile(ParentFolder,rootfolder);

load(fullfile(DataListFolder, [rootname '.mat']));     %DIR#-required list only
datalist = eval(rootname);
clear(rootname)
clear(DataListFolder)
for z = 1 : length(datalist)
    [sets.analysislist{z}, ~] = fileparts(datalist{z}); % runfolder address
    d = dir([sets.analysislist{z} '/*_analysis_downsampled.mat']);
    datalist{z} = fullfile(sets.analysislist{z}, d.name);
end


Folder2Save = fullfile(DataListFolder, 'Downsampled/SubsetFlies');
if exist(Folder2Save, 'dir') ~= 7, mkdir(Folder2Save); end
cd(Folder2Save)

%% remake datalist
% flies to keep
flyNumKeep = [175 179 180 157]; %NOTE - this is already sorted 100 then 86, for simplicity

zetas100 = ismember(dataTable.fly, flyNumKeep) & dataTable.set == 100;
zetas86 = ismember(dataTable.fly, flyNumKeep) & dataTable.set == 86;
zetas = zetas100 | zetas86;

NaZ_100 = sum(zetas100);
NaZ_86 = sum(zetas86);
NaZ = sum(zetas);

datalist = cat(2, datalist(zetas100), datalist(zetas86));
dataTable = dataTable(zetas,:);
clear zetas zetas100 zetas86 

save(fullfile(Folder2Save, 'datalist_selectWindFlies.mat'), 'datalist', 'dataTable', 'NaZ', 'NaZ_100', 'NaZ_86' )

%%
zetas100 =  dataTable.set == 100;
zetas86 =  dataTable.set == 86;

for z = 1 : NaZ
    [sets.analysislist{z}, ~] = fileparts(datalist{z}); % runfolder address
end
sets.datalist = datalist;

for z = 1 : NaZ
    aZ{z} = matfile(datalist{z});
end
[flyfolders, basenames] = extractRunfoldersBasenames(sets.datalist);
% sets.runs_basename = extractRunsBasename_comboFlies( sets, basenames, basenames);
sets.flyfolder = Folder2Save;
sets.runs_basename = 'allWINDS';
sets.useSignificantPxls     = 2;    % 0: use all the pixels.
                                    % 2: each z gets its own.
set(0,'DefaultFigureWindowStyle','docked')

%%
% [pixel2keep_100, iKeep_100, R_100] = makeSaveR_mf( aZ(1:NaZ_100), sets );
% [pixel2keep_86, iKeep_86, R_86] = makeSaveR_mf( aZ(NaZ_100+1:end), sets );
[pixel2keep_100, iKeep_100, R_100] = makeSaveR_mf( aZ(1:NaZ_100), sets, -2, 1.3 );
[pixel2keep_86, iKeep_86, R_86] = makeSaveR_mf( aZ(NaZ_100+1:end), sets, -2, 1.3 );
%%
% % if need to load
% R86load = load(fullfile(sets.flyfolder, 'R_matrix_unfiltered_Downsampled_100.mat'));
% R_86            = R86load.R;
% iKeep_86        = R86load.iKeep;
% pixel2keep_86   = R86load.pixel2keep;
% clear R86load

% build up
R = cat(1, R_100, R_86);
clear R_100 R_86

z = 1;
R_100.sizeMAP = size(aZ{z}.pixels2keep);
R_100.iKeep = iKeep_100;
R_100.basenames = basenames(zetas100);
R_100.pixel2keep = pixel2keep_100;
% R_100.R = R(1:sum(R_100.pixel2keep(:)),:);

z = NaZ_100+1;
R_86.sizeMAP = size(aZ{z}.pixels2keep);
R_86.iKeep = iKeep_86;
R_86.basenames = basenames(zetas86);
R_86.pixel2keep = pixel2keep_86;
% R_86.R = R(sum(R_100.pixel2keep(:))+1:end,:);

% save
filenameR = fullfile(Folder2Save, 'R_matrix_downsampled_unfiltered_winds_100&86.mat'); %only changed the path to selct flies folder, not filename.
save(filenameR, 'R_100', 'R_86', 'R', '-v7.3')

clear iKeep_86 iKeep_100 pixel2keep_86 pixel2keep_100

%% Save zs2 Rmatrix for linkage calculation, over flipped, sustained phase only
% makeFlippedR_winds;
makeFlippedR_winds_withbaseline_unfilt; %work also when there is no baseline window, except for a filename called in

%% compute linkage on orchestra ...
% do it

%% save extra stuff once and for all
save(fullfile(Folder2Save, 'datalist_selectWindFlies.mat'), 'basenames', 'flyfolders', 'sets', '-append')


