rootname = 'datalist_allWINDS';
rootfolder = 'WINDS_100';

ParentFolder = '/Users/galileo/Dropbox (HMS)/Data';
DataListFolder = fullfile(ParentFolder,rootfolder);
Folder2Save = fullfile(DataListFolder, 'Downsampled');
cd(Folder2Save)

load(fullfile(DataListFolder, [rootname '.mat']));     %DIR#-required list only
datalist = eval(rootname);
clear(rootname)
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
set(0,'DefaultFigureWindowStyle','docked')


flyfoldersUnique = unique(flyfolders, 'stable');
NaF = length(flyfoldersUnique);
for i=1:NaF, flyNumUnique(i) = str2double(flyfoldersUnique{i}(end-5:end-3)); end    % unsorted by FMatrix, but sortd by flyfoldersUnique!!!


%load R
load(fullfile(Folder2Save, 'R_matrix_downsampled_winds_100&86.mat')); 

clear d rootname rootfolder i z