function cmap = makeRainbowCmap(varargin)
if nargin
    nK = varargin{1};
end


map1 = brewermap(6, 'PuRd');
map1 = map1(3:end,:);
map2 = flipud(brewermap(5, 'YlOrRd')); %YlOrRd > OrRd if I don't need yellow elsewhere
map2 = map2(1:4,:);
% map3 = brewermap(5, 'YlGn'); %Greens
% map3 = map3(1:end-1,:);
map3 = brewermap(11, 'YlGn'); %Greens
map9 = map3(2:end,:);

map10(1,:) = map9(1,:).*1/2 + map9(2,:).*1/2;
map10(2,:) = map9(3,:).*1/2 + map9(4,:).*1/2;
map10(3,:) = map9(5,:).*1/2 + map9(6,:).*1/2;
map10(4,:) = map9(7,:).*1/2 + map9(8,:).*1/2;
map10(5,:) = map9(9,:).*1/2 + map9(10,:).*1/2;
map3 = map10(2:end,:);
% map3 = map9(3:2:9,:);


% % % map4 = flipud(brewermap(10, 'BrBG'));
% % % map4 = map4(2:4,:);
map4 = [0    0.5172    0.5862; 0    1.0000    0.7586]; %dark aqua, acqua from distinguish_col(11,10)
map4 = [0.1039    0.5947    0.5335;  0.2510    0.9020    0.7578];
% map4(1,:) = [0    0.2706    0.1608];        % [0    0.2667    0.1059];

map5 = [0.41, 0.41, 0.41; 0.12 0.12 0.12];

blue =  [0.1    0.45  0.95]                    %[0.0       0.4       0.9];
brown = brewermap(7, 'BrBG'); brown = brown(1,:);           %[0.5490    0.3176    0.0392]
yellow = brewermap(7, 'Set2'); yellow = yellow(6,:);        %[1.0000    0.8510    0.1843]
% violet = brewermap(6, 'Purples'); violet = violet(end,:);   %[0.3294    0.1529    0.5608]
violet = [0.45    0.14    0.5709];

cmap = cat(1, blue, violet, map1, map2, map3, map4, map5 ,brown);

figure; imagesc(1:size(map1,1)); colormap(map1); title('map1: PuRd')
figure; imagesc(1:size(map2,1)); colormap(map2); title('map2: YlOrRd')
figure; imagesc(1:size(map3,1)); colormap(map3); title('map3: YlGn')
figure; imagesc(1:size(map4,1)); colormap(map4); title('map4')
figure; imagesc(1:size(map5,1)); colormap(map5); title('map5')

%% re map 2 - linear RED-YELLOW 
map2r(:,1) = [0.7400    0.94    0.99    1.0000];  %linspace(0.8,1,4)';
map2r(:,2) = linspace(0,0.8,4)';
map2r(:,3) = linspace(0,0.2,4)';
figure; imagesc(1:size(map2r,1)); colormap(map2r); title('map2r: YlOrRd new')

%% re map 3 - YELLOW-GREEN
map3r(:,1) = linspace(0.71, 0,      4)';
map3r(:,2) = linspace(0.87,    0.3267, 4)';
map3r(:,3) = linspace(0.32,  0.1, 4)';
map3r(3:4,:) = map3(3:4,:);
figure; imagesc(1:size(map3r,1)); colormap(map3r); title('map3r: greens')

%% re map mapblues 
clear mapblues
blues = flipud(brewermap(9, 'GnBu'));
blues(6:end,:) = [];
aquadarker = [0.6,0.4]*map4;
mapblues(1,:) = aquadarker;
% mapblues(2,:) = blues(4,:);
mapblues(2,:) = blue; %mean(blues(2:4,:));
mapblues(3,:) = blues(1,:);
mapblues(3,3) = 0.6;
figure; imagesc(1:size(mapblues,1)); colormap(mapblues); title('mapblues')
%could average the middle two if you have an extra colro


blues2 = [aquadarker; ...
          0.150     0.500       1;...
          0.0714    0.2510      0.6500;...
          0         0           0.48]; %new version brown-less
figure; imagesc(1:size(blues2,1)); colormap(blues2); title('blues2')

blues3 = [aquadarker; ...
         0.2157    0.7569    0.9451; ...
         0    0.4471    0.7373; ...
         0.2157    0.1843    0.5529];
figure; imagesc(1:size(blues2,1)); colormap(blues2); title('blues2')
figure; imagesc(1:size(blues3,1)); colormap(blues3); title('blues3')

%% re mapviolets
clear mapviolets
mapviolets(1,:) = violet;
mapviolets(2:5,:) = map1;
figure; imagesc(1:size(mapviolets,1)); colormap(mapviolets); title('mapviolets')


%% combine
% cmap = cat(1, map5, map2r, map3r, mapblues, mapviolets, brown); %old
% version with brown

cmap = cat(1, map5, map2r, map3r, blues3, mapviolets);
figure; imagesc(1:size(cmap,1)); colormap(cmap); title('final cmap')
save('/Users/galileo/Dropbox (HMS)/Data/cmap19_spectral.mat', 'cmap', '-append')

%% make a 11- version
keep12 = [1 2 3 5 6 8 10 12 13 14 15 17];
cmap12 = cmap(keep12,:);
figure; imagesc(1:size(cmap12,1)); colormap(cmap12); title('final cmap12')
save('/Users/galileo/Dropbox (HMS)/Data/cmap19_spectral.mat', 'cmap12', '-append')

keep11 = [1 2 3 5 6 8 10 12 14 15 17];
cmap11 = cmap(keep11,:);
figure; imagesc(1:size(cmap11,1)); colormap(cmap11); title('final cmap11')
save('/Users/galileo/Dropbox (HMS)/Data/cmap19_spectral.mat', 'cmap11', '-append')


end

