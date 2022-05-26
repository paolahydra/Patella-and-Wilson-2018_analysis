%% clean stimuli (command)
load_panSpec_downsampled;
load('/Users/galileo/Dropbox (HMS)/Data/chunks_stimfinal.mat', 'chunks')
load('/Users/galileo/Dropbox (HMS)/Data_Raw_Metadata_BackedUp/BackedUp_Metadata/fly171_run01_metadata_2016-12-20_163728.mat', 'stimuli')

%load sample aZ
aZ = matfile(datalist_allJON{1});
T = aZ.T;
fastStimMeta = aZ.fastStimMeta;
% Allstimuli = aZ.ALLstimuli;

folder2save = '/Users/galileo/Dropbox (HMS)/figures/stimuli';
cd(folder2save)


%% chirps  %deal with amplitude later
cn = 'chirps_down_up';
chunks.(cn)

%temporarily change ramp_correct to 0.5 0.5
edit Chirp_down
edit Chirp_up

for c = 1:length(chunks.(cn).classnames)
    clName = chunks.(cn).classnames{c};
    clN = fastStimMeta.classNum(find(strcmp(cat(1,fastStimMeta.className), clName),1));
    stN = chunks.(cn).(clName);
    
    ts = T(clN).ts;
    stim = stimuli(stN).stim;
%     stim.startFrequency = 150;
%     stim.endFrequency = 150;
    stim = stim.stimulus;
    
    %command prototype
    figure; plot(ts, stim);
    axis tight
    set(gca, 'ylim', [-1 1])
    set(gca, 'TickDir', 'out')
    export_fig(sprintf('%s_prototype.eps',clName))
    
    %sensor with amplitude
    figure; hold on; title('sensor')
    DC_offset = mean(T(clN).sensor(1,1e2:1:.5e4));
    sensor_um = (T(clN).sensor(1,:) - DC_offset) * 3;
    plot(ts, sensor_um);
    set(gca, 'TickDir', 'out')
    set(gca, 'xlim', [ts(1), ts(end)])
    ylabel('um')
    export_fig(sprintf('%s_sensor.eps',clName))
    
    %add spectrogram
    figure
    spectrogram(stim,10000,500,10000,4e4,'yaxis');    % used to be 20000 - 500 - 20000. Trying to distort less the temporal dimension.
    colorbar off
    ax = gca;
    ax.TickDir = 'out';
    ax.YLim = [0,0.6];
    title('spectrogram over command')
    export_fig(sprintf('%s_spectrogram.eps',clName))
end

close all
clc
chunks

%% songs  %deal with amplitude later
cn = 'courtships';
chunks.(cn)

for c = 1:length(chunks.(cn).classnames)
    clName = chunks.(cn).classnames{c};
    clN = fastStimMeta.classNum(find(strcmp(cat(1,fastStimMeta.className), clName),1));
    stN = chunks.(cn).(clName);
    
    stim = stimuli(stN).stim.stimulus;
    DC_offset = mean(T(clN).sensor(1,1e2:1:.5e4));
    sensor_um = (T(clN).sensor(1,:) - DC_offset) * 3;
    
    ts = 0:1/4e4:17.65; 
    ts(length(sensor_um)+1:end) = [];
    ts = ts-1.8;
    
    %command prototype
    figure; plot(ts, stim);
    set(gca, 'xlim', [ts(1), ts(end)])
    set(gca, 'TickDir', 'out')
    export_fig(sprintf('%s_prototype.eps',clName))
    
    %sensor with amplitude
    figure; hold on; title('sensor')
    DC_offset = mean(T(clN).sensor(1,1e2:1:.5e4));
    sensor_um = (T(clN).sensor(1,:) - DC_offset) * 3;
    plot(ts, sensor_um);
    set(gca, 'TickDir', 'out')
    set(gca, 'xlim', [ts(1), ts(end)])
    ylabel('um')
    export_fig(sprintf('%s_sensor.eps',clName))
    
    %add spectrogram
    figure
    spectrogram(stim,10000,500,10000,4e4,'yaxis');    % used to be 20000 - 500 - 20000. Trying to distort less the temporal dimension.
    colorbar off
    ax = gca;
    ax.TickDir = 'out';
    ax.YLim = [0,0.6];
    title('spectrogram over command')
    export_fig(sprintf('%s_spectrogram.eps',clName))
    
%     figure
%     spectrogram(sensor_um,10000,500,10000,4e4,'yaxis');    % used to be 20000 - 500 - 20000. Trying to distort less the temporal dimension.
%     colorbar off
%     ax = gca;
%     ax.TickDir = 'out';
%     ax.YLim = [0,0.6];
%     title('spectrogram over sensor')
% %     export_fig(sprintf('%s_spectrogram.eps',clName))
end

close all
clc
chunks

%% steps  %deal with amplitude later
cn = 'steps';
chunks.(cn)


clName = chunks.(cn).classnames{1};
clN = fastStimMeta.classNum(find(strcmp(cat(1,fastStimMeta.className), clName),1));

for c = 2:length(chunks.(cn).DCoffsets)
      
    
%     %add spectrogram
%     figure
%     spectrogram(stim,10000,500,10000,4e4,'yaxis');    % used to be 20000 - 500 - 20000. Trying to distort less the temporal dimension.
%     colorbar off
%     ax = gca;
%     ax.TickDir = 'out';
%     ax.YLim = [0,0.6];
%     title('spectrogram over command')
%     export_fig(sprintf('%s_%d_spectrogram.eps',clName, c))
end

close all
clc
chunks

%% pips  %deal with amplitude later
cn = 'pips';
chunks.(cn)

edit PipStimulus.m

clName = chunks.(cn).classnames{1};
clN = fastStimMeta.classNum(find(strcmp(cat(1,fastStimMeta.className), clName),1));

for c = 1:length(chunks.(cn).amplitude_3)
    stN = chunks.(cn).amplitude_3(c);
    
    ts = T(clN).ts;
    stim = stimuli(stN).stim.stimulus;
    
    %command prototype
    figure; plot(ts, stim);
    set(gca, 'xlim', [ts(1), ts(end)])
    set(gca, 'ylim', [-1,1])
    set(gca, 'TickDir', 'out')
    export_fig(sprintf('%s_%d_prototype_ampl3.eps',clName, c))
    
    %sensor with amplitude
    figure; hold on; title('sensor')
    DC_offset = mean(T(clN).sensor(stN,1e2:1:.5e4));
    sensor_um = (T(clN).sensor(stN,:) - DC_offset) * 3;
    plot(ts, sensor_um);
    set(gca, 'TickDir', 'out')
    set(gca, 'xlim', [ts(1), ts(end)])
    set(gca, 'ylim', [-2,2])
    ylabel('um')
    export_fig(sprintf('%s_%d_sensor_ampl3.eps',clName, c))
    
    %add spectrogram
    figure
    spectrogram(stim,10000,500,10000,4e4,'yaxis');    % used to be 20000 - 500 - 20000. Trying to distort less the temporal dimension.
    colorbar off
    ax = gca;
    ax.TickDir = 'out';
    ax.YLim = [0,0.6];
    title('spectrogram over command')
    export_fig(sprintf('%s_%d_spectrogram_ampl3.eps',clName, c))
    close all
end



for c = 4
    stN = chunks.(cn).amplitude_2(c);
    
    ts = T(clN).ts;
    stim = stimuli(stN).stim.stimulus;
    
    %command prototype
    figure; plot(ts, stim);
    set(gca, 'xlim', [ts(1), ts(end)])
    set(gca, 'ylim', [-1,1])
    set(gca, 'TickDir', 'out')
    export_fig(sprintf('%s_%d_prototype_ampl2.eps',clName, c))
    
    %sensor with amplitude
    figure; hold on; title('sensor')
    DC_offset = mean(T(clN).sensor(stN,1e2:1:.5e4));
    sensor_um = (T(clN).sensor(stN,:) - DC_offset) * 3;
    plot(ts, sensor_um);
    set(gca, 'TickDir', 'out')
    set(gca, 'xlim', [ts(1), ts(end)])
    set(gca, 'ylim', [-2,2])
    ylabel('um')
    export_fig(sprintf('%s_%d_sensor_ampl2.eps',clName, c))
end
for c = 4
    stN = chunks.(cn).amplitude_1(c);
    
    ts = T(clN).ts;
    stim = stimuli(stN).stim.stimulus;
    
    %command prototype
    figure; plot(ts, stim);
    set(gca, 'xlim', [ts(1), ts(end)])
    set(gca, 'ylim', [-1,1])
    set(gca, 'TickDir', 'out')
    export_fig(sprintf('%s_%d_prototype_ampl1.eps',clName, c))
    
    %sensor with amplitude
    figure; hold on; title('sensor')
    DC_offset = mean(T(clN).sensor(stN,1e2:1:.5e4));
    sensor_um = (T(clN).sensor(stN,:) - DC_offset) * 3;
    plot(ts, sensor_um);
    set(gca, 'TickDir', 'out')
    set(gca, 'xlim', [ts(1), ts(end)])
    set(gca, 'ylim', [-2,2])
    ylabel('um')
    export_fig(sprintf('%s_%d_sensor_ampl1.eps',clName, c))
end

close all
clc
chunks

%% tones  %deal with amplitude later
cn = 'tones';
chunks.(cn)


clName = chunks.(cn).classnames{1};
clN = fastStimMeta.classNum(find(strcmp(cat(1,fastStimMeta.className), clName),1));

for c = 1:length(chunks.(cn).amplitude_2_carrPhase_1)
    stN = chunks.(cn).amplitude_2_carrPhase_1(c);
    
    ts = T(clN).ts;
    stim = stimuli(stN).stim.stimulus;
    
    %command prototype
    figure; plot(ts, stim);
    set(gca, 'xlim', [ts(1), ts(end)])
    set(gca, 'ylim', [-2,2])
    set(gca, 'TickDir', 'out')
    export_fig(sprintf('%s_%d_prototype.eps',clName, c))
    
    %sensor with amplitude
    figure; hold on; title('sensor')
    DC_offset = mean(T(clN).sensor(c,1e2:1:.5e4));
    sensor_um = (T(clN).sensor(1,:) - DC_offset) * 3;
    plot(ts, sensor_um);
    set(gca, 'TickDir', 'out')
    set(gca, 'xlim', [ts(1), ts(end)])
    set(gca, 'ylim', [-5,5])
    ylabel('um')
    export_fig(sprintf('%s_%d_sensor.eps',clName, c))
    
    %add spectrogram
    figure
    spectrogram(stim,10000,500,10000,4e4,'yaxis');    % used to be 20000 - 500 - 20000. Trying to distort less the temporal dimension.
    colorbar off
    ax = gca;
    ax.TickDir = 'out';
    ax.YLim = [0,0.6];
    title('spectrogram over command')
    export_fig(sprintf('%s_%d_spectrogram.eps',clName, c))
end

close all
clc


