
%% import indices of ripples, stimulations and packets

basepath = pwd;
basename = bz_BasenameFromBasepath(basepath);

ripFil = [basepath '/' basename '.evt.rip'];
rip_evs = LoadEvents(ripFil);
rip_st = rip_evs.time(cellfun(@any,regexp(rip_evs.description,'start')));

stimFil = [basepath '/' basename '.evt.stm'];
stim_evs = LoadEvents(stimFil);
stim_st = stim_evs.time(cellfun(@any,regexp(stim_evs.description,'start')));

udsFil = [basepath '/' basename '.evt.uds'];
uds_evs = LoadEvents(udsFil);
uds_st = uds_evs.time(cellfun(@any,regexp(uds_evs.description,'start')));
uds_end = uds_evs.time(cellfun(@any,regexp(uds_evs.description,'stop')));

%% 

[rip_status,rip_interval,rip_index] = InIntervals(rip_st, [uds_st uds_end]);

[stim_status,stim_interval,stim_index] = InIntervals(stim_st, [uds_st uds_end]);

rip_on = find(rip_index ~=0); 
rip_off = find(rip_index ==0); 

rip_on_idx = rip_st(rip_on); 
rip_off_idx = rip_st(rip_off); 

stim_on = find(stim_index ~=0); 
stim_off = find(stim_index ==0); 

stim_on_idx = stim_st(stim_on); 
stim_off_idx = stim_st(stim_off); 

stim_idx = {}; 
sinewave_stim_st = 1300; % in seconds, change for every session 
% dsc1904_181015_1 : 1500 -> off responses 
% dsc1904_181015_2 : 900  -> off responses
% dsc1904_181016_1 : 1240 -> no responses
% dsc1904_181016_2 : 1000

% dsc4307_181016_1 : 1300

stim_idx{1} = rip_on_idx; 
stim_idx{2} = rip_off_idx; 
stim_idx{3} = stim_on_idx(stim_on_idx>sinewave_stim_st);
stim_idx{4} = stim_on_idx(stim_on_idx<sinewave_stim_st);
stim_idx{5} = stim_off_idx(stim_off_idx>sinewave_stim_st);
stim_idx{6} = stim_off_idx(stim_off_idx<sinewave_stim_st);
%%
cd ..
chans = 50; % choose a superficial channel, usually this one 
conditions = length(stim_idx); 

basepath = pwd;
basename = bz_BasenameFromBasepath(basepath);

fileinfo = dir([basename '.lfp']);
[xml, ~] = LoadXml(basename); 
Fs = xml.lfpSampleRate;
num_channels = xml.nChannels;
num_samples = fileinfo.bytes/(num_channels * 2); % int16 = 2 bytes
%data = zeros(length(chans), 375, size(stim_idx,1)); % 375
data = cell(1,conditions); 

for l = 1:conditions
    data{l} = zeros(375, size(stim_idx{l},1)); % 375
end

warning off

for c = 1:conditions
    for stimi = 1:size(stim_idx{c},1)
        data{c}(:,stimi) = bz_LoadBinary(fileinfo.name,'frequency',Fs,'nChannels',num_channels, 'channels', chans+1, 'start', stim_idx{c}(stimi) - 0.1,'duration', 0.3).*0.195; % -0.1 -> 0.3 
    end
end


figure(55)
plot(mean(data{1},2)); 
hold on
plot(mean(data{2},2)); 
plot(mean(data{3},2)); 
plot(mean(data{4},2)); 
plot(mean(data{5},2)); 
plot(mean(data{6},2)); 

legend('rip on', 'rip off', 'stim on sine', 'stim on square', 'stim off sine', 'stim off square'); 

%% 

[~, mean_wave_on, period, ~, ~]= morlett_snippet(data{3}');

[~, mean_wave_off, ~, ~, ~]= morlett_snippet(data{5}');

close all

plot_morlett(mean_wave_on,period,1250, 1);

plot_morlett(mean_wave_off,period,1250, 1);
