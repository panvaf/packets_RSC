%% GET MATRIX OF SPIKE TIMES
include_MUA = false;

load CellParams.mat

spiketimesArray = cell2mat({CellParams.SpikeTimes}');

if include_MUA
    load MUA.cellinfo.mat
    spiketimesArray = [spiketimesArray; spiketimes];
end

%% DEFINE PARAMETERS FOR BINNING

basepath = pwd;   
basename = bz_BasenameFromBasepath(basepath); 
fileinfo = dir([basename '.dat']);
[xml, ~] = LoadXml(basename); 
Fs = xml.SampleRate;
num_channels = xml.nChannels;
num_samples = fileinfo.bytes/(num_channels * 2);
rec_length = num_samples/Fs;  % in seconds 

%% parameters
bin_size = 0.001;
window = .5; % in sec, size of window to look for percentage of silent bins
th_up = .7; % defines desyn state
th_low = .5; % defines syn state
nBins = floor(rec_length/bin_size);

%% Find percentage of bins that are silent

spiketimes_unique = unique(spiketimesArray);
spike_events = false([nBins 1]);

for i=1:length(spiketimes_unique)
    spike_events(ceil(spiketimes_unique(i)/bin_size)) = true;
end

perc_active = ones(nBins,1);
win_size = floor(window/bin_size);

for i = ceil(win_size/2):nBins - floor(win_size/2)
    perc_active(i) = sum(spike_events((i-floor(win_size/2)+1):(i+floor(win_size/2))))/win_size;
end

peak = max(perc_active(win_size:end-win_size));
th_low = th_low*peak;
th_up = th_up*peak;

time = bin_size:bin_size:rec_length;
figure
plot(time,perc_active)
title('State index')
xlabel('Time (s)')
ylabel('Percentage')
hold on
plot(time,th_up*ones(nBins,1))
plot(time,th_low*ones(nBins,1))


%% define syn and desyn states

syn = zeros(1,nBins);
desyn = zeros(1,nBins);

for i=1:size(perc_active,1)
    if perc_active(i) > th_up
        desyn(i) = 1;
    elseif perc_active(i) < th_low
        syn(i) = 1;
    end
end

f = find((syn) == 0) ;
syn = syn(f(1):f(end)); % has to start from 0 for the detection to work
transitions = diff(syn);
on_syn = find(transitions==1)';  % Final output
off_syn = find(transitions==-1)';

f = find((desyn) == 0) ;
desyn = desyn(f(1):f(end)); % has to start from 0 for the detection to work
transitions = diff(desyn);
on_desyn = find(transitions==1)';  % Final output
off_desyn = find(transitions==-1)';


INX_syn(:,1) = on_syn;
INX_syn(:,2) = off_syn;
INX_desyn(:,1) = on_desyn;
INX_desyn(:,2) = off_desyn;

%% 

Fs = 1/(rec_length/nBins);
basepath = pwd;
basename = bz_BasenameFromBasepath(basepath);
filename_desyn = [basepath filesep basename '.evt.desyn'];
filename_syn = [basepath filesep basename '.evt.syn'];

n = size(INX_syn,1); 
temp=(INX_syn./Fs)'; 
channelID = 1;

events.time = temp(:);  

for i = 1:2:2*n
	events.description{i,1} = [' start ' int2str(channelID)];
	events.description{i+1,1} = ['stop ' int2str(channelID)];
end

SaveEvents(filename_syn,events);

n = size(INX_desyn,1); 
temp=(INX_desyn./Fs)'; 
channelID = 1;

events.time = temp(:);  

for i = 1:2:2*n
	events.description{i,1} = [' start ' int2str(channelID)];
	events.description{i+1,1} = ['stop ' int2str(channelID)];
end

SaveEvents(filename_desyn,events);