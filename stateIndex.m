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
window = .05; % in sec, size of window to look for percentage of silent bins
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

time = bin_size:bin_size:rec_length;
figure
plot(time,perc_active)
title('State index')
xlabel('Time (s)')
ylabel('Percentage')