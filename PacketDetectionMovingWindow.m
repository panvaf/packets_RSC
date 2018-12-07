%% GET MATRIX OF SPIKE TIMES

load CellParams.mat
spiketimesArray = cell2mat({CellParams.SpikeTimes}');

%% DEFINE PARAMETERS FOR BINNIG

basepath = pwd;   
basename = bz_BasenameFromBasepath(basepath); 
fileinfo = dir([basename '.lfp']);
[xml, ~] = LoadXml(basename); 
Fs = xml.lfpSampleRate;
num_channels = xml.nChannels;
num_samples = fileinfo.bytes/(num_channels * 2);
rec_length = num_samples/Fs;

bin_size = 0.01; % in seconds 
nBins = floor(rec_length/bin_size); % 

%% define parameters for moving window

win_size = 0.1; % 50 ms min packet size in literature, integrate over longer time to reduce the effect of coincidences
win_size = win_size/bin_size;

step_size = 0.01; % 10 ms == 1 bin
step_size = step_size/bin_size; 

%% run a moving window 

syn = zeros(1,(((nBins + 1 - win_size)/step_size +1)));
n_spikes = zeros(1,(((nBins + 1 - win_size)/step_size +1)));
spikes_hist = hist(spiketimesArray,nBins);
thresh = 0; % initialize threshold, set to percentage of maximum number of spikes in a window
perc = 0.3; % percentage of maximum number of spikes in a window that corresponds to a packet

for nn = 1:((nBins + 1 - win_size)/step_size-1)
    n_spikes(nn) = sum(spikes_hist(nn:nn+win_size));
    if n_spikes(nn)>thresh
        thresh = n_spikes(nn);
    end
end

thresh = perc*thresh;

for nn = 1:((nBins + 1 - win_size)/step_size +1)
    if n_spikes(nn) >= thresh
        syn(nn)=1;
    end
end

%% detect onsets and offsets 

f = find((syn) == 0) ;
syn = syn(f(1):f(end)); % has to start from 0 for the detection to work
transitions = diff(syn);
on = find(transitions==1)';  % Final output
off = find(transitions==-1)';

plot(syn); % was reduce_plot(syn) requires tuckermcclure-matlab-plot-big toolbox / use plot instead 
hold on;
alls = sort(vertcat(on, off)); 
alls(:,2)= 1;
scatter(alls(:,1), alls(:,2), 'filled');
ylim([-0.1 2]);
alls(:,2)=[];
INX(:,1) = on; 
INX(:,2) = off; 

INX = INX + f(1) - 1; % add the time we cut so indices are aligned properly 

%% 

Fs = 1/(rec_length/nn) 
basepath = pwd;
basename = bz_BasenameFromBasepath(basepath);
filename = [basepath filesep basename '.evt.uds'];

n = size(INX,1); 
UDS=(INX./Fs)'; 
channelID = 1;

events.time = UDS(:); 

for i = 1:2:2*n
	events.description{i,1} = [' start ' int2str(channelID)];
	events.description{i+1,1} = ['stop ' int2str(channelID)];
end

SaveEvents(filename,events);