%% GET MATRIX OF SPIKE TIMES
include_MUA = false;
add_random_jitter = false; % better to have MUA on false when this is true
jitter_range = 1; % in s

load CellParams.mat

if add_random_jitter
    spiketimesArray = [];
    for i=1:size(CellParams,2)
        spiketimes = CellParams(i).SpikeTimes + jitter_range*rand();
        spiketimesArray = [spiketimesArray; spiketimes];
    end
else
    spiketimesArray = cell2mat({CellParams.SpikeTimes}');
end

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
rec_length = num_samples/Fs;
spiketimesArray(spiketimesArray>rec_length) = []; % needed because of added jitter

% in seconds 
bin_size = 0.001;
nBins = floor(rec_length/bin_size);

%% Parameters for detection, need to be tuned separately for each dataset
kernel_size = 0.04; % in sec, gaussian kernel
<<<<<<< HEAD
per_act = .8; % percentage of average local population activity that corresponds to a packet
win = 10*60; % size of local window, in seconds
isclose = 0.035; % in sec, merge events that are too close
min_duration = 0.03; % in sec, look for silent periods before and after packet
per_sil = .3; % percentage of maximum of packet that must be reached by a silent in-between state
min_dur = 0.04; % in sec, delete packets that are smaller than this
max_dur = 0.7; % in sec, delete packets that are bigger than this
=======
per_act = .6; % percentage of average local population activity that corresponds to a packet
win = 5*60; % size of local window, in seconds
isclose = 0.035; % in sec, merge events that are too close
min_duration = 0.05; % in sec, look for silent periods before and after packet
per_sil = .2; % percentage of maximum of packet that must be reached by a silent in-between state
min_dur = 0.04; % in sec, delete packets that are smaller than this
max_dur = 0.5; % in sec, delete packets that are bigger than this
>>>>>>> upstream/master

% initial params (Tuesday): .1,1.5,.025,.05,.4,.02
% got quite good diagrams (Wednesday) with 0.07, 1.3, 0.025, 0.05, 0.3, 0.02
% better detection with 0.07, 1.2, 0.025, 0.05, 0.2, 0.02, still some
% packets missing

% Noam's parameters (Thursday): 0.02, 1, .025, .05, .05, .02, .5 a lot of false
% positives

% 0.02, 1.2, .025, .05, .05, .02, .5 is basically a coincidence detector

% 0.04, 1, 0.035, .05, .08, .04, .5 is a good compromise between picking
% up small packets and avoiding coincidence detection. could lower th a
% little more if we want to be more robust at the cost of losing packets

% LR1: more strict threshold on silence because a lot of asynchronous
% activity? chose to increase the threshold for detection:
% 0.04, 1.2, 0.035, .05, .08, .04, .5

%% define kernel

% 50 ms min packet size in literature, integrate over longer time to reduce the effect of coincidences
kernel_size = floor(kernel_size/bin_size);
kernel = gausswin(kernel_size);

%% convolve with kernel 

win = floor(win/bin_size);
syn = zeros(1,nBins);
spikes_hist = hist(spiketimesArray,nBins);
clear spiketimesArray

pop = filter(kernel,1,spikes_hist);
% pop = zscore(pop); % does not increase discrimination in a meaningfull way

mean_win = movmean(pop,win);
thresh = per_act*mean_win; % set threshold to percentage local average

for i=1:size(pop,2)
    if pop(i) > thresh(i)
        syn(i) = 1;
    end
end

%% detect onsets and offsets 

f = find((syn) == 0) ;
syn = syn(f(1):f(end)); % has to start from 0 for the detection to work
transitions = diff(syn);
on = find(transitions==1)';  % Final output
off = find(transitions==-1)';

%% if two packets are too close, collate them

isclose = floor(isclose/bin_size);

for i=2:size(off)
    if on(i)-off(i-1)<isclose
        on(i) = 0;
        off(i-1) = 0;
    end
end

on(on==0) = [];
off(off==0) = [];

%% if a packet is not happening between periods of inactivity, delete it
% must always be done after packet collation! (because otherwise we might lose part of a packet)

min_duration = floor(min_duration/bin_size);
th = per_sil*mean_win;

for i=2:size(off)-1
    min_bef = min(pop((on(i)-min_duration):on(i)));
    min_aft = min(pop(off(i):(off(i)+min_duration)));
    if min_bef > th(on(i)) | min_aft > th(off(i))
        on(i) = 0;
        off(i) = 0;
    end
end

on(on==0) = [];
off(off==0) = [];

%% delete packets that are too small or too large

min_dur = floor(min_dur/bin_size);
max_dur = floor(max_dur/bin_size);

for i=1:size(off)
    dif = off(i) - on(i);
    if dif < min_dur | dif > max_dur
        on(i) = 0;
        off(i) = 0;
    end
end

on(on==0) = [];
off(off==0) = [];

%%
plot(syn); % was reduce_plot(syn) requires tuckermcclure-matlab-plot-big toolbox / use plot instead 
hold on;
alls = sort(vertcat(on, off)); 
alls(:,2)= 1;
scatter(alls(:,1), alls(:,2), 'filled');
ylim([-0.1 2]);
alls(:,2)=[];
INX(:,1) = on; 
INX(:,2) = off; 

INX = INX + f(1) - 1 - kernel_size/2; % add the time we cut so indices are aligned properly ,

% add offset of ~50ms in DSC4307_181016_1_RSC file

%% 
Fs = 1/(rec_length/nBins)
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