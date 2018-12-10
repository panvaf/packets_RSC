
%% GET MATRIX OF SPIKE TIMES
include_MUA = false;
add_random_jitter = false; % better to have MUA on false when this is true
jitter_range = 1; % in s

load CellParams.mat
CellParams([85, 91, 106]) = []; % for LR1 noisy cells (Int?)

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
kernel_size = 0.01; % in sec, gaussian kernel
isclose = 0.03; % in sec, merge events that are too close
min_duration = 0.02; % minumun duration of silence before and after packets
th = .5; % max allowed percentage of activity during off periods
min_dur = 0.04; % min duration of packets allowed
max_dur = 0.7; % max duration of packets allowed

% define kernel

kernel_size = kernel_size/bin_size;
kernel = gausswin(kernel_size);
%%

spikes_hist = hist(spiketimesArray,nBins);
clear spiketimesArray
pval = 0.1;
zval = abs(norminv(pval));

pop = filter(kernel,1,spikes_hist);
pop=smooth(pop, 300, 'sgolay')'; 

zmap = zscore(pop);
zmap(abs(zmap)<zval) = 0;
zmap(zmap<0) = 0; 

islands = bwconncomp(zmap);
a = {islands.PixelIdxList{1:end}};
b = []; 
for i = 1:length(a)
    temp = length(a{i}); 
    if temp >= min_dur/bin_size && temp <= max_dur/bin_size % ms of silent period
        b = [b i]; 
    end
end

new_islands = a(b);

m = [];
for j = 1:length(b)-1
    m(j,1) = new_islands{j}(1); 
    m(j,2) = new_islands{j}(end); 
end

D = abs(m(:,2) - m(:,1)); 
    
dur_ok = find(D>min_dur/bin_size & D<max_dur/bin_size);  

on = m(dur_ok,1); 
off = m(dur_ok,2); 

%% if two packets are too close, collate them


for i=2:size(off)
    if on(i)-off(i-1)<isclose/bin_size;
        on(i) = 0;
        off(i-1) = 0;
    end
end

on(on==0) = [];
off(off==0) = [];
%% delete packets that are not prominent 


for i=2:size(off)-1
    max_in = max(pop(on(i):off(i)));
    min_bef = mean(pop((on(i)-min_duration/bin_size):on(i)));
    min_aft = mean(pop(off(i):(off(i)+min_duration/bin_size)));
    if min_bef > th*max_in | min_aft > th*max_in
        on(i) = 0;
        off(i) = 0;
    end
end

on(on==0) = [];
off(off==0) = [];

%% 

INX(:,1) = on; 
INX(:,2) = off; 

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





