%% GET MATRIX OF SPIKE TIMES

load CellParams.mat

%% DEFINE PARAMETERS FOR BINNING

basepath = pwd;   
basename = bz_BasenameFromBasepath(basepath); 
fileinfo = dir([basename '.dat']);
[xml, ~] = LoadXml(basename); 
Fs = xml.SampleRate;
num_channels = xml.nChannels;
num_samples = fileinfo.bytes/(num_channels * 2);
rec_length = num_samples/Fs;


% in seconds 
bin_size = 0.001;
nBins = floor(rec_length/bin_size);

%% get LFP time series for spiking activity frequency window
fr_samp = 20000;
data = bz_LoadBinary(fileinfo.name,'frequency',fr_samp,'nChannels',32,'channels',8);
data = double(data');
[S, F, T, P] = spectrogram(data, 20, [], [], fr_samp);
% imagesc(F,T,log(1+P)) % plot spectogram
fr_low = 300;
fr_high = 700;
len = size(P,1);
LFP = sum(P(floor(len/fr_samp*fr_low):ceil(len/fr_samp*fr_high),:),1);

%% Parameters for detection, need to be tuned separately for each dataset
kernel_size = 0.06; % in sec, gaussian kernel
perc = 1.0; % percentage of average population activity that corresponds to a packet
isclose = 0.025; % in sec, merge events that are too close
min_duration = 0.05; % in sec, look for silent periods before and after packet
th = .15; % percentage of maximum of packet that must be reached by a quiet in-between state
min_dur = 0.02; % in sec, delete packets that are smaller than this

% initial params (Tuesday): .1,1.5,.025,.05,.4,.02
% got quite good diagrams (Wednesday) with 0.07, 1.3, 0.025, 0.05, 0.3, 0.02
% better detection with 0.07, 1.2, 0.025, 0.05, 0.2, 0.02, still some
% packets missing


%% define kernel

% 50 ms min packet size in literature, integrate over longer time to reduce the effect of coincidences
kernel_size = kernel_size/bin_size;
kernel = gausswin(kernel_size);

%% convolve with kernel 

pop = filter(kernel,1,LFP);
thresh = perc*mean(pop); % set threshold to percentage of maximum

syn = zeros(1,nBins);

for i=1:size(pop,2)
    if pop(i) > thresh
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

isclose = isclose/bin_size;

for i=2:size(off)
    if on(i)-off(i-1)<isclose
        on(i) = 0;
        off(i-1) = 0;
    end
end

on(on==0) = [];
off(off==0) = [];

%% if a packet is not prominent enough, delete it
% must always be done after packet collation! (because otherwise we might lose part of a packet)

min_duration = min_duration/bin_size;

for i=2:size(off)-1
    max_in = max(pop(on(i):off(i)));
    min_bef = min(pop((on(i)-min_duration):on(i)));
    min_aft = min(pop(off(i):(off(i)+min_duration)));
    if min_bef > th*max_in | min_aft > th*max_in
        on(i) = 0;
        off(i) = 0;
    end
end

on(on==0) = [];
off(off==0) = [];

%% delete packets that are too small

min_dur = min_dur/bin_size;

for i=1:size(off)
    if off(i) - on(i) < min_dur
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

INX = INX + f(1) - 1 - kernel_size/2 + 50; % add the time we cut so indices are aligned properly 

% systematic offset of ~50ms in DSC4307_181016_1_RSC file, why?

%% 

Fs = 1/(rec_length/nBins);
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