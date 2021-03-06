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
window = .02; % in sec, size of window to look for percentage of silent bins
mult = 100; % averaging window how much larger?
win = mult*window; % in sec, for non stationarity of state index
per_win = 5; % window to compute percentage of silent states
nn = length(CellParams);
th = nn/100; % threshold below which it is considered silent, should depend on number of neurons
per_win = floor(per_win/bin_size);
nBins = floor(rec_length/bin_size);
window = floor(window/bin_size);
win = floor(win/bin_size);
kernel = gausswin(window);
th_syn = 33;
th_desyn = 67;  % as percentage on sorted values of state index

%% Compute MUA

spikes_hist = hist(spiketimesArray,nBins);
clear spiketimesArray

MUA = filter(kernel,1,spikes_hist);
perc_active = ones(1,nBins);

% with moving average
for i = ceil(per_win/2):nBins - floor(per_win/2)
    perc_active(i) = sum(MUA((i-floor(per_win/2)+1):(i+floor(per_win/2)))>th)/per_win;
end

% perc_active = lowpass(perc_active,50,Fs); % no signal analysis toolbox?

a = prctile(perc_active,[th_syn th_desyn]);
th_syn = a(1); th_desyn = a(2);

%% define syn and desyn states

syn = zeros(1,nBins);
desyn = zeros(1,nBins);

for i=1:length(perc_active)
    if perc_active(i) > th_desyn
        desyn(i) = 1;
    elseif perc_active(i) < th_syn
        syn(i) = 1;
    end
end

f = find((syn) == 0) ;
syn = syn(f(1):f(end)); % has to start from 0 for the detection to work
transitions = diff(syn);
on_syn = find(transitions==1)';  % Final output
off_syn = find(transitions==-1)';

z = find((desyn) == 0);
desyn = desyn(z(1):z(end)); % has to start from 0 for the detection to work
transitions = diff(desyn);
on_desyn = find(transitions==1)';  % Final output
off_desyn = find(transitions==-1)';

% z(1) = 1; f(1) = 1; % deactivate feature

INX_syn(:,1) = on_syn;
INX_syn(:,2) = off_syn;
INX_desyn(:,1) = on_desyn;
INX_desyn(:,2) = off_desyn;
INX_syn = INX_syn + f(1) - 1 - window/2;
INX_desyn = INX_desyn + z(1) - 1 - window/2;

%% merge successive syn and desyn states to get larger windows

temp_syn = ones(size(INX_syn,1),3);
temp_syn(:,1:2) = INX_syn;
temp_desyn = zeros(size(INX_desyn,1),3);
temp_desyn(:,1:2) = INX_desyn;
full = [temp_syn; temp_desyn];
clear temp_syn temp_desyn

full = sortrows(full);

counter = 1; state = 1;
for i = 1:size(full,1)
    if state==full(i,3)
        full(counter,2) = full(i,2);
        full(i,1) = 0;
    else
        state = full(i,3);
        counter = i;
    end
end

del = full(:,1)==0;
full(del,:) = [];

n_syn = sum(full(:,3));
n_desyn = size(full,1) - n_syn;

full = sortrows(full,[3 1]);
INX_desyn = full(1:n_desyn,1:2);
INX_syn = full((n_desyn+1):end,1:2);
clear full

%% plot

time = bin_size:bin_size:rec_length;
figure
plot(time,MUA/sum(kernel))
hold on
plot(time,perc_active)
plot(time,th_syn*ones(size(time)))
plot(time,th_desyn*ones(size(time)))
title('State index')
xlabel('Time (s)')
ylabel('Percentage')
scatter(INX_syn(:,1)/1000,th_syn*ones(size(INX_syn,1),1))
scatter(INX_syn(:,2)/1000,th_syn*ones(size(INX_syn,1),1))
scatter(INX_desyn(:,1)/1000,th_desyn*ones(size(INX_desyn,1),1))
scatter(INX_desyn(:,2)/1000,th_desyn*ones(size(INX_desyn,1),1))
legend('MUA','State Index','Syn threshold','Desyn threshold','syn-start','syn-end','desyn-start','desyn-end')

%% 

Fs = 1/(rec_length/nBins);
basepath = pwd;
basename = bz_BasenameFromBasepath(basepath);
filename_desyn = [basepath filesep basename '.evt.des'];
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