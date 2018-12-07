%% GET MUA ACTIVITY FROM EACH SHANK

load CellParams.mat

ShankID = cell2mat({CellParams.ShankID}'); 

%nShanks = length(unique(ShankID)); 

spiketimes = ({CellParams.SpikeTimes}');

MUA = cell(1,max(ShankID));

for i = 1:length(ShankID)
      MUA{ShankID(i)} = [MUA{ShankID(i)}; cell2mat(spiketimes(ShankID(i)))];
end

for j = 1:length(MUA)
    MUA{j} = sort(MUA{j}); 
    MUA{j} = unique(MUA{j});
end

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
%% BIN MUA FOR EACH SHANK

for i = 1:length(MUA)
    [h(:,i),t] = binspikes(MUA{i},1/bin_size,[0 rec_length]);  % requires Chronux toolbox
    %h(:,i) = hist(MUA{i}, nBins); % otherwise use this line 
end

emp = [];
for j = 1:size(h,2)
    if ~any(h(:,j))
        emp = [emp j]; 
    end
end
h(:,emp) = []; % get rid of shanks where there weren't any spikes 

%% define parameters for moving window

win_size = 0.05; % 50 ms min packet size in literature
win_size = win_size/bin_size;

step_size = 0.01; % 10 ms == 1 bin
step_size = step_size/bin_size; 
thresh = 3; % number of simultanoeously active shanks

%% run a moving window 

syn = zeros(1,((floor(size(h,1) - win_size)/step_size +1)));

for nn = 1:(floor(size(h,1) - win_size)/step_size +1)
    ini = (nn-1)*step_size+1;
    fin = ini + win_size-1;
    
    for s = 1:size(h,2)
        ss(s) = sum(h(ini:fin,s));
    end
    
    temp = find(ss ~= 0); 
    
    if length(temp) >= thresh
        a=1;
    else 
        a=0;
    end
    syn(nn) = a;
end

%% detect onsets and offsets 

clear INX alls on off 

f = find((syn) == 0) ;

syn = syn(f(1):f(end)); % has to start from 0 for the detection to work

thr = 0.9;  % Specify threshold 
idxl = syn>=thr;
idxl(1) = 0;
idx = find(idxl);
yest = syn(idx-1)<thr;  % onset
yest2 = syn(idx+1)<thr; % offset 
on=idx(yest)';  % Final output
off=idx(yest2)';
plot(syn); % was reduce_plot(syn) requires tuckermcclure-matlab-plot-big toolbox / use plot instead 
hold on;
alls = sort(vertcat(on, off)); 
alls(:,2)= thr;
scatter(alls(:,1), alls(:,2), 'filled');
ylim([-0.1 2]);
alls(:,2)=[];
INX(:,1) = on; 
INX(:,2) = off; 

INX = INX + f(1); % add the time we cut so indices are aligned properly 

%% 

Fs = 1/(rec_length/nn) 
basepath = pwd;
basename = bz_BasenameFromBasepath(basepath);
filename = [basepath filesep basename '.evt.uds'];

n = size(INX,1); 
UDS=(INX./Fs)'; 
channelID = 1;
 
events.time = UDS(:); 

for i = 1:2:2*n,
	events.description{i,1} = [' start ' int2str(channelID)];
	events.description{i+1,1} = ['stop ' int2str(channelID)];
end

SaveEvents(filename,events);