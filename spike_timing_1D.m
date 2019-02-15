% compute mcc measure from Luczak 2009

load CellParams.mat
load StateIndex.mat
% spiketimesArray = cell2mat({CellParams.SpikeTimes}');
n_neu = size(CellParams,2);

%% only keep packets with a duratiom in specified interval

min_dur = .05; % in sec
max_dur = .1;

for i=1:size(UDS,2)
    dif = UDS(2,i) - UDS(1,i);
    if dif < min_dur | dif > max_dur
        UDS(2,i) = 0;
        UDS(1,i) = 0;
    end
end

UDS(:,UDS(1,:)==0) = [];
n_pac = size(UDS,2);

%% create 2D cell array neurons*packets

cellTot = cell(n_neu,n_pac);
for i=1:n_neu
    spiketimes = cell2mat({CellParams(i).SpikeTimes}');
    for j=1:n_pac
        % take only spikes in the corresponding packet
        temp = spiketimes(spiketimes<UDS(2,j) & spiketimes>UDS(1,j));
        % reference should be beginning of packet
        temp = temp - UDS(1,j);
        cellTot(i,j) = {temp};
        % could also merge in a 1D cell array, but this way it is tidier
    end
end

%% Compute crosscorrelograms and their center of mass

all = 1:n_neu;
com = zeros(1,n_neu);
for i = 1:n_neu
    spikes = [];
    % not very efficient but could not find my way with cell..
    for j = 1:n_pac
        spikes = [spikes cellTot{i,j}'];
    end
    other = setdiff(all,i);
    other_spikes = [];
    for k=1:n_neu-1
        num = other(k);
        for j = 1:n_pac
            other_spikes = [other_spikes cellTot{num,j}'];
        end
    end
    spikes = sort(spikes);
    other_spikes = sort(other_spikes);
    cch = CrossCorr(spikes, other_spikes, 0.002, 100);
    cch = cch./0.002./length(spikes);
    timevec = linspace(-.1 , .1, 101);
    figure
    bar(timevec,smooth(cch))
    com(i) = mean(cch.*timevec')/sum(cch);
end