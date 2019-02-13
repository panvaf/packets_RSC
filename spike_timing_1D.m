% compute ?cc measure from Luczak 2009

load CellParams.mat
load StateIndex.mat
spiketimesArray = cell2mat({CellParams.SpikeTimes}');
n_neu = size(CellParams.SpikeTimes,2);

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
        temp = temp - UDS(2,j);
        cellTot(i,j) = mat2cell(temp,1,length(temp));
        % could also merge in a 1D cell array, but this way it is tidier
    end
end

%% Compute crosscorrelograms and their center of mass

all = 1:n_neu;
com = zeros(1,n_neu);
for i = 1:n_neu
    spikes = cell2mat({cellTot(i,:)}');
    other = setdiff(all,i);
    other_spikes = cell2mat({cellTot(other,:)}');
    cch = CrossCorr(spikes, other_spikes, 0.002, 50);
    cch = cch./0.002./length(spikes);
    com(i) = mean(cch); % need to modify because it is distribution not observations
end