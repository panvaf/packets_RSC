% compute mcc measure from Luczak 2009

load CellParams.mat
load StateIndex.mat

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

%% find total numbers of spikes for each neuron to preallocate memory
spike_n = zeros(n_neu,1);

parfor i=1:n_neu
    spiketimes = cell2mat({CellParams(i).SpikeTimes}');
    for j=1:n_pac
        % take only spikes in the corresponding packet
        temp = spiketimes(spiketimes<UDS(2,j) & spiketimes>UDS(1,j));
        % reference should be beginning of packet
        spike_n(i) = spike_n(i) + length(temp);
    end
end


%% create 1D cell array with spikes from every neuron

cellTot = cell(n_neu);
parfor i=1:n_neu
    counter = 1;
    cur_spikes = zeros(spike_n(i),1);
    spiketimes = cell2mat({CellParams(i).SpikeTimes}');
    for j=1:n_pac
        % take only spikes in the corresponding packet
        temp = spiketimes(spiketimes<UDS(2,j) & spiketimes>UDS(1,j));
        % reference should be beginning of packet
        temp = temp - UDS(1,j);
        num = length(temp);
        if num
            cur_spikes(counter:(counter+num-1)) = temp';
        end
        counter = counter + num;
    end
    cellTot{i} = cur_spikes;
end

%% Compute crosscorrelograms and their center of mass

all = 1:n_neu;
com = zeros(1,n_neu);
crosscors = zeros(101,n_neu);
parfor i = 1:n_neu
    spikes = cellTot{i}';
    other = setdiff(all,i);
    other_spikes = zeros(sum(spike_n)-spike_n(i),1);
    counter = 1;
    for j=1:n_neu-1
        pointer = other(j);
        num = spike_n(pointer);
        other_spikes(counter:(counter+num-1)) = cellTot{pointer}';
        counter = counter + num;
    end
    spikes = sort(spikes);
    other_spikes = sort(other_spikes);
    cch = CrossCorr(spikes, other_spikes, 0.002, 100);
    cch = cch./0.002./length(spikes);
    timevec = linspace(-.1 , .1, 101);
    cor = timevec(end)./(timevec(end)*1.02 - abs(timevec));
    cch_cor = cch.*cor';
    crosscors(:,i) = cch_cor;
    % figure; bar(timevec,smooth(cch_cor))
    com(i) = mean(cch_cor.*timevec')/sum(cch_cor);
    % still would have to correct for varying packet length (big lags not favoured)
end

figure
plot(com*1000)
ylabel('CrossCor center of mass (ms)')
xlabel('Neuron #')