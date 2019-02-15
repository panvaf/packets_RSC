% nwaydecomp for the characterization of packets using split procedure
% code by original git repository, modifications by panvaf

%% load only essential stuff from fieldtrip

addpath /home/panteleimon/Documents/fieldtrip;
ft_defaults;

%% load data and bring them in appropriate format

dat = load('GG.069_Behavior.mat');
lindist = dat.whlrld(:,7);
% linearized rat distance in maze
fsample = 20000;
% sampling rate in Hz of recording, see documentation
begtrial = round(dat.SessionNP(:,2) * (fsample./512)); % convert to lindist sampling rate
endtrial = round(dat.SessionNP(:,3) * (fsample./512)); % convert to lindist sampling rate
ntrial = numel(begtrial);
data = cell(1,ntrial);
% correct end of trial to last known position
for itrial = 1:ntrial
    currldist = lindist(begtrial(itrial):endtrial(itrial));
    nzind = find(currldist,1,'last');
    endtrial(itrial) = begtrial(itrial) -1 + nzind;
end
% ensure start/end of trial are closest to beg/end of maze
for itrial = 1:ntrial
    currldist = lindist(begtrial(itrial):endtrial(itrial));
    endtrial(itrial) = begtrial(itrial) -1 + find(currldist>=currldist(end),1);
    begtrial(itrial) = begtrial(itrial) -1 + find(currldist>currldist(1),1);
end
begtrial = begtrial .* 512; % convert back to timestamp sampling rate (20kHz)
endtrial = endtrial .* 512; % convert back to timestamp sampling rate (20kHz)
% create spike trains
for itrial = 1:ntrial
    % fetch timestamps and neuren IDs for current trial
    ind = dat.spiket >= begtrial(itrial) & dat.spiket <= endtrial(itrial);
    currts = dat.spiket(ind);
    currts = currts - begtrial(itrial) + 1; % convert to trial specific samples
    neurid = dat.spikeind(ind);
    % create sparse spike train matrix
    nsample = (endtrial(itrial)-begtrial(itrial))+1;
    data{itrial} = sparse(neurid,currts,ones(size(currts)),dat.numclus,nsample);
end
% bookkeeping variables
trialtype = dat.SessionNP(:,4); % 1 (right) or 2 (left) lap
label = []; % neuronIDs
for inneuron = 1:dat.numclus
    label{inneuron} = ['neuron' num2str(inneuron)];
end
% select neurons with mean spike rate above 1Hz
trialtime = cellfun(@size,data,repmat({2},[1 ntrial])) ./ fsample;
nspikes = cellfun(@sum,data,repmat({2},[1 ntrial]),'uniformoutput',0);
nspikes = cat(2,nspikes{:});
spikerate = mean(bsxfun(@rdivide,nspikes,trialtime),2);
selind = spikerate>1;
for itrial = 1:ntrial
    data{itrial} = data{itrial}(selind,:);
end
label = label(selind);

clearvars -except data label
% from now on everything should be the same. I only need to bring our data
% to the format of "data"

%% split data

% obtain splits
nneuron = size(data{1},1); % data is the cell-array of spike trains we created above
ntrial = numel(data);
dataodd = cell(1,ntrial);
dataeven = cell(1,ntrial);
for itrial = 1:ntrial
% get spike time stamps/indices
    nsample = size(data{itrial},2);
    [neurid, ts] = ind2sub(size(data{itrial}),find(data{itrial}));
    oddind = 1:2:numel(ts);
    evenind = 2:2:numel(ts);
    % create new sparse spike train matrices
    dataodd{itrial} = sparse(neurid(oddind), ts(oddind), ones(size(oddind)),nneuron,nsample);
    dataeven{itrial} = sparse(neurid(evenind),ts(evenind),ones(size(evenind)),nneuron,nsample);
end

%% compute cross spectra

% set cross spectra choices
timwin = 0.020; % complex exponential length in seconds
freq = 50:50:1000; % frequency in Hz
fsample = 20000; % sampling rate of the spike trains in Hz
% set n's
ntrial = numel(dataodd);
nneuron = size(dataodd{1},1);
nfreq = numel(freq);
% pre-allocate
fourierodd = NaN(nneuron,nfreq,ntrial,nneuron,'single');
fourierodd = complex(fourierodd,fourierodd);
% convolve with complex exponential, obtain cross spectrum, obtain square root
for itrial = 1:ntrial
    for ifreq = 1:nfreq
        disp(['working on trial #' num2str(itrial) ' @' num2str(freq(ifreq)) 'Hz'])
        % construct complex exponential
        timeind = (-(round(timwin .* fsample)-1)/2 : (round(timwin .* fsample)-1)/2) ./fsample;
        wlt = exp(1i*timeind*2*pi*freq(ifreq));
        % get spectrum via sparse convolution
        spectrum = sconv2(dataodd{itrial},wlt,'same');
        % obtain cross spectrum, csd
        csd = full(spectrum*spectrum');
        % normalize csd by number of time-points (for variable length epochs)
        csd = csd ./ (size(dataodd{itrial},2) ./ fsample);
        % obtain square root of cross spectrum
        [V L] = eig(csd);
        L = diag(L);
        tol = max(size(csd))*eps(max(L));
        zeroL = L<tol;
        eigweigth = V(:,~zeroL)*diag(sqrt(L(~zeroL)));
        % save in fourier
        currm = size(eigweigth,2);
        fourierodd(:,ifreq,itrial,1:currm) = eigweigth;
    end
end

% do the same for even spikes

% set n's
ntrial = numel(dataeven);
nneuron = size(dataeven{1},1);
nfreq = numel(freq);
% pre-allocate
fouriereven = NaN(nneuron,nfreq,ntrial,nneuron,'single');
fouriereven = complex(fouriereven,fouriereven);
% convolve with complex exponential, obtain cross spectrum, obtain square root
for itrial = 1:ntrial
    for ifreq = 1:nfreq
        disp(['working on trial #' num2str(itrial) ' @' num2str(freq(ifreq)) 'Hz'])
        % construct complex exponential
        timeind = (-(round(timwin .* fsample)-1)/2 : (round(timwin .* fsample)-1)/2) ./fsample;
        wlt = exp(1i*timeind*2*pi*freq(ifreq));
        % get spectrum via sparse convolution
        spectrum = sconv2(dataeven{itrial},wlt,'same');
        % obtain cross spectrum, csd
        csd = full(spectrum*spectrum');
        % normalize csd by number of time-points (for variable length epochs)
        csd = csd ./ (size(dataeven{itrial},2) ./ fsample);
        % obtain square root of cross spectrum
        [V L] = eig(csd);
        L = diag(L);
        tol = max(size(csd))*eps(max(L));
        zeroL = L<tol;
        eigweigth = V(:,~zeroL)*diag(sqrt(L(~zeroL)));
        % save in fourier
        currm = size(eigweigth,2);
        fouriereven(:,ifreq,itrial,1:currm) = eigweigth;
    end
end

%% Normalize cross-spectra
N = 32; % an example N
[nneuron,nfreq,ntrial,ntime] = size(fourierodd);
% compute sum of power
power = zeros(nneuron,1);
for itrial = 1:ntrial
    currfour = double(squeeze(fourierodd(:,:,itrial,:)));
    power = power + nansum(nansum(abs(currfour).^2,3),2);
end
% compute scaling such that power is its kth root
scaling = (power .^ (1/N)) ./ power;
for itrial = 1:ntrial
    for ifreq = 1:nfreq
        currfour = double(squeeze(fourierodd(:,ifreq,itrial,:)));
        nonnanind = ~isnan(currfour(1,:));
        currfour = currfour(:,nonnanind);
        currfour = bsxfun(@times,currfour,sqrt(scaling));
        % save in fourier
        fourierodd(:,ifreq,itrial,nonnanind) = single(currfour);
    end
end

% do the same for even spikes

[nneuron,nfreq,ntrial,ntime] = size(fouriereven);
% compute sum of power
power = zeros(nneuron,1);
for itrial = 1:ntrial
    currfour = double(squeeze(fouriereven(:,:,itrial,:)));
    power = power + nansum(nansum(abs(currfour).^2,3),2);
end
% compute scaling such that power is its kth root
scaling = (power .^ (1/N)) ./ power;
for itrial = 1:ntrial
    for ifreq = 1:nfreq
        currfour = double(squeeze(fouriereven(:,ifreq,itrial,:)));
        nonnanind = ~isnan(currfour(1,:));
        currfour = currfour(:,nonnanind);
        currfour = bsxfun(@times,currfour,sqrt(scaling));
        % save in fourier
        fouriereven(:,ifreq,itrial,nonnanind) = single(currfour);
    end
end

%% run SPACE

% create input structure for nd_nwaydecomposition
fourierdata = [];
fourierdata.fourier = fourierodd;
% square root of cross spectra of all spikes
fourierdata.fouriersplit{1} = fourierodd;
% square root of cross spectra of odd spikes
fourierdata.fouriersplit{2} = fouriereven;
% square root of cross spectra of even spikes
fourierdata.freq = freq;
% frequencies in Hz
fourierdata.label = label;
% neuron IDs
fourierdata.dimord = 'chan_freq_epoch_tap'; % req. for SPACE due to other uses
% fourierdata.trialinfo = trialtype;
% left/right trial codes
fourierdata.cfg = [];
% req. for SPACE due to other uses
% extract spike timing networks
cfg = [];
cfg.model = 'spacetime';
cfg.datparam = 'fourier';
cfg.ncompestsrdatparam = 'fouriersplit';
cfg.Dmode = 'identity';
cfg.ncompest = 'splitrel';
cfg.ncompestsrcritjudge = 'meanoversplitscons';
cfg.ncompeststart = 10;
cfg.ncompeststep = 5;
cfg.ncompestend = 50;
cfg.ncompestsrcritval = [.7 0 0 .7 0];
cfg.numiter = 1000;
cfg.convcrit = 1e-6;

% max number of iteration
% stop criterion of algorithm: minimum
%
cfg.randstart = 50;
% number of random initializations
cfg.ncompestrandstart = 50;
cfg.distcomp.system = 'torque';
% number of random init. for split-rel. proc.
nwaycomp = nd_nwaydecomposition(cfg,fourierdata);