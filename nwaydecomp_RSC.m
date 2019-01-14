% nwaydecomp for the characterization of packets
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

%% compute cross spectra

% one variable taken from above
data; % cell-array of spike train matrix per trial
% set cross spectra choices
timwin = 0.020; % complex exponential length in seconds
freq = 50:50:1000; % frequency in Hz
fsample = 20000; % sampling rate of the spike trains in Hz
% set n's
ntrial = numel(data);
nneuron = size(data{1},1);
nfreq = numel(freq);
% pre-allocate
fourier = NaN(nneuron,nfreq,ntrial,nneuron,'single');
fourier = complex(fourier,fourier);
% convolve with complex exponential, obtain cross spectrum, obtain square root
for itrial = 1:ntrial
    for ifreq = 1:nfreq
        disp(['working on trial #' num2str(itrial) ' @' num2str(freq(ifreq)) 'Hz'])
        % construct complex exponential
        timeind = (-(round(timwin .* fsample)-1)/2 : (round(timwin .* fsample)-1)/2) ./fsample;
        wlt = exp(1i*timeind*2*pi*freq(ifreq));
        % get spectrum via sparse convolution
        spectrum = sconv2(data{itrial},wlt,'same');
        % obtain cross spectrum, csd
        csd = full(spectrum*spectrum');
        % normalize csd by number of time-points (for variable length epochs)
        csd = csd ./ (size(data{itrial},2) ./ fsample);
        % obtain square root of cross spectrum
        [V L] = eig(csd);
        L = diag(L);
        tol = max(size(csd))*eps(max(L));
        zeroL = L<tol;
        eigweigth = V(:,~zeroL)*diag(sqrt(L(~zeroL)));
        % save in fourier
        currm = size(eigweigth,2);
        fourier(:,ifreq,itrial,1:currm) = eigweigth;
    end
end
nwaycompallN = cell(5,1);

%% Normalize cross-spectra
Ns = [1 8 16 32 64];

N = 64; % an example N
[nneuron,nfreq,ntrial,ntime] = size(fourier);
% compute sum of power
power = zeros(nneuron,1);
for itrial = 1:ntrial
    currfour = double(squeeze(fourier(:,:,itrial,:)));
    power = power + nansum(nansum(abs(currfour).^2,3),2);
end
% compute scaling such that power is its kth root
scaling = (power .^ (1/N)) ./ power;
for itrial = 1:ntrial
    for ifreq = 1:nfreq
        currfour = double(squeeze(fourier(:,ifreq,itrial,:)));
        nonnanind = ~isnan(currfour(1,:));
        currfour = currfour(:,nonnanind);
        currfour = bsxfun(@times,currfour,sqrt(scaling));
        % save in fourier
        fourier(:,ifreq,itrial,nonnanind) = single(currfour);
    end
end

%% run SPACE

% create input structure for nd_nwaydecomposition
fourierdata = [];
fourierdata.fourier = fourier;
% square root of cross spectra
fourierdata.freq = freq;
% frequencies in Hz
fourierdata.label = label;
% neuron IDs
fourierdata.dimord = 'chan_freq_epoch_tap'; % req. for SPACE due to other uses
fourierdata.trialinfo = trialtype;
% left/right trial codes
fourierdata.cfg = [];
% extract spike timing networks
cfg = [];
cfg.model = 'spacetime';
cfg.datparam = 'fourier';
cfg.Dmode = 'identity';
cfg.ncomp = 5;
cfg.numiter = 1000;
cfg.convcrit = 1e-6;
% the model for extracting spike timing networks
% field containing square root of the cross spectra
% necessary, see background/SPACE ref papers
% number of networks to extract
% max number of iterations
% stop criterion of algorithm: minimum relative
% difference in fit between iterations
cfg.randstart = 25;
% number of random initializations
nwaycomp = nd_nwaydecomposition(cfg,fourierdata);

%% visualize neuron profiles

% plot extracted neuron profiles for all Ns
% first, we assume the Fourier arrays for all Ns are gathered in a cell-array
nwaycompallN{1} % N = 1
nwaycompallN{2} % N = 8
nwaycompallN{3} % N = 16
nwaycompallN{4} % N = 32
nwaycompallN{5} % N = 64
figure
N = [1 8 16 32 64];
nneuron = numel(nwaycompallN{1}.label);
for iN = 1:5
    % plot neuron profiles (contained in the 1st cell)
    subplot(1,5,iN)
    hold on
    for inetw = 1:5
        plot(nwaycompallN{iN}.comp{inetw}{1},'marker','.');
    end
    ylim = get(gca,'ylim');
    set(gca,'ylim',[0 ylim(2)*1.05],'xlim',[1 nneuron])
    xlabel('neurons')
    ylabel('loading')
    title(['neuron profiles N = ' num2str(N(iN))])
    axis square
end