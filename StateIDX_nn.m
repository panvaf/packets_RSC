
%% DEFINE PARAMETERS FOR BINNING
load CellParams.mat

include_MUA = true; 

basepath = pwd;   
basename = bz_BasenameFromBasepath(basepath); 
fileinfo = dir([basename '.dat']);
[xml, ~] = LoadXml(basename); 
Fs = xml.SampleRate;
num_channels = xml.nChannels;
num_samples = fileinfo.bytes/(num_channels * 2);
rec_length = num_samples/Fs;

spiketimes_rsc = cell2mat({CellParams.SpikeTimes}');

if include_MUA
    load MUA.cellinfo.mat
    spiketimes_rsc = [spiketimes_rsc; spiketimes];
end

spiketimes_rsc = unique(spiketimes_rsc);
spiketimes_rsc = sort(spiketimes_rsc);

% in seconds 
bin_size = 0.001;
nBins = floor(rec_length/bin_size);
%% 

% kernel_packet = 0.02;
% kernel_packet = floor(kernel_packet/bin_size);
% kernel = gausswin(kernel_packet);

sigma=0.01; %standard deviation of the kernerl in s 
edges = [-3*sigma:bin_size:3*sigma]; %Time ranges from -3*std to 3*std
kernel = normpdf(edges,0,sigma); %Gaussian kernel
kernel = kernel*bin_size; %multiply by the bin width so the probabilities sum to 1

spikes_hist = hist(spiketimes_rsc,nBins);
%clear spiketimesArray

pop = filter(kernel,1,spikes_hist);
%%
win_size = 3; % in s
step_size = 0.001; 
win_size = win_size/bin_size;
step_size = step_size/bin_size; 
half_win = win_size/2; 

% SI = []; 
% 
% for nn = 1:(floor(size(test,2) - win_size)/step_size +1)
%     ini = (nn-1)*step_size+1;
%     fin = ini + win_size-1;
%     
%     temp = find(test(1,ini:fin) == 0); 
%     SI(nn) = 1 - length(temp)/win_size; 
% end

SI = zeros(1, length(pop)); 

nsteps = length(pop) - floor(win_size/2)*2; 

thresh = prctile(pop, 5);

parfor nn = 1:nsteps
    if nn == 1
        mid = nn*floor(win_size/2)+1; 
    else
        mid = (floor(win_size/2)+1) + (nn-1)*step_size;
    end
             
    temp = find(pop(1,mid-half_win:mid+half_win) <= thresh); % find(pop(1,mid-half_win:mid+half_win) == 0)
    SI(nn+half_win) = 1 - length(temp)/win_size; 
end

frag_SI = smooth(SI(half_win+1:end-half_win),500); 
frag_dat = pop(half_win+1:end-half_win); 

prc_high = 90;
prc_low = 10;
low = prctile(frag_SI, prc_low); 
high =  prctile(frag_SI, prc_high);

t2 = linspace(half_win*bin_size, rec_length - half_win*bin_size, length(frag_SI));
figure
yyaxis left
reduce_plot(t2,frag_dat);
yyaxis right
reduce_plot(t2,smooth(frag_SI,500),'r'); 
hold on
plot(get(gca,'xlim'),[low low],'k:')
plot(get(gca,'xlim'),[high high],'k:')

desyn = false([length(frag_SI) 1]); 
desyn(frag_SI>= high) = true;
desyn = [false([half_win 1]); desyn; false([half_win 1])]; 

syn = false([length(frag_SI) 1]); 
syn(frag_SI<= low) = true; 
syn = [false([half_win 1]); syn; false([half_win 1])];

%% define syn and desyn states

transitions = diff(syn);
syn_on = find(transitions==1)';  % Final output
syn_off = find(transitions==-1)';

syn_on = syn_on./1000; 
syn_off = syn_off./1000; 

transitions2 = diff(desyn);
desyn_on = find(transitions2==1)';  % Final output
desyn_off = find(transitions2==-1)';
desyn_on = desyn_on./1000; 
desyn_off = desyn_off./1000;

%% Look for RSC spikes within syn and desyn states

mua_rsc_syn = [];
for r = 1:length(syn_on)
    temp = spiketimes_rsc(spiketimes_rsc >= syn_on(r) & spiketimes_rsc <= syn_off(r)); 
    mua_rsc_syn = [mua_rsc_syn; temp]; 
end

mua_rsc_desyn = [];
for m = 1:length(desyn_on)
    temp2 = spiketimes_rsc(spiketimes_rsc >= desyn_on(m) & spiketimes_rsc <= desyn_off(m)); 
    mua_rsc_desyn = [mua_rsc_desyn; temp2]; 
end

%% Load HPC data

load('CellParams.mat')
spiketimes_hpc = cell2mat({CellParams.SpikeTimes}');


if include_MUA
    load MUA.cellinfo.mat
    spiketimes_hpc = [spiketimes_hpc; spiketimes];
end

spiketimes_hpc = unique(spiketimes_hpc);
spiketimes_hpc = sort(spiketimes_hpc);

%% Look for HPC spikes within syn and desyn states
mua_hpc_syn = []; 
for r = 1:length(syn_on)
    temp = spiketimes_hpc(spiketimes_hpc >= syn_on(r) & spiketimes_hpc <= syn_off(r)); 
    mua_hpc_syn = [mua_hpc_syn; temp]; 
end

mua_hpc_desyn = [];
for m = 1:length(desyn_on)
    temp2 = spiketimes_hpc(spiketimes_hpc >= desyn_on(m) & spiketimes_hpc <= desyn_off(m)); 
    mua_hpc_desyn = [mua_hpc_desyn; temp2]; 
end

%%

total_time = max(mua_rsc_desyn);
syn_time = sum(syn_off-syn_on);
desyn_time = sum(desyn_off - desyn_on);
rate_syn = (length(mua_rsc_syn) + length(mua_hpc_syn))/syn_time;
rate_desyn = (length(mua_rsc_desyn) + length(mua_hpc_desyn))/desyn_time;

CCG_syn = CrossCorr(mua_rsc_syn,mua_hpc_syn, 0.001, 1000); 
CCG_syn = CCG_syn./length(mua_rsc_syn)./0.001./syn_time; % / rate_syn;

CCG_desyn = CrossCorr(mua_rsc_desyn,mua_hpc_desyn, 0.001, 1000); 
CCG_desyn = CCG_desyn./length(mua_rsc_desyn)./0.001./desyn_time; % / rate_desyn;

t = linspace(-.5, .5, 1001);
figure
plot(t, smooth(CCG_syn,30), 'k')
hold on
plot(t, smooth(CCG_desyn,30), '--k')
legend('Syn','Desyn')

%% save parameters and results

params.sigma = sigma; 
params.win_size = win_size;
params.high = prc_high; 
params.low = prc_low;
params.bin_size = bin_size; 

save SI_Anaysis syn_off syn_on desyn_off desyn_on SI mua_hpc_desyn mua_hpc_syn mua_rsc_desyn mua_rsc_syn CCG_desyn CCG_syn t rate_desyn rate_syn params