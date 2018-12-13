load('CellParams.mat')
%%
basepath = pwd;
basename = bz_BasenameFromBasepath(basepath);
synFil = [basepath '/' basename '.evt.syn'];
syn_evs = LoadEvents(synFil);
syn_start = syn_evs.time(cellfun(@any,regexp(syn_evs.description,'start')));
syn_stop = syn_evs.time(cellfun(@any,regexp(syn_evs.description,'stop')));

spiketimes_rsc = cell2mat({CellParams.SpikeTimes}');

mua_rsc_syn = []; 
mua_rsc_desyn = [];

for r = 1:length(syn_start)
    temp = spiketimes_rsc(spiketimes_rsc >= syn_start(r) & spiketimes_rsc <= syn_stop(r)); 
    mua_rsc_syn = [mua_rsc_syn; temp]; 
end


mua_rsc_desyn = setdiff(spiketimes_rsc, mua_rsc_syn); 


%% 
load('CellParams.mat')

% ShankID = cell2mat({CellParams.ShankID}');
% hp = find(ShankID == 3) ; 
% hp = hp(1)-1; 
% hp = find(ShankID == 4) ; 
% hp = hp(end)+1;  

%spiketimes_hpc = {CellParams(1:hp).SpikeTimes}';
% spiketimes_hpc = {CellParams(hp:end).SpikeTimes}';
spiketimes_hpc = cell2mat({CellParams.SpikeTimes}');

% mua_hpc = [];

% for h = 1:length(spiketimes_hpc)
%    mua_hpc = [mua_hpc; spiketimes_hpc{h}]; 
% end

mua_hpc_syn = []; 
mua_hpc_desyn = [] ;

for r = 1:length(syn_start)
    temp = spiketimes_hpc(spiketimes_hpc >= syn_start(r) & spiketimes_hpc <= syn_stop(r)); 
    mua_hpc_syn = [mua_hpc_syn; temp]; 
end


mua_hpc_desyn = setdiff(spiketimes_hpc, mua_hpc_syn); 

%%
total_time = max(mua_rsc_desyn);
syn_time = sum(syn_stop-syn_start);
desyn_time = total_time - syn_time;
rate_syn = (length(mua_rsc_syn) + length(mua_hpc_syn))/syn_time;
rate_desyn = (length(mua_rsc_desyn) + length(mua_hpc_desyn))/desyn_time;

CCG_syn = CrossCorr(mua_hpc_syn, mua_rsc_syn, 0.001, 10000); 
CCG_syn = CCG_syn./length(mua_hpc_syn)./0.001./rate_syn; 

CCG_desyn = CrossCorr(mua_hpc_desyn, mua_rsc_desyn, 0.001, 10000); 
CCG_desyn = CCG_desyn./length(mua_hpc_desyn)./0.001/rate_desyn; 

t = linspace(-5, 5, 10001);
figure(5), clf
plot(t, smooth(CCG_syn,30), 'k')
hold on
plot(t, smooth(CCG_desyn,30), '--k')
xlim([-0.1 0.1])
legend('Syn','Desyn')