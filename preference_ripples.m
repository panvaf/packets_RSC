% check for preferential occurence of ripples in syn vs desyn states

rip_t = 'peak';

basepath = pwd;
basename = bz_BasenameFromBasepath(basepath);

synFil = [basepath '/' basename '.evt.syn'];
syn_evs = LoadEvents(synFil);
syn(:,1) = syn_evs.time(cellfun(@any,regexp(syn_evs.description,'start')));
syn(:,2) = syn_evs.time(cellfun(@any,regexp(syn_evs.description,'stop')));

desynFil = [basepath '/' basename '.evt.des'];
desyn_evs = LoadEvents(desynFil);
desyn(:,1) = desyn_evs.time(cellfun(@any,regexp(desyn_evs.description,'start')));
desyn(:,2) = desyn_evs.time(cellfun(@any,regexp(desyn_evs.description,'stop')));

ripFil = [basepath '/' basename '.evt.rip'];
rip_evs = LoadEvents(ripFil);
rip = rip_evs.time(cellfun(@any,regexp(rip_evs.description,rip_t)));

time_syn = sum(syn(:,2) - syn(:,1));
time_desyn = sum(desyn(:,2) - desyn(:,1));

[in_syn,interval,index] = InIntervals(rip,syn);
[in_desyn,interval,index] = InIntervals(rip,desyn);

tot_syn = sum(in_syn);
tot_desyn = sum(in_desyn);

r_syn = tot_syn/time_syn;
r_desyn = tot_desyn/time_desyn;
norm = r_syn + r_desyn;

p_syn = r_syn/norm
p_desyn = r_desyn/norm

% DSC1914_181015_1_RSC: 0.7165 0.2835
% DSC1914_181015_2_RSC: 0.6541 0.3459
% DSC1914_181016_1_RSC: 0.7079 0.2921
% DSC1914_181016_2_RSC: 0.7682 0.2318
% DSC1914_181017_1_RSC: 0.3002 0.6998
% DSC1914_181017_2_RSC: 0.4952 0.5018
% LR1_RSC_180522_a:     0.5755 0.4245
% LR2_RSC_180524_143335_a: 0.7854 0.2146