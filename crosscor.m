basepath = pwd;
basename = bz_BasenameFromBasepath(basepath);

% define which time points will be correlated

% DSC1914_181015_1_RSC: stop/stop for rip, stop/stop for stm but robust either way
% DSC1914_181015_2_RSC: not much for rip, robust stm  either way
% DSC1914_181016_1_RSC: robust drop for rip, nothing on stm
% DSC1914_181016_2_RSC: not great, more pronounced for start/stop but also
% stop/stop for ripple, nothing on stm
% DSC4307_181016_1_RSC: not great, stop/stop for rip, not much on stm
% DSC4307_181016_2_RSC: stop/start for both but not enough data
% DSC4307_181017_1_RSC: start/start but not enough data
% DSC4307_181017_2_RSC: stop/stop good for both
% DSC4307_181018_1_RSC: start/start but stop/stop also good
% DSC4307_181018_2_RSC: not much here
% LR1_RSC_180522_a: stop/stop but also start/start big dip for rip, no stm
% LR2_RSC_180524_143335_a: stop/stop but also start/start and stop/start,
% no stm


uds_t = 'stop';
rip_t = 'start';
stm_t = 'start';

udsFil = [basepath '/' basename '.evt.uds'];
uds_evs = LoadEvents(udsFil);
uds_st = uds_evs.time(cellfun(@any,regexp(uds_evs.description,uds_t)));

ripFil = [basepath '/' basename '.evt.rip'];
rip_evs = LoadEvents(ripFil);
rip_st = rip_evs.time(cellfun(@any,regexp(rip_evs.description,rip_t)));

stmFil = [basepath '/' basename '.evt.stm'];
stm_evs = LoadEvents(stmFil);
stm_st = stm_evs.time(cellfun(@any,regexp(stm_evs.description,stm_t)));

cch1 = CrossCorr(rip_st, uds_st, 0.01, 100);
cch1 = cch1./0.01./length(rip_st);
cch2 = CrossCorr(stm_st, uds_st, 0.01, 100);
cch2 = cch2./0.01./length(stm_st);

timevec = linspace(-.5 , .5, 101);
bar(timevec,smooth(cch1))
figure
bar(timevec,smooth(cch2))

save(['/home/panteleimon/Documents/diagrams/' basename '.mat'], 'cch1','cch2')
saveas(figure(1),['/home/panteleimon/Documents/diagrams/' basename '_uds_' uds_t '_rip_' rip_t '_cch.fig'])
saveas(figure(2),['/home/panteleimon/Documents/diagrams/' basename '_uds_' uds_t '_stim_' stm_t '_cch.fig'])


% use InIntervals.m to check preferential occurence of ripples