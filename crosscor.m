basepath = pwd;
basename = bz_BasenameFromBasepath(basepath);

udsFil = [basepath '/' basename '.evt.uds'];
uds_evs = LoadEvents(udsFil);
uds_st = uds_evs.time(cellfun(@any,regexp(uds_evs.description,'start')));

ripFil = [basepath '/' basename '.evt.rip'];
rip_evs = LoadEvents(ripFil);
rip_st = rip_evs.time(cellfun(@any,regexp(rip_evs.description,'start')));

stmFil = [basepath '/' basename '.evt.stm'];
stm_evs = LoadEvents(stmFil);
stm_st = stm_evs.time(cellfun(@any,regexp(stm_evs.description,'start')));

cch1 = CrossCorr(rip_st, uds_st, 0.01, 100);
cch1 = cch1./0.01./length(rip_st);
cch2 = CrossCorr(stm_st, uds_st, 0.01, 100);
cch2 = cch2./0.01./length(stm_st);

timevec = linspace(-.5 , .5, 101);
bar(timevec,smooth(cch1))
figure
bar(timevec,smooth(cch2))

save(['/home/panteleimon/Documents/diagrams/' basename '.mat'], 'cch1','cch2')
saveas(figure(1),['/home/panteleimon/Documents/diagrams/' basename '_uds_rip_cch.fig'])
saveas(figure(2),['/home/panteleimon/Documents/diagrams/' basename '_uds_stim_cch.fig'])


% use InIntervals.m to check preferential occurence of ripples