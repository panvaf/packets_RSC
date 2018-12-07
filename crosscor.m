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



cch = CrossCorr(rip_st, uds_st, 0.01, 100);
cch = CrossCorr(stm_st, uds_st, 0.01, 100);


