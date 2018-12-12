intervals = false;
basepath = pwd;   
basename = bz_BasenameFromBasepath(basepath); 
ripFil = [basepath '/' basename '.evt.rip'];
rip_evs = LoadEvents(ripFil);

if intervals
rip_st = rip_evs.time(cellfun(@any,regexp(rip_evs.description,'peak')));
[status,interval,index] = InIntervals(rip_st, SYN');
syn_perc = sum(syn)/length(syn);
inside_perc =  sum(status)/length(status);
pref = inside_perc/syn_perc

% 1.64, 1.14, 1.73, 1.59, 1.11, 1.14, 1.22, 1.21, 1.35, 1.20, 1.07, 3.49,
% 1.4, 2.02

else
    
%load StateIndex.mat
rip_st(:,1) = rip_evs.time(cellfun(@any,regexp(rip_evs.description,'start')));
rip_st(:,2) = rip_evs.time(cellfun(@any,regexp(rip_evs.description,'stop')));
rip = zeros(1,length(state_index));
for i=1:size(rip_st,1)
    rip(round(rip_st(i,1)*Fs):round(rip_st(i,2)*Fs)) = 1;
end
max_lag = 30; % in s
N = round(max_lag*Fs);
r = xcorr(state_index,rip,N,'coeff');
time = linspace(-10,10,length(r));
plot(time,r)
end