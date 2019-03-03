% compute packet statistics

load StateIndex.mat

starts = UDS(1,:);
ends = UDS(2,:);
dur = ends - starts;
IPI = starts(2:end) - ends(1:end-1);
    
durbef_IPI = corr(IPI',dur(1:end-1)');
dur_IPIbef = corr(IPI',dur(2:end)');

cat = {'dur-IPI';'IPI-dur'};
val = [durbef_IPI dur_IPIbef];
barh(val)
set(gca,'yticklabel',cat)
title('Correlations between packet strength and Inter-Packet-Interval')

figure
histogram(dur,10,'Normalization','pdf')
title('Duration of packets')
xlabel('Duration (s)')