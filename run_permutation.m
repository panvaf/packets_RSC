reps = 1000;
hits = zeros(1,reps);

for i=1:reps
    hits(i) = PacketDetectionConvolution();
    disp(i)
end

save('/home/panteleimon/Documents/packets_RSC/hits.mat','hits')