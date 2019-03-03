% code written by Sam Mckenzie, subsequent modifications by panvaf

%% load data
load DSC1914_181015_1_RSC.spikes.cellinfo.mat
load StateIndex.mat

%% bin and smooth
mn = min(cellfun(@min,spikes.times));
mx = max(cellfun(@max,spikes.times));
% k = gaussian2Dfilter([100 1],5);
k = fspecial('gaussian',[100 1],5);
dt = .001;
ts = mn-dt:dt:mx+dt;
temp = cell2mat(cellfun(@(a) nanconvn(histc(a,ts),k)',spikes.times,'uni',0)');

%% only keep packets with a duratiom in specified interval
NEURAL = temp;
%{
min_dur = .05; % in sec
max_dur = .2;

for i=1:size(UDS,2)
    dif = UDS(2,i) - UDS(1,i);
    if dif < min_dur | dif > max_dur
        UDS(2,i) = 0;
        UDS(1,i) = 0;
    end
end

UDS(:,UDS(1,:)==0) = [];


%% choose time bins when packets have been detected
index = round((UDS-mn)/dt);
tot = sum(index(2,:)-index(1,:));
NEURAL = zeros(size(temp,1),tot);

ind = 1;
for i=1:size(index,2)
    len = index(2,i)-index(1,i);
    NEURAL(:,ind:ind+len) = temp(:,index(1,i):index(2,i));
    ind = ind + len;
end
%}
%% break data into training set and test set
splitN = floor(size(NEURAL,2)*.1);

trainNEURAL = NEURAL(:,1:splitN); 

testNEURAL = NEURAL(:,(splitN+1):end); 
%% plot one example factorization
rng(235); % fixed rng seed for reproduceability
X = trainNEURAL;
K = 10;
L = .1; % units of seconds
Lneural = ceil(L/dt);

shg
display('Running seqNMF on real neural data (from songbird HVC, recorded by Emily Mackevicius, Fee Lab)')
[W, H, ~,loadings,power]= seqNMF(X,'K',K,'L',Lneural,...
            'lambdaL1W', .1, 'lambda', .0001, 'maxiter', 100, 'showPlot', 1,...
            'lambdaOrthoW', 0); 
p = .05; % desired p value for factors


%%

events = LoadEvents('LR1_RSC_180522_a.evt.stm');
stim = events.time(cellfun(@any,regexp(events.description,'start')));
events = LoadEvents('LR1_HPC_180522_a.evt.rip');
rip = events.time(cellfun(@any,regexp(events.description,'start')));
  %%
  
 
[N,K,L] = size(W);
idx = repmat(b,1,201) + repmat(-100:100,length(b),1);
T = size(NEURAL,2);
WTX = zeros(K, T);

for l = 1 : L
    %X_shifted = circshift(X,-l+1,2);       
    X_shifted = circshift(NEURAL,[0,-l+1]);       
    WTX = WTX + W(:, :, l)' * X_shifted;
end   

%%

stim = stim(stim>ts(1));
[~,b]  =histc(stim,ts);
  
idx = repmat(b,1,4001)-50 + repmat(-2000:2000,length(b),1);
idx = idx(all(idx,2)>0,:);
stim = stim(stim>ts(1));
[~,b]  =histc(stim,ts);proj_stim= [];

for i = 1:5
   p = WTX(i,:);
   proj_stim{i} = p(idx);
end


rip = rip(rip>ts(1));
[~,b]  =histc(rip,ts);
  
idx = repmat(b,1,4001)-50 + repmat(-2000:2000,length(b),1);
idx = idx(all(idx,2)>0,:);
rip = rip(rip>ts(1));
[~,b]  =histc(rip,ts);proj_rip= [];

for i = 1:5
   p = WTX(i,:);
   proj_rip{i} = p(idx);
end

%%
figure
ax = tight_subplot(5,2);
ax = reshape(ax,2,5);
ax  = ax';
for i = 1:10
axes(ax(i,1))
plot(-2000:2000,nanmean(proj_stim{i}),'k')
xlim([-500 500])
axes(ax(i,2))

plot(-2000:2000,nanmean(proj_rip{i}),'r')
xlim([-500 500])
end
%%
[~,b]  =histc(stim,ts);
idx = repmat(b,1,1001) + repmat(-500:500,length(b),1);
idx = idx(all(idx,2)>0,:);
[a,b] = pastalkova(squeeze(W(:,1,:)));
s = nan(size(NEURAL,1),size(idx,2),size(idx,1));
for i = 1:size(NEURAL,1)
    x = NEURAL(i,:);
    s(i,:,:) = x(idx');
end

%{

%%

%%%%
%I HAVEN'T PLAYED WITH THIS YET -SM



[pvals,is_significant] = test_significance(testNEURAL,W,p);

W = W(:,is_significant,:); 
H = H(is_significant,:); 

% plot, sorting neurons by latency within each factor
[max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(W(:,:,:),1);
indSort = hybrid(:,3);
tstart = 180; % plot data starting at this timebin
figure; WHPlot(W(indSort,:,:),H(:,tstart:end), X(indSort,tstart:end), ...
    0,trainSONG(:,floor(tstart*SONGfs/VIDEOfs):end))
title('Significant seqNMF factors, with raw data')
figure; WHPlot(W(indSort,:,:),H(:,tstart:end), ...
    helper.reconstruct(W(indSort,:,:),H(:,tstart:end)),...
    0,trainSONG(:,floor(tstart*SONGfs/VIDEOfs):end))
title('SeqNMF reconstruction')

%% Procedure for choosing lambda
nLambdas = 20; % increase if you're patient
K = 10; 
X = trainNEURAL;
lambdas = sort([logspace(-1,-5,nLambdas)], 'ascend'); 
loadings = [];
regularization = []; 
cost = []; 
for li = 1:length(lambdas)
    [N,T] = size(X);
    [W, H, ~,loadings(li,:),power]= seqNMF(X,'K',K,'L',Lneural,...
        'lambdaL1W', .1, 'lambda', lambdas(li), 'maxiter', 100, 'showPlot', 0); 
    [cost(li),regularization(li),~] = helper.get_seqNMF_cost(X,W,H);
    display(['Testing lambda ' num2str(li) '/' num2str(length(lambdas))])
end
%% plot costs as a function of lambda
windowSize = 3; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
Rs = filtfilt(b,a,regularization); 
minRs = prctile(regularization,10); maxRs= prctile(regularization,90);
Rs = (Rs-minRs)/(maxRs-minRs); 
R = (regularization-minRs)/(maxRs-minRs); 
Cs = filtfilt(b,a,cost); 
minCs =  prctile(cost,10); maxCs =  prctile(cost,90); 
Cs = (Cs -minCs)/(maxCs-minCs); 
C = (cost -minCs)/(maxCs-minCs); 

clf; hold on
plot(lambdas,Rs, 'b')
plot(lambdas,Cs,'r')
scatter(lambdas, R, 'b', 'markerfacecolor', 'flat');
scatter(lambdas, C, 'r', 'markerfacecolor', 'flat');
xlabel('Lambda'); ylabel('Cost (au)')
set(legend('Correlation cost', 'Reconstruction cost'), 'Box', 'on')
set(gca, 'xscale', 'log', 'ytick', [], 'color', 'none')
set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])

%% chose lambda=.005; run multiple times, see number of sig factors
loadings = [];
pvals = []; 
is_significant = []; 
X = trainNEURAL;
nIter = 20; % increase if patient
display('Running seqNMF multiple times for lambda=0.005')

for iteri = 1:nIter
    [W, H, ~,loadings(iteri,:),power]= seqNMF(X,'K',K,'L',Lneural,...
            'lambdaL1W', .1, 'lambda', .005, 'maxiter', 100, 'showPlot', 0); 
    p = .05;
    [pvals(iteri,:),is_significant(iteri,:)] = test_significance(testNEURAL,W,p);
    W = W(:,is_significant(iteri,:)==1,:); 
    H = H(is_significant(iteri,:)==1,:); 
    [max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(W(:,:,:),1);
    indSort = hybrid(:,3);
    tstart = 300; 
    clf; WHPlot(W(indSort,:,:),H(:,tstart:end), X(indSort,tstart:end), 0,trainSONG(:,floor(tstart*SONGfs/VIDEOfs):end))
    display(['seqNMF run ' num2str(iteri) '/' num2str(nIter)])
end
figure; hold on
h = histogram(sum(is_significant,2), 'edgecolor', 'w', 'facecolor', .7*[1 1 1]); 
h.BinCounts = h.BinCounts/sum(h.BinCounts)*100; 
xlim([0 10]); 
xlabel('# significant factors')
ylabel('% seqNMF runs')

%% Plot factor-triggered song examples and rastors
addpath(genpath('misc_elm')); 
figure; HTriggeredSpec(H,trainSONG,VIDEOfs,SONGfs,Lsong); 

figure; HTriggeredRaster(H,trainNEURAL(indSort,:),Lneural);
%}
