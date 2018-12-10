% Difference in spike count between KS1 and KS2
ks1=loadKSdir('E:\',struct('excludeNoise',false));
ks2=loadKSdir('D:\ks2',struct('excludeNoise',false));

cluid_ks1 = 960;
cluid_ks2 = 151;

%Compare spike times
stId_ks1 = find(ks1.clu==cluid_ks1);
stId_ks2 = find(ks2.clu==cluid_ks2);
st_ks1 = ks1.st(stId_ks1);
st_ks2 = ks2.st(stId_ks2);

[n1,edges] = histcounts( st_ks1, 10000);
n2 = histcounts( st_ks2, edges);

% compare number of spikes over time
figure;
plot(n1); hold on;
plot(n2);

%are those extra spikes noise? find them and assign to a new cluster to
%visualise in phy
KS2_sp_without_match = [];
for ks2sp = 1:length(st_ks2) %For each spike in KS2
    thisSpikeID = stId_ks2(ks2sp);
    
    %Was this KS2 spike identified in KS1 for this cluster?
    if ~any( abs(st_ks2(ks2sp) - st_ks1) < 2/1000 ) %If not any KS1 spikes within 1ms of this KS2 spike
        KS2_sp_without_match = [KS2_sp_without_match; thisSpikeID];
    end    
end

% 
% Read spikeTimes and overwrite with reassigned KS2 clusters
clu = readNPY(fullfile('D:\ks2', 'spike_clusters.npy'));
movefile(fullfile('D:\ks2', 'spike_clusters.npy'),fullfile('D:\ks2', 'spike_clusters_original.npy'));
clu(KS2_sp_without_match) = 1000;
writeNPY(clu, fullfile('D:\ks2', 'spike_clusters.npy'));


%Identify the extra spike times which KS2 detected, and see whether they are
%listed in ANY other nearby clusters for KS1
