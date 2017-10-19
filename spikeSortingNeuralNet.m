%% Train a neural net to help with manual spike sorting

%Load some data
ksDir = '\\zserver.cortexlab.net\Data\Subjects\Hench\2017-06-18\ephys_K3\sorting';
params.excludeNoise = false;
s = loadKSdir(ksDir,params);

cluIDs = s.cids';
Y = s.cgs'; %The classification into noise, mua or good

error('unfinished');