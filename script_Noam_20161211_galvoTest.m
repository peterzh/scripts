
% script to load data for 


%% file locations

mouseName = 'Noam';
thisDate = '2016-12-11';
ephysTag = 'V1';
tlExpNum = 1;
galvoTestExpNum = 4;
masterTimebase = 'SC';

ksDir = fullfile('\\basket.cortexlab.net\data\nick\', mouseName, thisDate, ['ephys_' ephysTag]);

expRoot = fileparts(dat.expPath(mouseName, thisDate, 1, 'main', 'master'));
rawDir = fullfile(expRoot, ['ephys_' ephysTag]);

alignDir = fullfile(expRoot, 'alignments');

tlFile = dat.expFilePath(mouseName, thisDate, tlExpNum, 'timeline', 'master');
protocolFile = dat.expFilePath(mouseName, thisDate, galvoTestExpNum, 'parameters', 'master');

%% sync information

if ~strcmp(ephysTag, masterTimebase)
    bEphysToMaster = readNPY(fullfile(alignDir, ...
        sprintf('correct_ephys_%s_to_ephys_%s.npy', ephysTag, masterTimebase)));
else % this one is master, so use a dummy conversion
    bEphysToMaster = [1; 0];
end

bTLtoMaster = readNPY(fullfile(alignDir, ...
    sprintf('correct_timeline_%d_to_ephys_%s.npy', tlExpNum, masterTimebase)));

stimOnTL = readNPY(fullfile(alignDir, ...
    sprintf('mpep_%d_onsets_in_timeline_%d.npy', galvoTestExpNum, tlExpNum)));
    
stimOffTL = readNPY(fullfile(alignDir, ...
    sprintf('mpep_%d_offsets_in_timeline_%d.npy', galvoTestExpNum, tlExpNum)));

stimOn = applyCorrection(stimOnTL, bTLtoMaster);
stimOff = applyCorrection(stimOffTL, bTLtoMaster);


%% load protocol information

load(protocolFile);
Protocol = parameters.Protocol;

amplitudes = Protocol.pars(strcmp(Protocol.parnames, 'amp1'),:);
shape = Protocol.pars(find(strcmp(Protocol.parnames, 'shape'),1,'last'),:);

stimIDs = zeros(1, prod(size(Protocol.seqnums)));
for q = 1:size(Protocol.seqnums,1)
    stimIDs(Protocol.seqnums(q,:)) = q;
end

ampsByTrial = amplitudes(stimIDs);
shapesByTrial = shape(stimIDs);

%% detect onsets of laser directly from timeline
load(tlFile);
tt = Timeline.rawDAQTimestamps;
wave = Timeline.rawDAQData(:, strcmp({Timeline.hw.inputs.name}, 'waveOutput'));

thresh = [0.05 0.1]; 
[~, waveOnsets] = schmittTimes(tt, wave, thresh);

waveOnsets = waveOnsets(waveOnsets>3100 & waveOnsets<3500); % this is the experiment we're interested in
waveOnsets = waveOnsets(diff([0; waveOnsets])>0.2); % only keep the first onset in each trial

% convert times of wave onsets to ephys times
waveOnsetsEphys = applyCorrection(waveOnsets, bTLtoMaster);

%% see that onsets detected correctly, and labels correct

figure; plot(tt, wave)
hold on; plot(waveOnsets, zeros(size(waveOnsets)), 'r.')
for a = unique(ampsByTrial)
    if a>0        
        galvoOnsets = ampsByTrial(ampsByTrial>0)==a & shapesByTrial(ampsByTrial>0)==2;
        sineOnsets = ampsByTrial(ampsByTrial>0)==a & shapesByTrial(ampsByTrial>0)==1;
        plot(waveOnsets(galvoOnsets), a/1000, 'go');
        plot(waveOnsets(sineOnsets), a/1000, 'ko');
    end
end
xlim([3100 3500]);

%% load spikes

s = loadKSdir(ksDir);

s.st = applyCorrection(s.st, bEphysToMaster);

%% compute basic things from spikes
[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = ...
    templatePositionsAmplitudes(s.temps, s.winv, s.ycoords, s.spikeTemplates, s.tempScalingAmps);

spikeDurs = templateDuration(s.spikeTemplates+1);

%% make some psth's
inclSpikes = spikeAmps>50; 
cluByDepth = ceil(spikeDepths/80)*80; 


% trGroups = ampsByTrial(ampsByTrial>0);
% trGroups(shapesByTrial(ampsByTrial>0)==2) = trGroups(shapesByTrial(ampsByTrial>0)==2)+8000; 

trGroups = ampsByTrial;
trGroups(shapesByTrial==2 & ampsByTrial>0) = trGroups(shapesByTrial==2 & ampsByTrial>0)+8000; 

narrowSpikes = spikeDurs<0.4e-3*s.sample_rate; % spikes less than 0.4ms duration

% psthViewer(s.st, cluByDepth, stimOn, [-1 2], trGroups);
psthViewer(s.st(~narrowSpikes), cluByDepth(~narrowSpikes), stimOn, [-0.25 1.25], trGroups);
% psthViewer(s.st(~narrowSpikes), cluByDepth(~narrowSpikes), stimOn, [-0.25 1.25], trGroups);


%% psth across depth

% inclTrials = ampsByTrial>0 & shapesByTrial==2;
% inclTrials = ampsByTrial>0 & shapesByTrial==2;
inclTrials = ampsByTrial==3000 & shapesByTrial==2;

eventTimes = stimOn(inclTrials); eventName = 'stimOnset';
inclSpikes = spikeAmps>50 & ~narrowSpikes; % narrow
% inclSpikes = spikeAmps>50 & ~narrowSpikes; % broad

win = [-0.25 1.25];
timeBinSize = 0.005;

bslWin = [-0.25 -0.1];
depthBinSize = 100;

f = figure; 

[timeBins, depthBins, allP] = psthByDepth(s.st(inclSpikes), ...
    spikeDepths(inclSpikes), depthBinSize, timeBinSize, eventTimes, win, bslWin);
plotPSTHbyDepth(timeBins, depthBins, allP, eventName, 'norm', [])
caxis([-1 1]*30);
% [timeBins, depthBins, allP] = psthByDepth(s.st(inclSpikes), ...
%     spikeDepths(inclSpikes), depthBinSize, timeBinSize, eventTimes, win, []);
% plotPSTHbyDepth(timeBins, depthBins, allP, eventName, [], [])

%% laser delta by depth and laser amplitude

amps = unique(ampsByTrial); amps(1)=[];
shapes = unique(shapesByTrial);

depthBinSize = 100;
bslWin = [-0.25 -0.1];

shapeTitles = {'Sine wave','Galvo emulation'};
figure('color','w');
for sh = shapes
    
    subplot(1,2,sh); hold on;
    
    for IE = [0 1]
        for am = amps
            
            inclTrials = (ampsByTrial==am) & shapesByTrial==sh;
            eventTimes = stimOn(inclTrials);
            inclSpikes = spikeAmps>50 & (narrowSpikes==IE);
            
            [timeBins, depthBins, allP] = psthByDepth(s.st(inclSpikes), spikeDepths(inclSpikes), depthBinSize, timeBinSize, eventTimes, win, bslWin);
            ax=plot(mean(allP(:,timeBins>0 & timeBins<0.15),2),depthBins(1:end-1),'.:');
            ax.MarkerSize=20;
            
            col = 1 - am/max(amps);
            if IE==0
                ax.Color=[1 col col];
            else
                ax.Color=[col col 1];
            end
        end
        
    end
    
    xlabel('Post-pre laser psth [z score]'); ylabel('depth on electrode array (µm)');
    title(shapeTitles{sh});
end





