subjects = {'Cori','Radnitz','Muller','Moniz','Hench'};
rawDir = '\\zserver.cortexlab.net\Data\Subjects';
bin_spikes = 0.1; %Spike counting bin duration = 100msec
fr_smoothing = 0.2; %Window (in sec) to smooth firing rate estimate. Must be >= bin_spikes

for n = 1:length(subjects)
    subjDir = fullfile(rawDir,subjects{n});
    
    %Find all alf folders
    alf_dirs = dirPlus(subjDir, 'DirFilter', 'alf', 'RecurseInvalid', true, 'ReturnDirs', true);
    
    for f = 1:length(alf_dirs)
        alfDir = alf_dirs{f};
        
        parts = strsplit(alfDir,'\');
        eRef = strjoin(parts([6,5]),'_');
        
        %Load behavioural data
        feedback = readNPY(fullfile(alfDir, 'cwFeedback.type.npy'));
        choice = readNPY(fullfile(alfDir, 'cwResponse.choice.npy'));
        
        if mean(feedback==1) > 0.6 && max(choice)==3
            fprintf('GOOD: %s\n',eRef);
            contrastLeft = readNPY(fullfile(alfDir, 'cwStimOn.contrastLeft.npy'));
            contrastRight = readNPY(fullfile(alfDir, 'cwStimOn.contrastRight.npy'));
            inclTrials = readNPY(fullfile(alfDir, 'cwTrials.inclTrials.npy'));
            repNum = readNPY(fullfile(alfDir, 'cwTrials.repNum.npy'));
            cweA = table(contrastLeft, contrastRight, choice, feedback, inclTrials, repNum);
            
            % - cwtA - table of times of events in trials, containing stimOn, beeps,
            % and feedbackTime
            stimOn = readNPY(fullfile(alfDir, 'cwStimOn.times.npy'));
            beeps = readNPY(fullfile(alfDir, 'cwGoCue.times.npy'));
            responseTime = readNPY(fullfile(alfDir, 'cwResponse.times.npy'));
            cwtA = table(stimOn, beeps, responseTime);
            
            % - moveData - a struct with moveOnsets, moveOffsets, moveType
%             wheelMoves = readNPY(fullfile(alfDir, 'wheelMoves.intervals.npy'));
%             moveOnsets = wheelMoves(:,1); moveOffsets = wheelMoves(:,2);
%             moveType = readNPY(fullfile(alfDir, 'wheelMoves.type.npy'));
%             moveData.moveOnsets = moveOnsets; moveData.moveOffsets = moveOffsets; moveData.moveType = moveType;
            wheelPos = readNPY(fullfile(alfDir,'wheel.position.npy'));
            wheelT = readNPY(fullfile(alfDir,'wheel.timestamps.npy'));
            wheelT = linspace(wheelT(1,2),wheelT(2,2), length(wheelPos) );
            
            %Create own table
            behav = table;
            behav.stimulus = [cweA.contrastLeft cweA.contrastRight];
            behav.response = categorical(cweA.choice, 1:3, {'Left','Right','NoGo'});
            behav.repeatNum = cweA.repNum;
            behav.feedbackType = categorical(cweA.feedback, [-1 1], {'Incorrect','Correct'});
            behav.stimulusOnTime = cwtA.stimOn;
            behav.goCueTime = cwtA.beeps;
            behav.responseTime = cwtA.responseTime;
            behav = behav(cweA.inclTrials==1,:);
            
            %Resample the response Times for GO trials to get an estimate
            %for the time of the actual action to NoGo
            goCue2ResponseTime = behav.responseTime - behav.goCueTime;
            goCue2ResponseTime(behav.response=='NoGo') = [];
            goCue2NoGoResponseTime = randsample(goCue2ResponseTime,sum(behav.response=='NoGo'));
            behav.responseTime(behav.response=='NoGo') = behav.goCueTime(behav.response=='NoGo') + goCue2NoGoResponseTime;
            
            
            %Detect early GO movements (after stim onset but before go cue)
            firstMoveTime = nan(size(behav.response));
            figure('color','w'); subplot(2,1,1); hold on;
            
            RT_estimate = nan(size(behav.response));
            for tr = 1:height(behav)
                %                 idx = find(moveData.moveOnsets > behav.stimulusOnTime(tr),1,'first');
                if behav.response(tr) ~= 'NoGo'
                    idx = behav.stimulusOnTime(tr) < wheelT & wheelT < behav.stimulusOnTime(tr)+1;
                    pos1 = wheelPos(idx);
                    t1 = wheelT(idx);
                    
                    %Transform to common scaling
                    pos1 = pos1 - pos1(1);
                    pos1_norm = pos1/max(abs(pos1));
                    
                    hxx=plot(t1-t1(1),pos1);
                    
                    if behav.response(tr) == 'Left'
                        idx = find(pos1_norm > +0.2,1,'first');   
                        hxx.Color = [0.4940 0.1840 0.5560 0.5];
                    elseif behav.response(tr) == 'Right'
                        idx = find(pos1_norm < -0.2,1,'first');
                        hxx.Color = [0.4660 0.6740 0.1880 0.5];
                    end
                    

                    if isempty(idx)
                        firstMoveTime(tr,1) = behav.responseTime(tr);
                    else
                        firstMoveTime(tr,1) = t1(idx);
                    end
                    
                    RT_estimate(tr,1) = firstMoveTime(tr,1) - behav.stimulusOnTime(tr,1);
                end
            end
            set(gca,'xcolor','w');
            title('Rotary encoder aligned to stim onset');
            
            %Resample the firstMove time for GO trials, for NoGo trials
            stimOn2firstMove = firstMoveTime - behav.stimulusOnTime;
            stimOn2firstMove(isnan(stimOn2firstMove)) = [];
            stimOn2firstMoveNoGo = randsample(stimOn2firstMove,sum(behav.response=='NoGo'));
            firstMoveTime(behav.response=='NoGo') = behav.stimulusOnTime(behav.response=='NoGo') + stimOn2firstMoveNoGo;
            
            subplot(2,1,2);
            distributionPlot(RT_estimate(behav.response=='Left'),'widthDiv',[2 2],'color',[0.4940 0.1840 0.5560],'histOri','right','showMM',0,'xyOri','flipped');
            distributionPlot(RT_estimate(behav.response=='Right'),'widthDiv',[2 1],'color',[0.4660 0.6740 0.1880],'histOri','left','showMM',0,'xyOri','flipped');
            set(gca,'ycolor','w','xcolor','k'); title('Time histogram: stim onset -> first mov');

            linkaxes(get(gcf,'children'),'x');
            xlim([0 1]);
            xlabel('Time from stimulus onset');
            
            set(gcf,'Position',[0 0 600 1000]);
            print(gcf,['../figures/chronometry/' eRef],'-dpdf','-bestfit');
            

            behav.firstMoveTime = firstMoveTime;
            epoches = {'stimulusOnTime','goCueTime','firstMoveTime','responseTime'};
            
            %Fit GLM and plot
            D = struct;
            D.stimulus = behav.stimulus;
            D.response = double(behav.response);
            g = GLM(D).setModel('C50-subset').fit;
            fig = g.plotFit;
            set(fig,'Position',[0 0 600 1000]);
            print(fig,['../figures/behavModel/' eRef],'-dpdf','-bestfit');
            close all;
            
            ZL = g.ZL(g.parameterFits, g.Zinput(g.data));
            ZR = g.ZR(g.parameterFits, g.Zinput(g.data));
            behav.GLM_offset = [ZL ZR zeros(size(ZL))];
            
            % Now find ephys data
            ephys_dat_folders = dirPlus(alfDir,'ReturnDirs', true);

            spikes = table;
            population = table;
            firingRate_times = (behav.stimulusOnTime(1)-5) : bin_spikes : (behav.responseTime(end)+5);
            for e = 1:length(ephys_dat_folders)
                ephysDir = ephys_dat_folders{e};
                
                parts = strsplit(ephysDir,'\');
                ephysTag = parts{end};
                
                % - st - vector of spike times
                spikeTimes = readNPY(fullfile(ephysDir, 'spikes.times.npy'));
                %Trim away spike Times beyond the time period of the experiment
                exptIdx = firingRate_times(1) < spikeTimes & spikeTimes < firingRate_times(end);
                spikeTimes = spikeTimes(exptIdx);
                
                % - clu - vector of cluster identities
                spikeCluster = readNPY(fullfile(ephysDir, 'spikes.clusters.npy'));
                spikeCluster = spikeCluster(exptIdx);
                
                sp=table;
                clusterIDs = readNPY(fullfile(ephysDir, 'clusters.ids.npy'));
                clusterType = readNPY(fullfile(ephysDir, 'clusters.groups.npy'));
                sp.clusterType = categorical(clusterType,[0 1 2 3],{'Noise','MUA','Good','Unsorted'});
                clusterDepth = readNPY(fullfile(ephysDir, 'clusters.depths.npy'));
                 
                %For each cluster, generate estimate of firing rate by binning the spike
                %events
                sp.firingRates = nan(length(clusterIDs), length(firingRate_times));
                for clu = 1:length(clusterIDs)
                    cluID = clusterIDs(clu);
                    clu_spikeTimes = spikeTimes(spikeCluster==cluID);
                    
                    spikeCounts = hist(clu_spikeTimes, firingRate_times);
                    sp.firingRates(clu,:) = spikeCounts/bin_spikes;
                end
                
                %Mark the brain structure each cluster is within
                bordersFile = fullfile(alfDir, ephysTag, ['borders_' ephysTag '.tsv']);
                borders = readtable(bordersFile ,'Delimiter','\t', 'FileType', 'text');
                
                sp.clusterRegion = cell(size(clusterIDs));
                for area = 1:height(borders)
                    between = borders.lowerBorder(area) <= clusterDepth & clusterDepth <= borders.upperBorder(area);
                    sp.clusterRegion(between) = borders.name(area);
                    
                    %Also generate population rate for that region
                    pop = table;
                    pop.region = borders.name(area);
                    popSpikes = spikeTimes( any(spikeCluster == clusterIDs(between)',2) );
                    spikeCounts = hist(popSpikes, firingRate_times);
                    pop.firingRates = spikeCounts/bin_spikes;
                    population = [population; pop];
                end
                
                %Remove any clusters with no identified brain region
                noArea = find(cellfun(@(s)isempty(s),sp.clusterRegion));
                sp(noArea,:) = [];

                spikes = [spikes; sp];
            end
            
            spikes.Properties.UserData = struct('fr_times',firingRate_times,'bin_spikes',bin_spikes);
            population.Properties.UserData = struct('fr_times',firingRate_times,'bin_spikes',bin_spikes);
            save(['../preproc/' eRef '.mat'],'behav','spikes','population','epoches');
        else
            fprintf('BAD: %s\n',eRef);
        end
        
    end
    
end