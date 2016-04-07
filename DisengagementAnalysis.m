%Script to try to uncover bouts of task disengagement in the mice. One
%approach is to plot the reaction times for each trial, separated by the
%contrast types. Regions of disengagement should have larger RT. Similarly
%one can plot pupil diameters and see whether they predict disengagement
%well. A more model-based approach is to fit a GLM to a session, and then
%plot the predictive performance of the model on each trial, separatly
%plotting for each choice made. The idea is that in periods of
%disengagement, the GLM parameter fits do not reflect the internal driver
%of the behaviour and therefore during these trials the model will do a bad
%job. i.e. the neg loglik will be high on those trials.

expRefs = {'2015-12-03_1_EJ007';
    '2015-10-13_2_EJ007';
    '2015-12-01_2_EJ009'};

expRefs = dat.listExps('Laveran');
expRefs = expRefs(98:end);
% % % 
% expRefs = {'2013-07-03_2_AR062';
% '2013-07-04_1_AR062';
% '2013-07-04_2_AR062';
% '2013-07-05_1_AR062';
% '2013-07-03_5_AR061';
% '2013-07-04_1_AR061';
% '2013-07-05_1_AR061';
% '2013-07-05_2_AR061'};

for s = 1:length(expRefs)
    try
        D = struct;
        
        try
            block = dat.loadBlock(expRefs{s});
            trials = block.trial;
            
            for t=1:block.numCompletedTrials
                D.contrast_cond(t,:) = trials(t).condition.visCueContrast';
                D.response(t,1) = trials(t).responseMadeID';
                D.repeatNum(t,1) = trials(t).condition.repeatNum;
                D.feedbackType(t,1) = trials(t).feedbackType;
                D.RT(t,1) = trials(t).responseMadeTime - trials(t).interactiveStartedTime;
                D.trialTimes(t,1) = trials(t).trialStartedTime;
            end
        catch %For new signals data format
            block = dat.loadBlock(expRefs{s});
            
            for t=1:length(block.events.responseValues)
                c = block.paramsValues(t).contrast*sign(block.paramsValues(t).targetAzimuth);
                if c==0
                    D.contrast_cond(t,:) = [0 0];
                elseif c<0
                    D.contrast_cond(t,:) = [-c 0];
                elseif c>0
                    D.contrast_cond(t,:) = [0 c];
                end
                
                r=block.events.responseValues(t);
                if r==-1
                    D.response(t,1) = 1;
                elseif r==0
                    D.response(t,1) = 3;
                elseif r==1
                    D.response(t,1) = 2;
                end
                
                D.repeatNum(t,1) = block.events.repeatNumValues(t);
                D.feedbackType(t,1) = block.events.feedbackValues(t)==1;
                D.RT(t,1) = block.events.responseTimes(t) - block.events.interactiveOnTimes(t);
                D.trialTimes(t,1) = block.events.trialNumTimes(t);
            end
        end
        
        g = GLM(getrow(D,D.repeatNum==1));
        f=figure;
        
        if g.ContrastDimensions==1
            g=g.setModel('C50-subset').fit;
            subplot(2,3,3);
            g.plotFit;
        else
            g=g.setModel('C50').fit;
        end
        
        ph = g.calculatePhat(g.parameterFits,g.Zinput(D)); %Test model on all data
        loglik=log2(ph(sub2ind(size(ph), [1:length(D.response)]', D.response)));
        
        
        h(1)=subplot(2,3,[1 2]); hold on;
        for r = 1:max(D.response)
            trials = find(D.response==r);
            plot(D.trialTimes(trials),-loglik(D.response==r),'.','markersize',10);
        end
        hold off; xlabel('trial start time'); ylabel('negative log lik');
        legend({'Chose L','Chose R','NG'}); title(expRefs{s},'interpreter','none');
        
        h(2)=subplot(2,3,[4 5]); hold on;
        for r = 1:max(D.response)
            trials = find(D.response==r);
            medianRT = median(D.RT(D.response==r));
            plot(D.trialTimes(trials),D.RT(D.response==r)-medianRT,'.','markersize',10);
        end
        hold off; xlabel('trial start time'); ylabel('RT - median RT');
        
        linkaxes(h,'x');
        
        savefig(f,['B:\figures\Disengagement\' expRefs{s} '.fig']);
        close all;
    catch
        warning(['error for' expRefs{s}]);
    end
end
