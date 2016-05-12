function expRefs = LaserOmnibusCheckSessions(name,expt)
%Go through laser sessions and check they have the correct parameters
%Maybe convert this to a function to pull out unilateral and bilateral data
%on different occasions
subject = name;
expRefs = vertcat(dat.listExps(subject));

good = ones(length(expRefs),1);
for b = 1:length(expRefs)
    
    try
        p = load(dat.expFilePath(expRefs{b},'parameters','m'));
        laserGLM(expRefs{b}); %check if any lasermanip problems
        
        switch(expt)
            case 'sparse_bilateral_2D'
                if p.parameters.numRepeats(23)~=1500
                    good(b)=0;
                end
                
                if p.parameters.rewardOnStimulus(2,end)~=2.5
                    good(b)=0;
                end
                
                block = dat.loadBlock(expRefs{b});
                if block.numCompletedTrials<100
                    good(b)=0;
                end
                
                stim = p.parameters.visCueContrast;
                if ~any(min(stim,[],1)>0)
                    good(b)=0;
                end
                
            case 'sparse_unilateral_1D'
                stim = p.parameters.visCueContrast;
                if any(min(stim,[],1)>0)
                    good(b)=0;
                end
                
%                 p.parameters.rewardOnStimulus
% max(p.parameters.rewardOnStimulus(2,:))
%                 if max(p.parameters.rewardOnStimulus(2,:)) ~= 1
%                     good(b)=0;
%                 end
%                 
                l = laserGLM(expRefs{b});
                if length(l.data.response)<100
                    good(b)=0;
                end
                
                if size(l.inactivationSite,1) < 6
                    good(b)=0;
                end
                
               
                
            case 'sparse_unilateral_2D'
                stim = p.parameters.visCueContrast;                
                if ~any(min(stim,[],1)>0)
                    good(b)=0;
                end
                
                %                 p.parameters.rewardOnStimulus
                % max(p.parameters.rewardOnStimulus(2,:))
                %                 if max(p.parameters.rewardOnStimulus(2,:)) ~= 1
                %                     good(b)=0;
                %                 end
                %
                l = laserGLM(expRefs{b});
                if length(l.data.response)<100
                    good(b)=0;
                end
                
                if size(l.inactivationSite,1) < 6
                    good(b)=0;
                end
                
                 if ~any(l.data.laser(:,2)<0)
                     good(b)=0;
                 end
                 
                 if good(b)==1
                     expRefs{b}
                 end
                
            otherwise
                error('choose something');
        end
        
    catch
        good(b)=0;
    end
end

expRefs = expRefs(good==1);
end