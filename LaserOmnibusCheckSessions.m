%Go through laser sessions and check they have the correct parameters
%Maybe convert this to a function to pull out unilateral and bilateral data
%on different occasions
subject = 'Murphy';
expRefs = vertcat(dat.listExps(subject));

good = ones(length(expRefs),1);
for b = 1:length(expRefs)
    
    try
        p = load(dat.expFilePath(expRefs{b},'parameters','m'));
        
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
%         trials = block.trial;
    
    catch
        good(b)=0;
    end
end

expRefs(good==1)