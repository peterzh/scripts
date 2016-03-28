%Load the file
%%load(..)

D=struct;
numSessions = length(NeuronsFR_DB_EJ008_Batch3);
for s = 1:numSessions
    numTrials = size(NeuronsFR_DB_EJ008_Batch3(s).data,1);
    stim = NeuronsFR_DB_EJ008_Batch3(s).data(:,2);
    stim = sort([stim zeros(numTrials,1)],2);
    stim = abs(stim);
    D.stimulus = stim;
    
    D.action = (NeuronsFR_DB_EJ008_Batch3(s).data(:,3)==1)+1;
    
    D.feedbackType = NeuronsFR_DB_EJ008_Batch3(s).data(:,4);
    
end