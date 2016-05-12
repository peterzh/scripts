%Shows long term progress of mouse performance data
name = 'Hopkins';
expRefs = dat.listExps(name);

%extract timestamp and performance on that session
%first time
formatIn = 'dd-mmm-yyyy';
i=1; time=[]; performance=[];
for session = 1:length(expRefs)
    try
        block = dat.loadBlock(expRefs{session});
        
        time(i) = datenum(block.startDateTimeStr(1:11),formatIn);
        performance(i) = mean([block.trial.feedbackType]==1);
        i=i+1;
    catch
    end
end
%plot
figure; plot(time-time(1),performance*100,'o--'); title(name); xlabel('days'); ylabel('% correct');