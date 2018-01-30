function bigTable = supermegaultrasessionsearcher()
sessionListFile = '\\basket.cortexlab.net\homes\peterzh\sessionList.csv';

try
    bigTable = readtable(sessionListFile);
    error('hi');
catch
    subjects = {'Nyx','Vin','Keynes','Spemann','Whipple','Morgan','Murphy','Chomsky','Heinz'};
    
    baseDirs = {'\\zserver.cortexlab.net\Data2\Subjects';
        '\\zserver.cortexlab.net\Data\expInfo';
        '\\zserver.cortexlab.net\Data\Subjects'};
    
    bigTable = table;
    for subj = 1:length(subjects)
        
        expRefs = {};
        for i = 1:3
            subjDir = fullfile(baseDirs{i},subjects{subj});
            blockFiles = dirPlus(subjDir,'filefilter','\_Block.mat$','struct',true);
            if ~isempty(blockFiles)
                eRefs = cellfun(@(bf) strrep(bf, '_Block.mat', ''), {blockFiles.name}, 'uni' ,0)';
                expRefs = [expRefs; eRefs];
            end
        end
        expRefs = unique(expRefs);
        
        
        for sess = 1:length(expRefs)
            %         try
            t = table;
            t.mouse = categorical(subj,1:length(subjects),subjects);
            t.expRef = expRefs(sess);
            
            
            try
                D = loadData(expRefs{sess});
            catch
                D.response = [];
            end
            
            if ~isempty(D.response)
                %                 t.data = D;
                
                t.numberTrials = length(D.response);
                t.taskType = categorical(max(D.response),2:3,{'2AFC','2AUC'});
                t.stimType = categorical( any(min(D.stimulus,[],2)>0) , [0,1], {'Detection','Discrimination'});
                t.performance = mean(D.feedbackType==1);
                t.laserPower = max(D.laserPower);
                
                isGoodPerf = t.performance>0.6;
                isManyTrials = t.numberTrials > 200;
                
                if isGoodPerf && isManyTrials
                    t.exclude = 0;
                else
                    t.exclude = 1;
                end
                
                bigTable = [bigTable; t];
            end
            %         catch
            %         end
        end
        
        
        
        
    end
    
    writetable(bigTable, sessionListFile, 'FileType', 'spreadsheet');
end

    
    

end