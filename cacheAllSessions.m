function cacheAllSessions()

subjects = {'Nyx','Beadle','Bovet','Vin','Keynes','Heinz','Spemann','Whipple','Morgan','Murphy','Chomsky'};

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
        
        
        [D,meta] = loadData( expRefs{sess} );
        
        if ~isempty(D.response) & length(D.response)>10
            %Categorise this session
            t = table;
            t.mouse = subjects(subj);
            t.expRef = expRefs(sess);
            t.numberTrials = length(D.response);
            t.taskType = categorical(max(D.response),2:3,{'2AFC','2AUC'});
            t.stimType = categorical( any(min(D.stimulus,[],2)>0) , [0,1], {'Detection','Discrimination'});
            t.performance = mean(D.feedbackType==1);
            
            if isfield(D,'laserPower') && max(D.laserPower)>0
                t.laserPower = max(D.laserPower);
                t.laserDuration = max(D.laserDuration);
                
                laserCoords = unique(D.laserCoord(~isnan(D.laserCoord(:,1)),:),'rows');
                t.numLaserCoords = size(laserCoords,1);
                
                regions = unique(D.laserRegion(~isnan(D.laserCoord(:,1))));
                t.laserRegions = {strjoin(cellstr(regions),' ')};
            else
                t.laserPower = 0;
                t.laserDuration = 0;
                t.numLaserCoords = 0;
                t.laserRegions = {''};
            end
            
            t.notes = {meta.notes};
            if t.performance < 0.6 || t.numberTrials < 200 || contains(t.notes,'exclude') || contains(t.notes,'Exclude')
                t.exclude = 1;
            else
                t.exclude = 0;
            end
            
            
%             disp(t);
            bigTable = [bigTable; t];
        end
    end
    
    
    
    
end

writetable(bigTable,'C:\Users\Peter\Desktop\LocalSessionData\SESSION_TABLE.csv');




end