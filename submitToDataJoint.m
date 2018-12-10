function submitToDataJoint(subj)

%get mice names
subjects = fetch(subj);

baseDirs = {'\\zserver.cortexlab.net\Data2\Subjects';
    '\\zserver.cortexlab.net\Data\expInfo';
    '\\zserver.cortexlab.net\Data\Subjects';
    '\\zubjects.cortexlab.net\Subjects'};

%For each subject
for subj = 1:length(subjects)
   disp(subjects(subj));
    
    %Get the expRefs for these existing sessions
    existing_expRefs = fetchn(d.Session & subjects(subj),'concat(session_date,"_",session_num,"_",mouse_name)->expRef');
    
    
    %Now get all expRefs on the server
    expRefs = {};
    for i = 1:4
        subjDir = fullfile(baseDirs{i},subjects(subj).mouse_name);
        blockFiles = dirPlus(subjDir,'filefilter','\_Block.mat$','struct',true);
        if ~isempty(blockFiles)
            eRefs = cellfun(@(bf) strrep(bf, '_Block.mat', ''), {blockFiles.name}, 'uni' ,0)';
            expRefs = [expRefs; eRefs];
        end
    end
    expRefs = unique(expRefs);
    expRefs = setdiff(expRefs,existing_expRefs);
    
    %Get only the expRefs which aren't on the server
    
    %Go through each session, submit to server
    for sess = 1:length(expRefs)
        [subject, date, seq] = dat.parseExpRef(expRefs{sess});
        try
            [D,meta] = loadData(expRefs{sess},'force');
            close all;
            try
                D.laserRegion = cellstr(D.laserRegion);
                D.prev_laserRegion = cellstr(D.prev_laserRegion);
            catch
            end
            
            key = struct('mouse_name',subject,...
                'session_date',datestr(date,29),...
                'session_num',seq,...
                'num_trials', length(D.response),...
                'performance', mean(D.feedbackType==1),...
                'stimulus_type', char( categorical( any(min(D.stimulus,[],2)>0) , [0,1], {'Detection','Discrimination'}) ),...
                'choice_type', char( categorical(max(D.response),2:3,{'2AFC','2AUC'}) ) );
            if key.num_trials > 0
                key.rig = meta.rig;
                key.data = D;
                key.project_id = identifyProject(key);
                insert(d.Session, key);
            end
            
        catch me
            disp(me);
        end
    end
    
end

parpopulate(d.Trial,subjects );
parpopulate(d.SessionPerformance,subjects);
parpopulate(d.SessionChronometric,subjects);
parpopulate(d.SessionPsychometric, proj(d.Session & 'choice_type="2AFC"'), proj(d.GLM & 'model_type="Binomial"'), subjects );
parpopulate(d.SessionPsychometric, proj(d.Session & 'choice_type="2AUC"'), proj(d.GLM & 'model_type="Multinomial"'), subjects );
end


function proj_id = identifyProject(key)

if any(contains(key.rig,{'zredone','zredtwo','zredthree','zgreyfour','zym3'}))
    proj_id = 'training';
    return;
end

if isfield(key.data,'laserType') & any(key.data.laserType>0)
    
    %                 %Number of coordinates
    %                 numCoords = size( unique(key.data.laserCoord(key.data.laserType~=0,:),'rows'), 1);
    
    if any(contains(key.rig,{'zym1','zym2'})) %Blue rigs
        
        %pulse or long-duration?
        if max(key.data.laserDuration) == 1.5
            prefix = 'galvo_';
        elseif max(key.data.laserDuration) == 0.025
            prefix = 'galvoPulse_';
        else
            keyboard;
        end
        
    elseif strcmp(key.rig,'zgood') %Nick's rig
        prefix = 'sparse_';
    end
    
    %Unilateral or bilateral?
    if max(key.data.laserType)==1
        proj_id = [prefix 'unilateral'];
    elseif max(key.data.laserType)==2
        proj_id = [prefix 'bilateral'];
    else
        proj_id = [key.choice_type '_' key.stimulus_type];
    end
else
    proj_id = [key.choice_type '_' key.stimulus_type];
end

end