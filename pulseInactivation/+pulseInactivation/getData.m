pulseInactivation.config;

%% Go through each subject, identify every data folder on the server associated with that subject
expRefs = cell(1,length(subjects));
for subj = 1:length(subjects)
    
    %Find all Timeline.mat files
    tl_files = dirPlus( fullfile(rawDataDir, subjects{subj}),...
        'FileFilter', '\_Timeline.mat', 'PrependPath', false);
    
    %Get the expRef associated with each
    parts = cellfun( @(file)  strsplit(file,'_') , tl_files, 'uni', 0);
    expRefs{subj} = cellfun( @(file) strjoin(file(1:3),'_') , parts, 'uni', 0);
end

%% Load each dataset, preprocess, check if meet requirements for the experiment
sessionList = table;
sessionList.expRef = cat(1,expRefs{:});
sessionList.good = zeros(length(allExpRefs),1);
for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    
    fprintf('%s\n',eRef);
    
    try
        D = loadData(eRef);
        
        is2D = any(min(D.stimulus,[],2) > 0);
        isManyTrials = length(D.response)>150;
        isRightDate = datenum(eRef(1:10),'yyyy-mm-dd') >= datenum('2017-07-11','yyyy-mm-dd');
        isUnilateral = max(D.laserType) == 1;
        isPulsed = max(D.laserDuration) == 0.025;
        isRandomLaser = length(unique(D.laserOnset))>20;
        
        if is2D && isManyTrials && isRightDate && isUnilateral && isPulsed && isRandomLaser
            
            if makePlots
                E = getrow(D,D.laserType==0);
                g = GLM(E).setModel('C50-subset').fit;
                fig = g.plotFit;
                
                set(fig,'Position',[0 0 600 1000]);
                print(fig,fullfile(figuresDir,'GLM_nonLaserFits',eRef),'-dpdf','-bestfit');
                close all;
            end
            
            sessionList.good(sess) = 1;
        end
    catch
    end
end

%Save sessionList
save(fullfile(preprocDir,'sessionList.mat'),'sessionList');

%% Extract data and combine into a single data format
load(fullfile(preprocDir,'sessionList.mat'));
filteredSessionList = sessionList(sessionList.good==1,:);

ALLD = struct;
for sess = 1:height(filteredSessionList)
    D = loadData(filteredSessionList.expRef{sess});
    D.sessionID = ones(size(D.response))*sess;
    
    ALLD = addstruct(ALLD,D);
end

E = getrow(ALLD,ALLD.laserType==0);
g = GLM(E).setModel('C50-subset').fit;
fig = g.plotFit;
                
                set(fig,'Position',[0 0 600 1000]);
                print(fig,fullfile(figuresDir,'GLM_nonLaserFits',eRef),'-dpdf','-bestfit');
                close all;
