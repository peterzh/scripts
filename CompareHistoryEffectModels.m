
%% %%%%%%%%%%%%%%% ONE SESSION HISTORY EFFECT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expRef = '2015-09-21_1_Hopkins';
g = GLM(expRef);
models = {'Offset','C^N-subset','C^N-subset-hist-response','C^N-subset-hist-contrast','C^N-subset-hist-success'};

for m = 1:length(models)
    g1 = g.setModel(models{m});
    
    switch(models{m})
        case 'C^N-subset-hist-response'
            g1.data.hist = [circshift(g1.data.response,1)==1 circshift(g1.data.response,1)==2 circshift(g1.data.response,1)==3];
            g1.data.hist(1,:) = zeros(1,3);
        case 'C^N-subset-hist-contrast'
            C1d = g1.data.contrast_cond(:,2) - g1.data.contrast_cond(:,1);
            g1.data.hist = [circshift(C1d,1)<0 circshift(C1d,1)==0 circshift(C1d,1)>0];
            g1.data.hist(1,:) = zeros(1,3);
        case 'C^N-subset-hist-success'
            C1d = g1.data.contrast_cond(:,2) - g1.data.contrast_cond(:,1);
            instructedR = C1d;
            instructedR(C1d<0) = 1;
            instructedR(C1d>0) = 2;
            instructedR(C1d==0) = 3;

            g1.data.hist = circshift(g.data.response==instructedR,1);
            g1.data.hist(1)=0;
    end
%     keyboard;
    g1 = g1.fitCV;
    save(fullfile(saveDir,[expRef '_model-'  models{m} '.mat']),'g1');
   
end

%% Load and plot
LL=[];
for m = 1:length(models)
    load(fullfile(saveDir,[expRef '_model-'  models{m} '.mat']));
    LL(m)=-sum(log2(g1.p_hat));
end

bar(LL);
set(gca,'XTickLabel',models);
ylabel('-SUM log_2 likelihood (small=better)');


%% %%%%%%%%%%%%%%% ALL SESSION HISTORY EFFECT %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subject = 'Hopkins';
expRefs = dat.listExps(subject);
saveDir = '\\basket.cortexlab.net\homes\peterzh\NeuroGLM-experimental\ModelFiles';

alldat = struct;
for f = 1:length(expRefs)
    try
        try
            g = laserGLM(expRefs{f});
        catch
            g = GLM(expRefs{f});
            g.data.laserIdx = zeros(length(g.data.response),1);
            g.data.laser = nan(length(g.data.response),2);
        end
        
        if length(g.data.laserIdx) ~= length(g.data.response)
            error('Mismatch laser + behav data length');
        end
        
        alldat = addstruct(alldat,g.data);
        disp(expRefs{f});

    catch
%         warning('error');
    end
end

D = getrow(alldat,alldat.laserIdx==0);
g = GLM(D);

%% Try many models
models = {'C^N-subset','C^N-subset-hist-response','C^N-subset-hist-contrast','C^N-subset-hist-success'};
LLs = [];
for m = 1:length(models)
    
    g1 = g.setModel(models{m});
    switch(models{m})
        case 'C^N-subset-hist-response'
            g1.data.hist = [circshift(g1.data.response,1)==1 circshift(g1.data.response,1)==2 circshift(g1.data.response,1)==3];
            g1.data.hist(1,:) = zeros(1,3);
        case 'C^N-subset-hist-contrast'
            C1d = g1.data.contrast_cond(:,2) - g1.data.contrast_cond(:,1);
            g1.data.hist = [circshift(C1d,1)<0 circshift(C1d,1)==0 circshift(C1d,1)>0];
            g1.data.hist(1,:) = zeros(1,3);
        case 'C^N-subset-hist-success'
            C1d = g1.data.contrast_cond(:,2) - g1.data.contrast_cond(:,1);
            instructedR = C1d;
            instructedR(C1d<0) = 1;
            instructedR(C1d>0) = 2;
            instructedR(C1d==0) = 3;

            g1.data.hist = circshift(g.data.response==instructedR,1);
            g1.data.hist(1)=0;
    end
    
    g1 = g1.fitCV(50);
    save(fullfile(saveDir,[subject '_AllSessions_model-'  g1.modelString '.mat']),'g1');
    disp('done');
end

%% Load and plot models

LL=[];
for m = 1:length(models)
    load(fullfile(saveDir,[subject '_AllSessions_model-'  models{m} '.mat']));
    LL(:,m) = -log2(g1.p_hat);
end

boxplot(LL,'sym','r','labels',{'C^N-subset','Response_{t-1}','Contrast_{t-1}','Success_{t-1}'});
set(findobj(gca,'tag','Whisker'),'Visible','off');
set(gca,'fontsize',17);
ylabel('Model likelihood [bits per trial]');
hold on; plot(mean(LL),'.','markersize',10);hold off;

