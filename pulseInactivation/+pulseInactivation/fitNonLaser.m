pulseInactivation.config;

%% Load each preproc dataset, fit GLM on non-laser data

preprocFiles = dir( fullfile(preprocDir,'*.mat') );

ALL_E = [];
for file = 1:length(preprocFiles)
    
    load( fullfile(preprocFiles(file).folder, preprocFiles(file).name) );
    D = D(D.laserType==0,:); %Only get nonlaser data
    D = D(D.repeatNum==1,:); %Only get repeatNum==1 data
    
    E = struct;
    E.stimulus = D.stimulus;
    E.response = nan(size(D.response));
    E.response(D.response=='Left')=1;
    E.response(D.response=='Right')=2;
    E.response(D.response=='NoGo')=3;
    
    g = GLM(E).setModel('C50-subset').fit;
    
    fig = g.plotFit;
    set(fig,'Position',[0 0 600 1000]);
    print(fig,[figuresDir '\' preprocFiles(file).name(1:end-4) '_GLM'],'-dpdf','-bestfit');
    close all;
    
    ALL_E = [ALL_E; E];
end

g = GLM(ALL_E).setModel('C50-subset').fit;
fig = g.plotFit;
set(fig,'Position',[0 0 600 1000]);
print(fig,[figuresDir '\' 'ALL_GLM'],'-dpdf','-bestfit');