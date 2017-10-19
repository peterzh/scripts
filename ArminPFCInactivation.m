%% Load data

% sessions = {'2017-09-14_1_ALK053';
%             '2017-09-15_1_ALK053';
%             '2017-09-18_1_ALK053';
%             '2017-09-19_1_ALK053';
%             '2017-09-20_1_ALK053';
%             '2017-09-21_1_ALK053';
%             '2017-09-14_2_ALK058';
%             '2017-09-15_4_ALK058';
%             '2017-09-18_1_ALK058';
%             '2017-09-19_1_ALK058';
%             '2017-09-20_1_ALK058';
%             '2017-09-21_1_ALK058'};
        
% sessions = {'2017-09-26_4_ALK053';
%             '2017-09-27_2_ALK053';
%             '2017-09-28_1_ALK053';
%             '2017-09-26_3_ALK058';
%             '2017-09-27_1_ALK058';
%             '2017-09-28_1_ALK058'};
        
%         sessions = {'2017-09-26_4_ALK053';
%             '2017-09-27_2_ALK053';
%             '2017-09-28_1_ALK053';
%             '2017-09-14_1_ALK053';
%             '2017-09-15_1_ALK053';
%             '2017-09-18_1_ALK053';
%             '2017-09-19_1_ALK053';
%             '2017-09-20_1_ALK053';
%             '2017-09-21_1_ALK053';
%             };
        
        sessions = {
            '2017-10-02_2_ALK053';
            '2017-10-03_2_ALK053';
            '2017-10-04_1_ALK053';
            '2017-10-02_2_ALK058';
            '2017-10-03_1_ALK058';
            };
        

D = struct;
numSessions = length(sessions);
for sess = 1:numSessions
    %Load
    d = loadData(sessions{sess});
    
    d.stim1D = d.stimulus(:,2) - d.stimulus(:,1);
    
    %Remove repeatNum>1 trials
    d = structfun(@(f) f(d.repeatNum==1,:), d, 'uni', 0);
    
    %Add session ID
    d.sessionID = ones(size(d.response))*sess;
    
    %Add Mouse ID
    switch(sessions{sess}(end-5:end))
        case 'ALK053'
            d.mouseID = ones(size(d.response))*1;
        case 'ALK058'
            d.mouseID = ones(size(d.response))*2;
    end
    
    D = addstruct(D,d);
    
end

D.sessionID = categorical(D.sessionID,1:length(sessions),sessions);

g = gramm('x',D.RT,'color',D.sessionID,'subset',D.response~=3);
g.facet_grid(num2cell(num2str(D.mouseID)),[]);
g.stat_density();
g.set_names('x','Reaction Time','color','Session');

figure;
g.draw();


%% Plot per-session
testCont = [linspace(1,0,100)' zeros(100,1); zeros(100,1) linspace(0,1,100)'];

figure('color','w');
for sess = 1:numSessions
%     figure('name',sessions{sess},'color','w');
    
    d = getrow(D,double(D.sessionID)==sess & D.laserType==0);
    
    %Plot model fit on non-laser data
    g=GLM(d).setModel('C50-subset').fit;
    ph = g.calculatePhat(g.parameterFits,testCont);
    
    ii=1;
    for r = [1 3 2]
        subplot(numSessions,3,ii + 3*(sess-1) ); hold on;
        plot(testCont(:,2)-testCont(:,1),ph(:,r),'k','linewidth',0.5);
        xlim([-1 1]); ylim([0 1]); 
        ii = ii+1;
        
        if r == 1
            ylabel({sessions{sess}(end-5:end), sessions{sess}(1:end-7)});
%             set(gca,'interpreter','tex');
        end
    end
    
    %Plot real non-laser data
    uC = unique(d.stim1D);
    for c = 1:length(uC)
        r = d.response(d.stim1D==uC(c));
        [ph,pci] = binofit(sum([r==1 r==2 r==3],1),length(r));
        
        ii = 1;
        for r=[1 3 2]
            subplot(numSessions,3,ii + 3*(sess-1) ); hold on;
            l=line([1 1]*uC(c),pci(r,:));
            plot(uC(c),ph(r),'k.','markersize',20);
            set(l,'Color','k','Linewidth',0.5);
            ii = ii+1;
        end
    end


    d = getrow(D,double(D.sessionID)==sess & D.laserType==1);
    
    %Plot model fit on laser data
    g=GLM(d).setModel('C50-subset').fit;
    ph = g.calculatePhat(g.parameterFits,testCont);
    
    ii=1;
    for r = [1 3 2]
        subplot(numSessions,3,ii + 3*(sess-1) ); hold on;
        plot(testCont(:,2)-testCont(:,1),ph(:,r),'r','linewidth',0.5);
        xlim([-1 1]); ylim([0 1]); 
        ii = ii+1;
    end
    
    %Plot real non-laser data
    uC = unique(d.stim1D);
    for c = 1:length(uC)
        r = d.response(d.stim1D==uC(c));
        [ph,pci] = binofit(sum([r==1 r==2 r==3],1),length(r));
        
        ii = 1;
        for r=[1 3 2]
            subplot(numSessions,3,ii + 3*(sess-1) ); hold on;
            l=line([1 1]*uC(c),pci(r,:));
            plot(uC(c),ph(r),'r.','markersize',20);
            set(l,'Color','r','Linewidth',0.5);
            ii = ii+1;
        end
    end
    
end

axs = get(gcf,'children');
set(axs(4:end),'xtick','','xcolor','w');
% set(axs,'ytick','','ycolor','w');
idx = (1:3*numSessions); idx(3:3:end) = [];
set(axs(idx),'ytick','','ycolor','w');

%% Plot pooled-session
testCont = [linspace(1,0,100)' zeros(100,1); zeros(100,1) linspace(0,1,100)'];

figure('name','pooled session','color','w');
d = getrow(D,D.laserType==0); 
% d = getrow(D,D.laserType==0 & double(D.sessionID)<=6); %ALK053
% d = getrow(D,D.laserType==0 & double(D.sessionID)>6); %ALK058

%Plot model fit on non-laser data
g=GLM(d).setModel('C50-subset').fit;
ph = g.calculatePhat(g.parameterFits,testCont);

ii=1;
for r = [1 3 2]
    subplot(1,3,ii); hold on;
    plot(testCont(:,2)-testCont(:,1),ph(:,r),'k','linewidth',0.5);
    xlim([-1 1]); ylim([0 1]);
    ii = ii+1;
end

%Plot real non-laser data
uC = unique(d.stim1D);
for c = 1:length(uC)
    r = d.response(d.stim1D==uC(c));
    [ph,pci] = binofit(sum([r==1 r==2 r==3],1),length(r));
    
    ii = 1;
    for r=[1 3 2]
        subplot(1,3,ii); hold on;
        l=line([1 1]*uC(c),pci(r,:));
        plot(uC(c),ph(r),'k.','markersize',20);
        set(l,'Color','k','Linewidth',0.5);
        ii = ii+1;
    end
end


d = getrow(D,D.laserType==1); 
% d = getrow(D,D.laserType==1 & double(D.sessionID)<=6); %ALK053
% d = getrow(D,D.laserType==1 & double(D.sessionID)>6); %ALK058

%Plot model fit on laser data
g=GLM(d).setModel('C50-subset').fit;
ph = g.calculatePhat(g.parameterFits,testCont);

ii=1;
for r = [1 3 2]
    subplot(1,3,ii); hold on;
    plot(testCont(:,2)-testCont(:,1),ph(:,r),'r','linewidth',0.5);
    xlim([-1 1]); ylim([0 1]);
    ii = ii+1;
end

%Plot real non-laser data
uC = unique(d.stim1D);
for c = 1:length(uC)
    r = d.response(d.stim1D==uC(c));
    [ph,pci] = binofit(sum([r==1 r==2 r==3],1),length(r));
    
    ii = 1;
    for r=[1 3 2]
        subplot(1,3,ii); hold on;
        l=line([1 1]*uC(c),pci(r,:));
        plot(uC(c),ph(r),'r.','markersize',20);
        set(l,'Color','r','Linewidth',0.5);
        ii = ii+1;
    end
end

axs = get(gcf,'children');
set(axs(1:2),'ytick','','ycolor','w');

    

