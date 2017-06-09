%% For each session, fit laser and nolaser data separately. Plot param changes
% get data
figLabels = {'RightVisual','LeftVisual','RightBarrel','LeftBarrel'};
subjects = {'Hopkins','Eijkman'};%,'Morgan','Whipple','Murphy','Spemann'};
expRefs = {};
for region = 1:length(figLabels)
    refs = LaserIdentifySessions(subjects,figLabels{region});
    refs(:,3)={region};
    expRefs=[expRefs;refs];
end

%% 
P_Laser = nan(size(expRefs,1),6);
P_noLaser = nan(size(expRefs,1),6);

ZL20_laser = nan(size(expRefs,1),1);
ZR20_laser = nan(size(expRefs,1),1);

ZL20_nolaser = nan(size(expRefs,1),1);
ZR20_nolaser = nan(size(expRefs,1),1);

pDiff = nan(size(expRefs,1),2);
NUMTRIALS=[];
for b=1:size(expRefs,1)
    disp([num2str(b) '/' num2str(size(expRefs,1))]);
    expRef = expRefs{b,1};
    laserID = expRefs{b,2};
    
    gl = laserGLM(expRef);
    
    g = GLM(expRef);
    
    % Make plots
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Overlap C50 models V2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     figure('name',['expRef: ' g.expRef ' Inactivation Site: [' num2str(gl.inactivationSite(laserID,:)) '] Model: ' MODEL]);
    g = g.setModel('C50-subset');
    g.regularise = @(b)(0.05*sum(b.^2));
%     g.regularise = @(b)(0.001*sum(b.^2));
    
    % NO LASER CALCULATIONS
    NUMTRIALS(b) = sum(gl.data.laserIdx==0 | gl.data.laserIdx==laserID);
    g.data = getrow(gl.data,gl.data.laserIdx==0);
    g = g.fit;
    P_noLaser(b,:) = g.parameterFits;
%     
%     %Identify contrast at 50% of one choice
    testC = linspace(0,1,1000);
    ph=g.calculatePhat(g.parameterFits,[linspace(0,1,1000)' zeros(1000,1)]);
    CL = testC(find(ph(:,1)>0.5,1,'first'));
    ph=g.calculatePhat(g.parameterFits,[zeros(1000,1) linspace(0,1,1000)']);
    CR = testC(find(ph(:,2)>0.5,1,'first'));
    
    if isempty(CL)
        CL = 0.08;
    end
    
    if isempty(CR)
        CR = 0.08;
    end
    
%     ZL20_nolaser(b) = g.ZL(P_noLaser(b,:),[0.2 0]);
%     ZR20_nolaser(b) = g.ZR(P_noLaser(b,:),[0 0.2]);
    ZL20_nolaser(b) = g.ZL(P_noLaser(b,:),[CL 0]);
    ZR20_nolaser(b) = g.ZR(P_noLaser(b,:),[0 CR]);
    %LASER CALCULATION
    g.data = getrow(gl.data,gl.data.laserIdx==laserID);
    g = g.fit;
    P_Laser(b,:) = g.parameterFits;
    
%     ZL20_laser(b) = g.ZL(P_Laser(b,:),[0.2 0]);
%     ZR20_laser(b) = g.ZR(P_Laser(b,:),[0 0.2]);
    ZL20_laser(b) = g.ZL(P_Laser(b,:),[CL 0]);
    ZR20_laser(b) = g.ZR(P_Laser(b,:),[0 CR]);
    
%     pDL = g.calculatePhat(P_noLaser(b,:),[0.2 0]) - g.calculatePhat(P_Laser(b,:),[0.2 0]);
%     pDR = g.calculatePhat(P_noLaser(b,:),[0 0.2]) - g.calculatePhat(P_Laser(b,:),[0 0.2]);
    pDL = g.calculatePhat(P_noLaser(b,:),[CL 0]) - g.calculatePhat(P_Laser(b,:),[CL 0]);
    pDR = g.calculatePhat(P_noLaser(b,:),[0 CR]) - g.calculatePhat(P_Laser(b,:),[0 CR]);

    pDiff(b,1) = pDL(1);
    pDiff(b,2) = pDR(2);
end

% s(sum((s<RANGE(1) | s>RANGE(2)),2)>0,:)=nan
mouseID = cell2mat(cellfun(@(c)(strcmp(c(14),'H')),expRefs(:,1),'uni',0)); %1 Hopkins, 0 Eijkman

figure;
areaID = cell2mat(expRefs(:,3));
% pdiff = abs(P_Laser - P_noLaser);
s = [ZR20_nolaser - ZR20_laser, ZL20_nolaser - ZL20_laser];

%normality test

%t-Tests
% [H,P,CI,STATS] = ttest2( s(areaID==2,1) , s(areaID==4,1) , 'tail','both');
% disp(['delta ZR: ' figLabels{2} ' [M=' num2str(mean(s(areaID==2,1))) ' SD=' num2str(std(s(areaID==2,1)))   '] vs ' figLabels{4}  ' [M=' num2str(mean(s(areaID==4,1))) ' SD=' num2str(std(s(areaID==4,1)))   '] ' 'dof=' num2str(STATS.df) ' t=' num2str(STATS.tstat) ' p=' num2str(P)])
[P,H,STATS]=ranksum( s(areaID==2,1) , s(areaID==4,1) );
disp(['delta ZR: ' figLabels{2} ' [M=' num2str(mean(s(areaID==2,1))) ' SD=' num2str(std(s(areaID==2,1)))   '] vs ' figLabels{4}  ' [M=' num2str(mean(s(areaID==4,1))) ' SD=' num2str(std(s(areaID==4,1)))   '] ' ' z=' num2str(STATS.zval) ' p=' num2str(P)])



% 
% [H,P,CI,STATS] = ttest2( s(areaID==1,2) , s(areaID==3,2) , 'tail','both');
% disp(['delta ZL: ' figLabels{1} ' [M=' num2str(mean(s(areaID==1,2))) ' SD=' num2str(std(s(areaID==1,2)))   '] vs ' figLabels{3}  ' [M=' num2str(mean(s(areaID==3,2))) ' SD=' num2str(std(s(areaID==3,2)))   '] ' 'dof=' num2str(STATS.df) ' t=' num2str(STATS.tstat) ' p=' num2str(P)])

[P,H,STATS]=ranksum( s(areaID==1,2) , s(areaID==3,2) );
disp(['delta ZL: ' figLabels{1} ' [M=' num2str(mean(s(areaID==1,2))) ' SD=' num2str(std(s(areaID==1,2)))   '] vs ' figLabels{3}  ' [M=' num2str(mean(s(areaID==3,2))) ' SD=' num2str(std(s(areaID==3,2)))   '] ' ' z=' num2str(STATS.zval) ' p=' num2str(P)])




% s = [ZL20_laser./ZL20_nolaser, ZR20_laser./ZR20_nolaser];

% s = pDiff;
% s(sum((s<RANGE(1) | s>RANGE(2)),2)>0,:)=nan

% s = [sum(pdiff(:,1:2),2) sum(pdiff(:,3:4),2)];

plot(s(areaID==1,1),s(areaID==1,2),'rs',...
    s(areaID==2,1),s(areaID==2,2),'bs',...
    s(areaID==3,1),s(areaID==3,2),'ks',...
    s(areaID==4,1),s(areaID==4,2),'ks');
ax=get(gca,'children');
set(ax(1),'markerfacecolor',[1 1 1]*1,'markeredgecolor',[1 1 1]*0.6,'linewidth',2);
set(ax(2),'markerfacecolor',[1 1 1]*1,'markeredgecolor',[1 1 1]*0.2,'linewidth',2);
set(ax(3),'markerfacecolor',[146/255 0 0],'markeredgecolor',[146/255 0 0]);
set(ax(4),'markerfacecolor',[0 146/255 146/255],'markeredgecolor',[0 146/255 146/255]);
set(ax,'markersize',8);
legend(figLabels)
axis equal;

hold on
ezplot('y=x');
title('');
xlabel('ZR_{nolaser} - ZR_{laser}'); ylabel('ZL_{nolaser} - ZL_{laser}');

% xlabel('ZL20_{laser} / ZL20_{nolaser}'); ylabel('ZR20_{laser} / ZR20_{nolaser}');

% xlabel('%L20_{nolaser} - %L20_{laser}'); ylabel('%R20_{nolaser} - %R20_{laser}');
% ylim([0 plotLim]);

m = nanmedian(s(areaID==1,:));
% m = nanmean(s(areaID==1,:));
plot(m(1),m(2),'gd','markersize',20,'markeredgecolor',[0 146/255 146/255]);

m = nanmedian(s(areaID==2,:));
% m = nanmean(s(areaID==2,:));
plot(m(1),m(2),'rd','markersize',20,'markeredgecolor',[146/255 0 0]);

m = nanmedian(s(areaID==3,:));
% m = nanmean(s(areaID==3,:));
plot(m(1),m(2),'ks','markersize',20,'markeredgecolor',[1 1 1]*0.2);

m = nanmedian(s(areaID==4,:));
% m = nanmean(s(areaID==4,:));
plot(m(1),m(2),'ko','markersize',20,'markeredgecolor',[1 1 1]*0.6);
%Add error bars
hold off;