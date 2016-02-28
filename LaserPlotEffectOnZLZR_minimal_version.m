% For each session, fit laser and nolaser data separately. Plot param changes
expRefs = {'2015-09-17_1_Hopkins',1,1;
    '2015-09-22_2_Hopkins',1,1;
    '2015-10-05_1_Hopkins',1,1;
    '2015-10-09_1_Hopkins',1,1;
    '2015-10-14_1_Hopkins',1,1;
    %            '2015-10-16_1_Hopkins',1,1;
    '2015-10-19_4_Hopkins',1,1;
    '2015-10-20_2_Hopkins',2,1;
    '2015-10-20_2_Hopkins',3,1;
    '2015-09-16_1_Hopkins',1,2;
    '2015-10-03_1_Hopkins',1,2;
    '2015-10-06_1_Hopkins',1,2;
    '2015-10-07_1_Hopkins',1,2;
    '2015-10-08_1_Hopkins',1,2;
    '2015-10-08_1_Hopkins',2,2;
    '2015-10-14_1_Hopkins',2,2;
    '2015-10-22_1_Hopkins',1,2;
    '2015-10-03_1_Hopkins',3,3;
    '2015-10-13_1_Hopkins',2,3;
    '2015-10-22_1_Hopkins',3,3;
    '2015-10-01_1_Hopkins',2,4;
    '2015-10-02_1_Hopkins',2,4;
    '2015-10-13_1_Hopkins',3,4;
    '2015-10-22_1_Hopkins',2,4;
    '2015-09-30_1_Eijkman',1,1;
    '2015-10-04_1_Eijkman',1,1;
    '2015-10-15_1_Eijkman',1,1;
    '2015-10-23_1_Eijkman',3,1;
    '2015-09-28_1_Eijkman',1,2;
    '2015-10-02_1_Eijkman',1,2;
    '2015-10-07_1_Eijkman' ,1,2;
    '2015-10-08_1_Eijkman' ,1,2;
    '2015-10-10_1_Eijkman' ,1,2;
    '2015-10-23_1_Eijkman' ,1,2;
    '2015-10-23_1_Eijkman' ,2,2;
    '2015-11-18_1_Eijkman' ,1,2;
    '2015-11-19_1_Eijkman' ,1,2;
    '2015-09-27_1_Eijkman', 1,3;
    '2015-10-01_3_Eijkman' ,2,3;
    '2015-10-22_1_Eijkman' ,2,3;
    '2015-09-29_1_Eijkman',  1,4;
    '2015-09-30_1_Eijkman',  2,4;
    '2015-10-15_1_Eijkman', 2,4;
    '2015-10-16_1_Eijkman', 2,4;
    };

figLabel = {'RightV1','LeftV1','RightBarrel','LeftBarrel'};

P_Laser = nan(size(expRefs,1),6);
P_noLaser = nan(size(expRefs,1),6);

ZL20_laser = nan(size(expRefs,1),1);
ZR20_laser = nan(size(expRefs,1),1);

ZL20_nolaser = nan(size(expRefs,1),1);
ZR20_nolaser = nan(size(expRefs,1),1);
for b=1:size(expRefs,1)
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
    
    % NO LASER CALCULATIONS
    g.data = getrow(gl.data,gl.data.laserIdx==0);
    g = g.fit;
    P_noLaser(b,:) = g.parameterFits;
    
    ZL20_nolaser(b) = g.ZL(P_noLaser(b,:),[0.2 0]);
    ZR20_nolaser(b) = g.ZR(P_noLaser(b,:),[0 0.2]);
    %LASER CALCULATION
    g.data = getrow(gl.data,gl.data.laserIdx==laserID);
    g = g.fit;
    P_Laser(b,:) = g.parameterFits;
    
    ZL20_laser(b) = g.ZL(P_Laser(b,:),[0.2 0]);
    ZR20_laser(b) = g.ZR(P_Laser(b,:),[0 0.2]);
end

figure;
areaID = cell2mat(expRefs(:,3));
% pdiff = abs(P_Laser - P_noLaser);
s = [ZL20_nolaser - ZL20_laser, ZR20_nolaser - ZR20_laser];
% s = [sum(pdiff(:,1:2),2) sum(pdiff(:,3:4),2)];
plot(s(areaID==1,1),s(areaID==1,2),'ro',...
    s(areaID==2,1),s(areaID==2,2),'bo',...
    s(areaID==3,1),s(areaID==3,2),'ko',...
    s(areaID==4,1),s(areaID==4,2),'ko');
ax=get(gca,'children');
set(ax(1),'markerfacecolor',[1 1 1]*1,'markeredgecolor',[1 1 1]*0.6,'linewidth',2);
set(ax(2),'markerfacecolor',[1 1 1]*1,'markeredgecolor',[1 1 1]*0.2,'linewidth',2);
set(ax(3),'markerfacecolor',[146/255 0 0],'markeredgecolor',[146/255 0 0]);
set(ax(4),'markerfacecolor',[0 146/255 146/255],'markeredgecolor',[0 146/255 146/255]);
set(ax,'markersize',8);
legend(figLabel)
axis equal;

hold on
ezplot('y=x',[-3 4]);
title('');
xlabel('ZL20_{nolaser} - ZL20_{laser}'); ylabel('ZR20_{nolaser} - ZR20_{laser}');
% ylim([0 plotLim]);

m = median(s(areaID==1,:));
plot(m(1),m(2),'gd','markersize',20);

m = median(s(areaID==2,:));
plot(m(1),m(2),'rd','markersize',20);

m = median(s(areaID==3,:));
plot(m(1),m(2),'ks','markersize',20);

m = median(s(areaID==4,:));
plot(m(1),m(2),'ko','markersize',20);
%Add error bars
hold off;