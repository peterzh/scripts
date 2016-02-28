
%% Bar plot summarising the effect of the laser on GLM in V1 and barrel cortical areas
figDir = 'B:\figures\GLM+Laser';
MODEL = 'C50-subset';

expRefs = {'2015-09-17_1_Hopkins',1,1;
           '2015-09-22_2_Hopkins',1,1;
           '2015-10-05_1_Hopkins',1,1;
           '2015-10-09_1_Hopkins',1,1;
           '2015-10-14_1_Hopkins',1,1;
           '2015-10-16_1_Hopkins',1,1;
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
            '2015-10-22_1_Hopkins',2,4};
        
figLabel = {'HopkinsRightV1','HopkinsLeftV1','HopkinsRightBarrel','HopkinsLeftBarrel'};
% 
% expRefs = {'2015-09-30_1_Eijkman',1,1;
% '2015-10-04_1_Eijkman',1,1;
% '2015-10-15_1_Eijkman',1,1;
% '2015-10-23_1_Eijkman',3,1;
% 
%     '2015-09-28_1_Eijkman',1,2;
% '2015-10-02_1_Eijkman',1,2;        
% '2015-10-07_1_Eijkman' ,1,2;           
% '2015-10-08_1_Eijkman' ,1,2;          
% '2015-10-10_1_Eijkman' ,1,2;           
% '2015-10-23_1_Eijkman' ,1,2;           
% '2015-10-23_1_Eijkman' ,2,2;           
% '2015-11-18_1_Eijkman' ,1,2; 
% '2015-11-19_1_Eijkman' ,1,2;
% 
% '2015-09-27_1_Eijkman', 1,3;
% '2015-10-01_3_Eijkman' ,2,3;
% '2015-10-22_1_Eijkman' ,2,3;
% 
% '2015-09-29_1_Eijkman',  1,4;
% '2015-09-30_1_Eijkman',  2,4;
% '2015-10-15_1_Eijkman', 2,4;
% '2015-10-16_1_Eijkman', 2,4;
% };
% figLabel = {'EijkmanRightV1','EijkmanLeftV1','EijkmanRightBarrel','EijkmanLeftBarrel'};

D=struct;
% P_Laser = nan(size(expRefs,1),6);
% P_noLaser = nan(size(expRefs,1),6);
for b=1:size(expRefs,1)
    expRef = expRefs{b,1};
    laserID = expRefs{b,2};
    area = expRefs{b,3};
    gl = laserGLM(expRef);
    ROW = gl.data;
    ROW = getrow(ROW,ROW.laserIdx==0 | ROW.laserIdx==laserID);
    ROW.area = ones(length(ROW.laser),1)*area;
    D=addstruct(D,ROW);
end
D.laserIdx(D.laserIdx>1) = 1;

g = GLM(D).setModel(MODEL);
% g.parameterBounds(:,end) = [0.2 0.2];
% g.regularise  = @(b)(10*sum(b.^2));
contrast1D = D.contrast_cond(:,2)-D.contrast_cond(:,1);
Cmax = max(abs(contrast1D));

linewidth=2;
pDiff=[];
figure;
i=1;
for a = 1:length(figLabel)
    g.data = getrow(D,D.area==a & D.laserIdx==0);
    gNoLaser = g.fit;
    
    g.data = getrow(D,D.area==a & D.laserIdx==1);
    gLaser = g.fit;
    
    zl=subplot(length(figLabel),3,i);
    s=ezplot(@(x)(gNoLaser.ZL(gNoLaser.parameterFits,[-x 0])),[0 -Cmax]);
    set(s,'linestyle','--','Color','black','linewidth',linewidth);
    hold on;
    bL = gNoLaser.parameterFits(1);
    kL = gNoLaser.parameterFits(2);
    line([0 -Cmax],[bL bL]);
    line([0 -Cmax],[bL+kL bL+kL]);
    l=get(zl,'children');
    set(l([1,2]),'lineStyle','--');
    
    s=ezplot(@(x)(gLaser.ZL(gLaser.parameterFits,[-x 0])),[0 -Cmax]);
    set(s,'linestyle','-','Color','black','linewidth',linewidth);
    xlabel('C_L');
    ylabel('Z_L');
    xlim([-Cmax 0]);
    axis square;
    set(gca,'box','off');
    title(figLabel{a});
    bL = gLaser.parameterFits(1);
    kL = gLaser.parameterFits(2);
    line([0 -Cmax],[bL bL]);
    line([0 -Cmax],[bL+kL bL+kL]);
    hold off;
    
    zr=subplot(length(figLabel),3,i+1);
    s=ezplot(@(x)(gNoLaser.ZR(gNoLaser.parameterFits,[0 x])),[0 Cmax]);
    set(s,'linestyle','--','Color','black','linewidth',linewidth);
    hold on;
    bL = gNoLaser.parameterFits(3);
    kL = gNoLaser.parameterFits(4);
    line([0 Cmax],[bL bL]);
    line([0 Cmax],[bL+kL bL+kL]);
    l=get(zr,'children');
    set(l([1,2]),'lineStyle','--');
    
    s=ezplot(@(x)(gLaser.ZR(gLaser.parameterFits,[0 x])),[0 Cmax]);
    set(s,'linestyle','-','Color','black','linewidth',linewidth);
    xlabel('C_R');
    ylabel('Z_R');
    xlim([-Cmax 0]);
    axis square;
    set(gca,'box','off');
    title('');
    bL = gLaser.parameterFits(3);
    kL = gLaser.parameterFits(4);
    line([0 Cmax],[bL bL]);
    line([0 Cmax],[bL+kL bL+kL]);
    hold off;
    
    xlim([0 Cmax]);
    linkaxes([zl zr],'y');
    set([zr zl],'box','off')
    set(gcf,'color','w');
    ylim([-5 5]);
    set([zl zr],'fontsize',10);
    
    b(a) = subplot(length(figLabel),3,i+2);
    pChangeZL = abs(gLaser.parameterFits(1) - gNoLaser.parameterFits(1)) + abs(gLaser.parameterFits(2) - gNoLaser.parameterFits(2));
    pChangeZR = abs(gLaser.parameterFits(3) - gNoLaser.parameterFits(3)) + abs(gLaser.parameterFits(4) - gNoLaser.parameterFits(4));
%     bar(gLaser.parameterFits(1:4)-gNoLaser.parameterFits(1:4));
    bar([pChangeZL pChangeZR]);
    set(gca,'xticklabels',{'change ZL','change ZR'});
    ylabel('|Offset* - Offset| + |Scale* - Scale|')
    i=i+3;
    
%     bar(gLaser.parameterFits - gNoLaser.parameterFits);
%     title(figLabel{a});
%     pDiff(a,:) = gLaser.parameterFits - gNoLaser.parameterFits;
end
linkaxes(b,'y');
set(b,'box','off');