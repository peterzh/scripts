%% %%%%%%%%%%%%%  FOR PAPER: create publishable plot showing how the model can capture laser effects
% Do it by splitting the data into laser and nolaser and fitting

figDir = 'B:\figures\GLM+Laser';
MODEL = 'C50-subset';


% %Eijkman Left V1
% expRefs = {'2015-09-28_1_Eijkman',1;
%            '2015-10-02_1_Eijkman',1;
%            '2015-10-07_1_Eijkman',1;
%            '2015-10-08_1_Eijkman',1;
%            '2015-10-10_1_Eijkman',1;
%            '2015-10-23_1_Eijkman',1;
%            '2015-10-23_1_Eijkman',2;
%            '2015-11-18_1_Eijkman',1; %many contrast levels
%            '2015-11-19_1_Eijkman',1}; %many contrast levels
% figLabel = 'LeftV1';

% 
% %Eijkman Right V1
% expRefs = {'2015-09-30_1_Eijkman',1;
%            '2015-10-04_1_Eijkman',1;
%            '2015-10-15_1_Eijkman',1;
%            '2015-10-23_1_Eijkman',3};
% figLabel = 'RightV1';
% % % 
% % % 
% % %Hopkins Left V1
% expRefs = {'2015-09-16_1_Hopkins',1;
%            '2015-10-03_1_Hopkins',1;
%            '2015-10-06_1_Hopkins',1;
%            '2015-10-07_1_Hopkins',1;
% %            '2015-10-08_1_Hopkins',1;
%            '2015-10-08_1_Hopkins',2;
%            '2015-10-14_1_Hopkins',2
%            '2015-10-22_1_Hopkins',1
%            '2015-10-30_1_Hopkins',1};
% % figLabel = 'LeftV1';
% % % 
% % %Hopkins Right V1
expRefs = {'2015-09-17_1_Hopkins',1;
           '2015-09-22_2_Hopkins',1;
           '2015-10-05_1_Hopkins',1;
           '2015-10-09_1_Hopkins',1;
           '2015-10-14_1_Hopkins',1;
%            '2015-10-16_1_Hopkins',1;
           '2015-10-19_4_Hopkins',1;
           '2015-10-20_2_Hopkins',2;
           '2015-10-20_2_Hopkins',3};
figLabel = 'RightV1';

% expRefs = {'2015-09-16_1_Hopkins',1;
%            '2015-10-03_1_Hopkins',1;
%            '2015-10-06_1_Hopkins',1;
%            '2015-10-07_1_Hopkins',1;
%            '2015-10-08_1_Hopkins',1;
%            '2015-10-08_1_Hopkins',2;
%            '2015-10-14_1_Hopkins',2;
%            '2015-10-22_1_Hopkins',1};
% %            '2015-10-30_1_Hopkins',1};
% figLabel = 'LeftV1';


    grp_cols =      [ 0    0.4470    0.7410;
        0.8500    0.3250    0.0980;
        0.4940    0.1840    0.5560;
        0.4660    0.6740    0.1880;
        0.3010    0.7450    0.9330;
        0.6350    0.0780    0.1840];

% %Hopkins left barrel
% expRefs = {'2015-10-01_1_Hopkins' ,2;
% '2015-10-02_1_Hopkins' ,2;
% '2015-10-13_1_Hopkins' ,3;
% '2015-10-22_1_Hopkins' ,2};
% figLabel = 'LeftBarrel';

P_Laser = nan(size(expRefs,1),6);
P_noLaser = nan(size(expRefs,1),6);
for b=1:size(expRefs,1)

    expRef = expRefs{b,1};
    laserID = expRefs{b,2};
    
    gl = laserGLM(expRef);
    if size(gl.inactivationSite,1) > 1
        warning(['Block contains multiple inactivation sites. Using ' num2str( gl.inactivationSite(laserID,:) )]);
    end
    
    g = GLM(expRef);
%     g.regularise = @(b)(sum(b.^2)*0.1);
    
    % Make plots
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Overlap C50 models V2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('name',['expRef: ' g.expRef ' Inactivation Site: [' num2str(gl.inactivationSite(laserID,:)) '] Model: ' MODEL]);
    g = g.setModel(MODEL);
    
    % NO LASER CALCULATIONS
    g.data = getrow(gl.data,gl.data.laserIdx==0);
    g = g.fit;
    
    P_noLaser(b,:) = g.parameterFits;
    
    %extract datapoints with errors
    contrast1D = g.data.contrast_cond(:,2) - g.data.contrast_cond(:,1);
    uniqueC1D_noLaser = unique(contrast1D);
    prop_noLaser=[];
    prop_ci_noLaser=nan(3,2,length(uniqueC1D_noLaser));
    binom_se_noLaser=[];
    for c = 1:length(uniqueC1D_noLaser)
        D = getrow(g.data,contrast1D == uniqueC1D_noLaser(c));
%         p = sum([D.response==1 D.response==2 D.response==3])/length(D.response);
        [p,pci]=binofit(sum([D.response==1 D.response==2 D.response==3]),length(D.response),0.05);
        prop_noLaser = [prop_noLaser;p];
        
        %two different ways of doing the error
        prop_ci_noLaser(:,:,c) = pci;
        bse = sqrt((p.*(1-p)/length(D.response)));
        binom_se_noLaser = [binom_se_noLaser;bse];
    end
    
    %calculate prediction phat
    maxC = max(max(g.data.contrast_cond));
    evalC = [linspace(maxC,0,100)', zeros(100,1);
        zeros(100,1), linspace(0,maxC,100)'];
    evalC1d = evalC(:,2) - evalC(:,1);
    
    otherInputs = g.Zinput(g.data);
    otherInputs(:,1:2)=[];
    
    if isempty(otherInputs)
        inputs = evalC;
    else
        inputs = [evalC, zeros(length(evalC),size(otherInputs,2))];
    end
    
    phat_noLaser = g.calculatePhat(g.parameterFits,inputs);
    
    % LASER CALCULATIONS
    g.data = getrow(gl.data,gl.data.laserIdx==laserID);
    g = g.fit;
    
    P_Laser(b,:) = g.parameterFits;
    
    %extract datapoints with errors
    contrast1D = g.data.contrast_cond(:,2) - g.data.contrast_cond(:,1);
    uniqueC1D_Laser = unique(contrast1D);
    prop_Laser=[];
    binom_se_Laser=[];
    prop_ci_Laser=nan(3,2,length(uniqueC1D_noLaser));
    for c = 1:length(uniqueC1D_Laser)
        D = getrow(g.data,contrast1D == uniqueC1D_Laser(c));
        %         p = sum([D.response==1 D.response==2 D.response==3])/length(D.response);
        [p,pci]=binofit(sum([D.response==1 D.response==2 D.response==3]),length(D.response),0.05);
        prop_Laser = [prop_Laser;p];
        bse = sqrt((p.*(1-p)/length(D.response)));
        binom_se_Laser = [binom_se_Laser;bse];
        prop_ci_Laser(:,:,c) = pci;

    end
    
    %Calculate prediction phats
    phat_Laser = g.calculatePhat(g.parameterFits,inputs);
    
    labels = {'p left','p right','p nogo'};
    for r = [1,2,3]
        subplot(3,2,2*r-1);

        hold on;
        %plot data with 95% confidence intervals 
        err_noLaser = [prop_noLaser(:,r)-squeeze(prop_ci_noLaser(r,1,:)), squeeze(prop_ci_noLaser(r,2,:))-prop_noLaser(:,r)];
        err_Laser = [prop_Laser(:,r)-squeeze(prop_ci_Laser(r,1,:)), squeeze(prop_ci_Laser(r,2,:))-prop_Laser(:,r)];
        h=errorbar(uniqueC1D_noLaser,prop_noLaser(:,r),err_noLaser(:,1),err_noLaser(:,2),'o','MarkerSize',5,'Color',grp_cols(r,:));
        
        
        errorbar(uniqueC1D_Laser,prop_Laser(:,r),err_Laser(:,1),err_Laser(:,2),'.','MarkerSize',20,'Color',grp_cols(r,:));

%         %plot data with binomial standard error (smaller)
%         err_noLaser = binom_se_noLaser;
%         err_Laser = binom_se_Laser;
%         errorbar(uniqueC1D_noLaser,prop_noLaser(:,r),err_noLaser(:,r),'o','MarkerSize',5,'Color',grp_cols(r,:));
%         errorbar(uniqueC1D_Laser,prop_Laser(:,r),err_noLaser(:,r),'.','MarkerSize',20,'Color',grp_cols(r,:));
        
        %plot fits
        plot(evalC1d,phat_noLaser(:,r),'--','Color',grp_cols(r,:));
        plot(evalC1d,phat_Laser(:,r),'Color',grp_cols(r,:));
        
        ylim([0,1]);
        ylabel(labels{r});
        xlabel('Contrast');
        hold off;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Overlap C50 models V3 plotting Z functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subplot(2,2,2);
    Cmax = max(abs(contrast1D));
    g = g.setModel(MODEL);
    
    %non-laser Z fcn
    g.data = getrow(gl.data,gl.data.laserIdx==0);
    g = g.fit;
    s=ezplot(@(x)(g.ZR(g.parameterFits,[0 x])),[0 Cmax]);
%     set(s,'linestyle','--','Color',grp_cols(4,:));
    set(s,'linestyle','--','Color','black');

    hold on;
    s=ezplot(@(x)(g.ZL(g.parameterFits,[-x 0])),[0 -Cmax]);
    set(s,'linestyle','--','Color','black');
    idxMid = find(uniqueC1D_noLaser==0);
    CR = uniqueC1D_noLaser(idxMid:end);
    CL = uniqueC1D_noLaser(1:idxMid);
    scatter(CR,g.ZR(g.parameterFits,[zeros(length(CR),1) CR]),'black','s');
    scatter(CL,g.ZL(g.parameterFits,[-CL zeros(length(CL),1)]),'black','s');

    
    %laser Z fcn
    g.data = getrow(gl.data,gl.data.laserIdx==laserID);
    g = g.fit;
    s=ezplot(@(x)(g.ZR(g.parameterFits,[0 x])),[0 Cmax]);
    set(s,'linestyle','-','Color','black');
    hold on; 
    s=ezplot(@(x)(g.ZL(g.parameterFits,[-x 0])),[0 -Cmax]);
    set(s,'linestyle','-','Color','black');
    idxMid = find(uniqueC1D_Laser==0);
    CR = uniqueC1D_Laser(idxMid:end);
    CL = uniqueC1D_Laser(1:idxMid);
    scatter(CR,g.ZR(g.parameterFits,[zeros(length(CR),1) CR]),'black','s','filled');
    scatter(CL,g.ZL(g.parameterFits,[-CL zeros(length(CL),1)]),'black','s','filled');

    title('');
    xlabel('Contrast');
    ylabel('Log odds');
    hold off;
    
    axis auto
    xlim([-Cmax Cmax]);
    axis square;
    set(gca,'box','off');
    
%     keyboard;
    set(gcf,'color','w');
    print(fullfile(figDir,[expRef '_' figLabel '_' num2str(laserID) '.pdf' ]),'-dpdf','-painters');
end

% numP = size(BIAS_noLaser,1);

% subplot(1,2,1);
% plot([BIAS_noLaser, SCALE_noLaser],[BIAS_Laser, SCALE_Laser],'.','MarkerSize',20);
% hold on; 
% legend('Bias ZL','Bias ZR','Sens ZL','Sens ZR');
% ezplot(gca,'y=x'); 
% hold off;
% xlabel('Param noLaser');
% ylabel('Param Laser');
% axis(gca,'equal','auto');
% title(figLabel);

% plot(P_noLaser,P_Laser,'.','MarkerSize',20);
% hold on;
% legend(g.parameterLabels);
% ezplot(gca,'y=x'); 
% hold off;
% xlabel('Param noLaser');
% ylabel('Param Laser');
% axis(gca,'equal','auto');
% title(figLabel);


%% %%%%%%%%%%%%%  FOR PAPER: create summary plot
% Do it by splitting the data into laser and nolaser and fitting
    grp_cols =      [ 0    0.4470    0.7410;
        0.8500    0.3250    0.0980;
        0.4940    0.1840    0.5560;
        0.4660    0.6740    0.1880;
        0.3010    0.7450    0.9330;
        0.6350    0.0780    0.1840];
    
figDir = 'B:\figures\GLM+Laser';
MODEL = 'C50-subset';

expRefs = {'2015-10-03_1_Hopkins', 3;
'2015-10-13_1_Hopkins', 2;
'2015-10-22_1_Hopkins', 3;}
figLabel = 'RightBarrel';

D=struct;
for b=1:size(expRefs,1)
    expRef = expRefs{b,1};
    laserID = expRefs{b,2};
    gl = laserGLM(expRef);
    ROW = gl.data;
    ROW = getrow(ROW,ROW.laserIdx==0 | ROW.laserIdx==laserID);
    D=addstruct(D,ROW);
end
D.laserIdx(D.laserIdx>1) = 1;

g = GLM(D);
% g.regularise = @(b)(sum(b.^2)*0);
g = g.setModel(MODEL);

% NO LASER CALCULATIONS
g.data = getrow(g.data,g.data.laserIdx==0);
gNoLaser = g.fit;

%extract datapoints with errors
contrast1D = gNoLaser.data.contrast_cond(:,2) - gNoLaser.data.contrast_cond(:,1);
uniqueC1D_noLaser = unique(contrast1D);
prop_noLaser=[];
prop_ci_noLaser=nan(3,2,length(uniqueC1D_noLaser));
binom_se_noLaser=[];
for c = 1:length(uniqueC1D_noLaser)
    E = getrow(gNoLaser.data,contrast1D == uniqueC1D_noLaser(c));
    %         p = sum([D.response==1 D.response==2 D.response==3])/length(D.response);
    [p,pci]=binofit(sum([E.response==1 E.response==2 E.response==3]),length(E.response),0.05);
    prop_noLaser = [prop_noLaser;p];
    
    %two different ways of doing the error
    prop_ci_noLaser(:,:,c) = pci;
    bse = sqrt((p.*(1-p)/length(E.response)));
    binom_se_noLaser = [binom_se_noLaser;bse];
end

%calculate prediction phat
maxC = max(max(gNoLaser.data.contrast_cond));
evalC = [linspace(maxC,0,100)', zeros(100,1);
    zeros(100,1), linspace(0,maxC,100)'];
evalC1d = evalC(:,2) - evalC(:,1);

otherInputs = g.Zinput(gNoLaser.data);
otherInputs(:,1:2)=[];

if isempty(otherInputs)
    inputs = evalC;
else
    inputs = [evalC, zeros(length(evalC),size(otherInputs,2))];
end

phat_noLaser = gNoLaser.calculatePhat(gNoLaser.parameterFits,inputs);

% LASER CALCULATIONS
g.data = getrow(D,D.laserIdx==1);
gLaser = g.fit;
    
%extract datapoints with errors
contrast1D = gLaser.data.contrast_cond(:,2) - gLaser.data.contrast_cond(:,1);
uniqueC1D_Laser = unique(contrast1D);
prop_Laser=[];
binom_se_Laser=[];
prop_ci_Laser=nan(3,2,length(uniqueC1D_noLaser));
for c = 1:length(uniqueC1D_Laser)
    E = getrow(gLaser.data,contrast1D == uniqueC1D_Laser(c));
    %         p = sum([D.response==1 D.response==2 D.response==3])/length(D.response);
    [p,pci]=binofit(sum([E.response==1 E.response==2 E.response==3]),length(E.response),0.05);
    prop_Laser = [prop_Laser;p];
    bse = sqrt((p.*(1-p)/length(E.response)));
    binom_se_Laser = [binom_se_Laser;bse];
    prop_ci_Laser(:,:,c) = pci;
    
end

%Calculate prediction phats
phat_Laser = gLaser.calculatePhat(gLaser.parameterFits,inputs);

labels = {'p left','p right','p nogo'};
figure;
for r = [1,2,3]
    subplot(3,2,2*r-1);
    
    hold on;
    %plot data with 95% confidence intervals
    err_noLaser = [prop_noLaser(:,r)-squeeze(prop_ci_noLaser(r,1,:)), squeeze(prop_ci_noLaser(r,2,:))-prop_noLaser(:,r)];
    err_Laser = [prop_Laser(:,r)-squeeze(prop_ci_Laser(r,1,:)), squeeze(prop_ci_Laser(r,2,:))-prop_Laser(:,r)];
    h=errorbar(uniqueC1D_noLaser,prop_noLaser(:,r),err_noLaser(:,1),err_noLaser(:,2),'o','MarkerSize',5,'Color',grp_cols(r,:));
    
    
    errorbar(uniqueC1D_Laser,prop_Laser(:,r),err_Laser(:,1),err_Laser(:,2),'.','MarkerSize',20,'Color',grp_cols(r,:));
    
    %         %plot data with binomial standard error (smaller)
    %         err_noLaser = binom_se_noLaser;
    %         err_Laser = binom_se_Laser;
    %         errorbar(uniqueC1D_noLaser,prop_noLaser(:,r),err_noLaser(:,r),'o','MarkerSize',5,'Color',grp_cols(r,:));
    %         errorbar(uniqueC1D_Laser,prop_Laser(:,r),err_noLaser(:,r),'.','MarkerSize',20,'Color',grp_cols(r,:));
    
    %plot fits
    plot(evalC1d,phat_noLaser(:,r),'--','Color',grp_cols(r,:));
    plot(evalC1d,phat_Laser(:,r),'Color',grp_cols(r,:));
    
    ylim([0,1]);
    ylabel(labels{r});
    xlabel('Contrast');
    hold off;
end

zl=subplot(2,4,3); %ZL
Cmax = max(abs(contrast1D));
s=ezplot(@(x)(gNoLaser.ZL(gNoLaser.parameterFits,[-x 0])),[0 -Cmax]); 
set(s,'linestyle','--','Color','black');
hold on;
bL = gNoLaser.parameterFits(1);
kL = gNoLaser.parameterFits(2);
line([0 -Cmax],[bL bL]);
line([0 -Cmax],[bL+kL bL+kL]);
l=get(zl,'children');
set(l([1,2]),'lineStyle','--');

s=ezplot(@(x)(gLaser.ZL(gLaser.parameterFits,[-x 0])),[0 -Cmax]);
set(s,'linestyle','-','Color','black');
xlabel('C_L');
ylabel('Z_L');
xlim([-Cmax 0]);
axis square;
set(gca,'box','off');
title('');
bL = gLaser.parameterFits(1);
kL = gLaser.parameterFits(2);
line([0 -Cmax],[bL bL]);
line([0 -Cmax],[bL+kL bL+kL]);
hold off;

zr=subplot(2,4,4); %ZR
Cmax = max(abs(contrast1D));
s=ezplot(@(x)(gNoLaser.ZR(gNoLaser.parameterFits,[0 x])),[0 Cmax]); 
set(s,'linestyle','--','Color','black');
hold on;
bL = gNoLaser.parameterFits(3);
kL = gNoLaser.parameterFits(4);
line([0 Cmax],[bL bL]);
line([0 Cmax],[bL+kL bL+kL]);
l=get(zr,'children');
set(l([1,2]),'lineStyle','--');

s=ezplot(@(x)(gLaser.ZR(gLaser.parameterFits,[0 x])),[0 Cmax]);
set(s,'linestyle','-','Color','black');
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
%     keyboard;
print(fullfile(figDir,[figLabel '.pdf' ]),'-dpdf','-painters');
savefig(fullfile(figDir,[figLabel '.fig' ]));
