%% Load all sessions into one, marking each session separately.

subject = {'Murphy','Spemann','Morgan','Whipple'};
% subject = {'Whipple'};
expRefs = dat.listExps(subject)';
expRefs = vertcat(expRefs{:}); 

expRefs = {'2016-03-23_1_Murphy';
           '2016-03-24_1_Murphy';
           '2016-03-25_4_Murphy';
           '2016-03-24_1_Spemann'};
       
%Combine all sessions
D=struct;
i=1;
for s = 1:length(expRefs)
    disp([num2str(s) '/' num2str(length(expRefs))]);
    try
        d = laserGLM(expRefs{s}).data;
        d = structfun(@(x)(x(6:(end-14),:)),d,'uni',0); %trim first 5 trials and last 15
        if length(d.response) < 100
            error('not enough trials');
        end
        
        if max(d.laserIdx) < 50
            disp(expRefs{s});
            error('wrong data')
        end
        
        e = getrow(d,d.laserIdx==0);
        g=GLM(e).setModel('C50-subset').fit;
        c50sub.n(i) = g.parameterFits(5);
        c50sub.c50(i) = g.parameterFits(6);

%         g=GLM(e).setModel('C^N-subset').fit;
%         nsub.n(i) = g.parameterFits(5);
        
%         d.sessionID = ones(length(d.response),1)*s;
        d.sessionID = ones(length(d.response),1)*i;
        i=i+1;
        D = addstruct(D,d);
    catch
        warning(expRefs{s});
        expRefs{s}=[];
    end
end
expRefs(cellfun(@isempty,expRefs))=[];

% figure; subplot(1,2,1); hist(c50sub.n); xlabel('n');
% subplot(1,2,2); hist(c50sub.c50); xlabel('c50');
D = rmfield(D,'laserIdx');
D = getrow(D,D.repeatNum==1);
l = laserGLM(D);
% save('C:\Users\Peter\Desktop\scripts\LaserOmnibus_data.mat', 'l');

%% Define and fit model 
numSessions = max(l.data.sessionID);
numSites = max(l.data.laserIdx);
numTrials = size(l.data.response,1);

Y = l.data.response;

%construct design matrix X (for loop for ease of understanding not speed!)
X = sparse(numTrials,3*numSessions + 3*numSites);

% N = 0.4;
% cfn = @(c)(c.^N);
% c50=1;
% N=1;
% cfn = @(c)((c.^N)/(c.^N + c50^N));
contrast_rep_fcn = @(c,n,c50)((c.^n)/(c.^n + c50^n));
for i = 1:numTrials

    thisSession = l.data.sessionID(i);
    cfn = @(c)contrast_rep_fcn(c,c50sub.n(thisSession),c50sub.c50(thisSession));
    X(i,3*thisSession - 2) = 1;
    X(i,3*thisSession - 1) = cfn(l.data.contrast_cond(i,1));
    X(i,3*thisSession - 0) = cfn(l.data.contrast_cond(i,2));
    
    thisSite = l.data.laserIdx(i);
    if thisSite ~= 0
        X(i,3*numSessions + 3*thisSite - 2) = 1;
        X(i,3*numSessions + 3*thisSite - 1) = cfn(l.data.contrast_cond(i,1));
        X(i,3*numSessions + 3*thisSite - 0) = cfn(l.data.contrast_cond(i,2));
        
    end    
end

opts=struct;
opts.intr=0; %don't add a global intercept
opts.alpha=1; %lasso regularisation
% opts.alpha=0;

% penalty = ones(1,size(X,2)); penalty(1:3:end)=0;
% opts.penalty_factor=penalty; %Try penalising only the contrast terms
fit=cvglmnet(X,Y,'multinomial',glmnetSet(opts));
cvglmnetPlot(fit);
b=cvglmnetCoef(fit);
b=[b{1}-b{3} b{2}-b{3}];
b(1,:) = [];

% Organise parameter estimate values
sessionP = b(1:3*numSessions,:);
sessionP_Biases = sessionP(1:3:end,:);
sessionP_CLeft = sessionP(2:3:end,:);
sessionP_CRight = sessionP(3:3:end,:);

sitesP = b(3*numSessions+1:end,:);
sitesP_Biases = sitesP(1:3:end,:);
sitesP_CLeft = sitesP(2:3:end,:);
sitesP_CRight = sitesP(3:3:end,:);

%% Plot model parameters over sessions and sites

%Plot the session by session values of the parameters to see whether they
%change much
figure;
subplot(3,1,1); bar(sessionP_Biases); title('Biases for each session');
subplot(3,1,2); bar(sessionP_CLeft); title('CL scaling for each session');
subplot(3,1,3); bar(sessionP_CRight); title('CR scaling for each session');
set(gca,'XTick',1:length(expRefs),'XTickLabel',expRefs,'XTickLabelRotation',90);

%Plot the site by site values of the parameters to see whether they
%change much
figure;
subplot(3,1,1); bar(sitesP_Biases); title('Biases for each site');
subplot(3,1,2); bar(sitesP_CLeft); title('CL scaling for each site');
subplot(3,1,3); bar(sitesP_CRight); title('CR scaling for each site');

%% nice plot: overlay onto kirkcaldie brain
toDisplay = {sitesP_Biases,sitesP_CLeft,sitesP_CRight};%, sitesP_CLeft, sitesP_CRight};
labels = {'bias','left contrast sensitivity','right contrast sensitivity'};
% toDisplay = {sitesP_CRight};
% labels = {'right contrast sensitivity'}
dotSize=150;
dotShape='o';
kimg=imread('B:\stuff\kirkcaldie_brain_BW.PNG');

a=1;
figure;
side={'log(\pi_L/\pi_{NG})','log(\pi_R/\pi_{NG})'};
for d = 1:length(toDisplay)
    for lr=1:2
        subplot(length(toDisplay),2,a);
        imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal');
        set(imX,'alphadata',0.7);
        hold on;
        s(lr)=scatter(l.inactivationSite(:,2),l.inactivationSite(:,1),dotSize,toDisplay{d}(:,lr),dotShape,'filled'); axis equal; colorbar;
        ylim([-6 4]);
        caxis([-1 1]*max(abs(toDisplay{d}(:))));
%         set(s(lr),'markeredgecolor',[1 1 1]*1,'linewidth',0);
        hold off;
        title(['Additive effect of laser on ' side{lr} ' ' labels{d}])
        xlim([-5 5]);
        a=a+1;
    end
end
% 
cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
    linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
colormap(cmap);

%% % Plot number of repetitions in each site
figure;
subplot(1,2,1);
imagesc(X);
hold on;
pos=0;
for s = 1:length(expRefs)
    eRef = expRefs{s};
    numTrials = sum(D.sessionID==s);

    tx=text(50,pos + 0.5*numTrials,[eRef ' n=' num2str(numTrials)],'interpreter','none');
    tx.Color=[1 1 1];
    pos = pos + numTrials;
end
hold off;

subplot(1,2,2);
tab = tabulate(l.data.laserIdx);
tab = tab(2:end,2);
scatter(l.inactivationSite(:,2),l.inactivationSite(:,1),dotSize,tab,dotShape,'filled'); axis equal; colorbar; colormap(gca,'parula')
axis([-5 5 -4 3]);
% set(s,'markeredgecolor',[0.8 0.8 0.8]);
title('Number of trials at each site');

%% Split-half reliability
numTrials = length(Y);
C = cvpartition(numTrials,'KFold',3);

for split = 1:length(C.TrainSize)
    Ys = Y(C.training(split));
    Xs = X(C.training(split),:);
    fit=cvglmnet(X,Y,'multinomial',glmnetSet(opts));
    b=cvglmnetCoef(fit);
    b=[b{1}-b{3} b{2}-b{3}];
    b(1,:) = [];
    
    sitesP = b(3*numSessions+1:end,:);
    sitesP_Biases_split(:,:,split) = sitesP(1:3:end,:);
    sitesP_CLeft_split(:,:,split) = sitesP(2:3:end,:);
    sitesP_CRight_split(:,:,split) = sitesP(3:3:end,:);
end

figure;
h(1)=subplot(3,1,1);
Diff = std(sitesP_Biases_split,[],3);
hist(Diff); title('std in biases for all sites');
h(2)=subplot(3,1,2);
Diff = std(sitesP_CLeft_split,[],3);
hist(Diff); title('std in left contrast sens for all sites');
h(3)=subplot(3,1,3);
Diff = std(sitesP_CRight_split,[],3);
hist(Diff); title('std in right contrast sens for all sites');
% linkaxes(h,'x');

%% Plot all psychometric curves (select noLaser session and plot all laser additive curves)
% Only possible for 1D contrast trials since the contrast can be projected
% onto 1 axis.

sessionID = 1;
g = g.setModel('C50');
testC = [linspace(1,0,100)' zeros(100,1);
         zeros(100,1) linspace(0,1,100)'];

paramsNL = [sessionP_Biases(sessionID,:);
          sessionP_CLeft(sessionID,:);
          sessionP_CRight(sessionID,:)];
% paramsNL = [mean(sessionP_Biases);
%             mean(sessionP_CLeft);
%             mean(sessionP_CRight)];
% 
% paramsNL = [0 0;
%             3 0;
%             0 3];
        
shape = [c50sub.n(sessionID) c50sub.c50(sessionID)];
% shape = [mean(c50sub.n) mean(c50sub.c50)];
      
figure;
kimg=imread('B:\stuff\kirkcaldie_brain_BW.PNG');
imX=image(-5:1:5,4:-1:-6,kimg); set(gca,'ydir','normal');
set(gca,'Position',[0 0 1 1]);

for site = 1:size(l.inactivationSite,1)
    paramsL = [sitesP_Biases(site,:);
               sitesP_CLeft(site,:);
               sitesP_CRight(site,:)];
    paramsL = paramsL + paramsNL;
    posIn = l.inactivationSite(site,:)/10 + [0.58 0.465];
    axes;

    hold on;
    phat_NL = g.calculatePhat([paramsNL(:);shape'],testC);
    phat_NL = phat_NL(:,1:2);
    plot(diff(testC,[],2),phat_NL,':','linewidth',1.5);
    
    set(gca, 'ColorOrderIndex', 1);
    
    phat_L = g.calculatePhat([paramsL(:);shape'],testC);
    phat_L = phat_L(:,1:2);
    
%     plot(diff(testC,[],2),phat_L,'-');
    
    if sum(sum(abs(phat_NL-phat_L))) < 0.001 %Prevents weird display bug when phat_NL = phat_L exactly
        phat_L = phat_L + 0.001;
    end
    
    fill([diff(testC,[],2); flipud(diff(testC,[],2))],[phat_NL(:,1);flipud(phat_L(:,1))],[0 0.4470 0.7410],'linestyle','none','facealpha',0.5);
    fill([diff(testC,[],2); flipud(diff(testC,[],2))],[phat_NL(:,2);flipud(phat_L(:,2))],[0.8500 0.3250 0.0980],'linestyle','none','facealpha',0.5);
    hold off;  
    
    set(gca,'Position',[posIn(2) posIn(1) 0.07 0.07],'box','off','XColor',[0 0 0 0],'YColor',[0 0 0 0]);
end
% set(gcf,'color','white');

%% (NO GLM) simple map of change in % choice at C=0 over the laser sites, median over sessions

for s = 1:length(expRefs)
    noLaser_responses = l.data.response(l.data.laserIdx==0 & diff(l.data.contrast_cond,[],2)==0 & l.data.sessionID==s);
    noLaser_p = arrayfun(@(r)(sum(noLaser_responses==r)/length(noLaser_responses)),1:3);
    
    Laser_p=[];
    for site = 1:size(l.inactivationSite,1)
        Laser_responses = l.data.response(l.data.laserIdx==site & diff(l.data.contrast_cond,[],2)==0 & l.data.sessionID==s);
        Laser_p(site,:) = arrayfun(@(r)(sum(Laser_responses==r)/length(Laser_responses)),1:3);
    end
    
    pDiff(:,:,s) = bsxfun(@minus,Laser_p,noLaser_p);
end
pDiff = nanmedian(pDiff,3);

figure;
titles = {'pL_{laser} - pL_{noLaser}','pR_{laser} - pR_{noLaser}','pNG_{laser} - pNG_{noLaser}'};
for lr=1:3
    subplot(1,3,lr);
    scatter(l.inactivationSite(:,2),l.inactivationSite(:,1),300,pDiff(:,lr),'o','filled'); axis equal; colorbar;
    ylim([-6 4]);
    xlim([-4 4]);
    caxis([-1 1]*max(abs(pDiff(:))));
    title(titles{lr});
end

cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
    linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
colormap(cmap);