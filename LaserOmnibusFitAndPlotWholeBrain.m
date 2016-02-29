%% Load all sessions into one, marking each session separately.

subject = {'Murphy','Spemann','Morgan','Whipple'};
% subject = {'Whipple'};
expRefs = dat.listExps(subject)';
expRefs = vertcat(expRefs{:});
% 
% expRefs = {
%     '2016-02-15_1_Whipple';
% '2016-02-16_1_Whipple';
% '2016-02-17_1_Whipple';
% '2016-02-19_1_Whipple';
% % '2016-02-22_1_Whipple';
% '2016-02-23_2_Whipple';
% '2016-02-24_1_Whipple';
% '2016-02-25_1_Whipple';
% '2016-02-15_2_Morgan';
% '2016-02-18_1_Morgan';
% '2016-02-23_1_Morgan';
% '2016-02-24_1_Morgan';
% '2016-02-16_1_Spemann';
% '2016-02-17_1_Spemann';
% '2016-02-18_1_Spemann';
% '2016-02-19_1_Spemann';
% % '2016-02-22_2_Spemann';
% '2016-02-23_2_Spemann';
% '2016-02-24_1_Spemann';
% '2016-02-25_1_Spemann';
% '2016-02-24_1_Murphy';
% '2016-02-25_1_Murphy';
% '2016-02-25_2_Murphy'
% 
% };
%        

%Combine all sessions
D=struct;
i=1;
figure;
for s = 1:length(expRefs)
    disp([num2str(s) '/' num2str(length(expRefs))]);
    try
        d = laserGLM(expRefs{s}).data;
        d = structfun(@(x)(x(6:end,:)),d,'uni',0); %trim first 5 trials
        if length(d.response) < 100
            error('not enough trials');
        end
        
        e = getrow(d,d.laserIdx==0);
        g=GLM(e).setModel('C50-subset').fit;
%         g.plotFit; drawnow;
        c50sub.n(i) = g.parameterFits(5);
        c50sub.c50(i) = g.parameterFits(6);
        subplot(7,7,i); g.plotFit;drawnow;
        title(gca,'');
        xlabel(gca,'');
        ylabel(gca,'');
        
        g=GLM(e).setModel('C^N-subset').fit;
        nsub.n(i) = g.parameterFits(5);
        
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

figure; subplot(1,2,1); hist(c50sub.n); xlabel('n');
subplot(1,2,2); hist(c50sub.c50); xlabel('c50');
D = rmfield(D,'laserIdx');
% D = getrow(D,D.repeatNum==1);
l = laserGLM(D);
% save('C:\Users\Peter\Desktop\scripts\LaserOmnibus_data.mat', 'l');

%% Define and fit model 
numSessions = max(l.data.sessionID);
numSites = max(l.data.laserIdx);
numTrials = size(l.data.response,1);

Y = l.data.response;

%construct design matrix X (for loop for ease of understanding not speed!)
% XsessionOffsets = zeros(numTrials,numSessions);
% XsessionContrasts = zeros(numTrials,numSessions);
% XsiteOffsets = zeros(numTrials,numSites);
% XsiteContrasts = zeros(numTrials,numSites);

% X = zeros(numTrials,3*numSessions + 3*numSites);
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

%% Plot model parameters over sessions and sites
% Organise parameter estimate values
sessionP = b(1:3*numSessions,:);
sessionP_Biases = sessionP(1:3:end,:);
sessionP_CLeft = sessionP(2:3:end,:);
sessionP_CRight = sessionP(3:3:end,:);

%Plot the session by session values of the parameters to see whether they
%change much
figure;
subplot(3,1,1); bar(sessionP_Biases); title('Biases for each session');
subplot(3,1,2); bar(sessionP_CLeft); title('CL scaling for each session');
subplot(3,1,3); bar(sessionP_CRight); title('CR scaling for each session');
set(gca,'XTick',1:length(expRefs),'XTickLabel',expRefs,'XTickLabelRotation',90);

sitesP = b(3*numSessions+1:end,:);
sitesP_Biases = sitesP(1:3:end,:);
sitesP_CLeft = sitesP(2:3:end,:);
sitesP_CRight = sitesP(3:3:end,:);

%Plot the site by site values of the parameters to see whether they
%change much
figure;
subplot(3,1,1); bar(sitesP_Biases); title('Biases for each site');
subplot(3,1,2); bar(sitesP_CLeft); title('CL scaling for each site');
subplot(3,1,3); bar(sitesP_CRight); title('CR scaling for each site');

%Plot the map of these site by site parameter values

%% nice plot: overlay onto kirkcaldie brain
sitesP = b(3*numSessions+1:end,:);
sitesP_Biases = sitesP(1:3:end,:);
sitesP_CLeft = sitesP(2:3:end,:);
sitesP_CRight = sitesP(3:3:end,:);

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
