%% Load all sessions into one, marking each session separately.

subject = {'Murphy',116;
            'Whipple',199;
            'Morgan',198};

clear l;

for s = 1:size(subject,1)
    name = subject{s,1};
    firstSession = subject{s,2};
    
    expRefs = dat.listExps(name)';
    %     expRefs = vertcat(expRefs{firstSession:end});
    expRefs = expRefs(firstSession:end);
    
    %Combine all sessions
    D=struct;
    i=1;
    c50sub=[]; nsub=[];
    for session = 1:length(expRefs)
        disp([num2str(session) '/' num2str(length(expRefs))]);
        try
            d = laserGLM(expRefs{session},'lick').data;
            d = structfun(@(x)(x(6:(end-14),:)),d,'uni',0); %trim first 5 trials and last 15
            if length(d.response) < 100
                error('not enough trials');
            end
            
%             if max(d.laserIdx) < 50
%                 disp(expRefs{session});
%                 error('wrong data')
%             end
            
            e = getrow(d,d.laserIdx==0);
            g=GLM(e).setModel('C50-subset').fit;
            figure;g.plotFit;
            c50sub.n(i) = g.parameterFits(5);
            c50sub.c50(i) = g.parameterFits(6);
            
                    g=GLM(e).setModel('C^N-subset').fit;
                    nsub.n(i) = g.parameterFits(5);
            
            %         d.sessionID = ones(length(d.response),1)*s;
            d.sessionID = ones(length(d.response),1)*i;
            i=i+1;
            D = addstruct(D,d);
        catch
            warning(expRefs{session});
            expRefs{session}=[];
        end
    end
    expRefs(cellfun(@isempty,expRefs))=[];
    
    % figure; subplot(1,2,1); hist(c50sub.n); xlabel('n');
    % subplot(1,2,2); hist(c50sub.c50); xlabel('c50');
    D = rmfield(D,'laserIdx');
    D = getrow(D,D.repeatNum<5);
    l(s) = laserGLM(D);
    % save('C:\Users\Peter\Desktop\scripts\LaserOmnibus_data.mat', 'l');
end

%% Define and fit model 
sites = struct('Biases',{},'CLeft',{},'CRight',{},'laserSite',{});
sess = struct('Biases',{},'CLeft',{},'CRight',{});
for s = 1:size(subject,1)
    numSessions = max(l(s).data.sessionID);
    numSites = max(l(s).data.laserIdx);
    numTrials = size(l(s).data.response,1);
    
    Y = l(s).data.response;
    
    %construct design matrix X (for loop for ease of understanding not speed!)
    X = sparse(numTrials,3*numSessions + 3*numSites);
    
    % N = 0.4;
    % cfn = @(c)(c.^N);
    % c50=1;
    % N=1;
    % cfn = @(c)((c.^N)/(c.^N + c50^N));
    contrast_rep_fcn = @(c,n,c50)((c.^n)/(c.^n + c50^n));
    for i = 1:numTrials
        
        thisSession = l(s).data.sessionID(i);
        cfn = @(c)contrast_rep_fcn(c,c50sub.n(thisSession),c50sub.c50(thisSession));
        X(i,3*thisSession - 2) = 1;
        X(i,3*thisSession - 1) = cfn(l(s).data.contrast_cond(i,1));
        X(i,3*thisSession - 0) = cfn(l(s).data.contrast_cond(i,2));
        
        thisSite = l(s).data.laserIdx(i);
        if thisSite ~= 0
            X(i,3*numSessions + 3*thisSite - 2) = 1;
            X(i,3*numSessions + 3*thisSite - 1) = cfn(l(s).data.contrast_cond(i,1));
            X(i,3*numSessions + 3*thisSite - 0) = cfn(l(s).data.contrast_cond(i,2));
            
        end
    end
    
    opts=struct;
    opts.intr=0; %don't add a global intercept
    % opts.alpha=1; %lasso regularisation (L1)
    % opts.alpha=0; %ridge (L2)
    opts.alpha=0.5; %elasticnet (both)
    
    % penalty = ones(1,size(X,2)); penalty(1:3:end)=0;
    % opts.penalty_factor=penalty; %Try penalising only the contrast terms
    fit=cvglmnet(X,Y,'multinomial',glmnetSet(opts));
    cvglmnetPlot(fit);
    b=cvglmnetCoef(fit,0);
    b=[b{1}-b{3} b{2}-b{3}];
    b(1,:) = [];
    
    % Organise parameter estimate values
    sessionP = b(1:3*numSessions,:);
    sessionP_Biases = sessionP(1:3:end,:);
    sessionP_CLeft = sessionP(2:3:end,:);
    sessionP_CRight = sessionP(3:3:end,:);
    sess(s).Biases = sessionP(1:3:end,:);
    sess(s).CLeft = sessionP(2:3:end,:);
    sess(s).CRight = sessionP(3:3:end,:);
    
    sitesP = b(3*numSessions+1:end,:);
    sitesP_Biases = sitesP(1:3:end,:);
    sitesP_CLeft = sitesP(2:3:end,:);
    sitesP_CRight = sitesP(3:3:end,:);
    sites(s).Biases = sitesP(1:3:end,:);
    sites(s).CLeft = sitesP(2:3:end,:);
    sites(s).CRight = sitesP(3:3:end,:);
    sites(s).laserSite = l(s).inactivationSite;
    
    % plot design matrix and repetitions per site
    figure('name',subject{s,1});
    subplot(1,2,1);
    imagesc(X);
%     hold on;
%     pos=0;
%     for session = 1:numSessions
%         eRef = expRefs{session};
%         numTrials = sum(L.data.sessionID==session);
%         
%         tx=text(50,pos + 0.5*numTrials,[eRef ' n=' num2str(numTrials)],'interpreter','none');
%         tx.Color=[1 1 1];
%         pos = pos + numTrials;
%     end
%     hold off;
    
    subplot(1,2,2);
    tab = tabulate(l(s).data.laserIdx);
    tab = tab(2:end,2);
    scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),200,tab,'o','filled'); axis equal; colorbar; colormap(gca,'parula')
    axis([-5 5 -4 3]);
    % set(s,'markeredgecolor',[0.8 0.8 0.8]);
    title('Number of trials at each site');
end

%% Plot model parameters over sessions and sites
for s = 1:size(subject,1)
    name = subject{s,1};
    %Plot the session by session values of the parameters to see whether they
    %change much
    figure('name',name);
    subplot(3,1,1); bar(sess(s).Biases); title('Biases for each session');
    subplot(3,1,2); bar(sess(s).CLeft); title('CL scaling for each session');
    subplot(3,1,3); bar(sess(s).CRight); title('CR scaling for each session');
    set(gca,'XTick',1:length(expRefs),'XTickLabel',expRefs,'XTickLabelRotation',90);
    
    %Plot the site by site values of the parameters to see whether they
    %change much
    figure('name',name);
    subplot(3,1,1); bar(sites(s).Biases); title('Biases for each site');
    subplot(3,1,2); bar(sites(s).CLeft); title('CL scaling for each site');
    subplot(3,1,3); bar(sites(s).CRight); title('CR scaling for each site');
end

%% nice plot: overlay onto kirkcaldie brain
for s = 1:size(subject,1)
    toDisplay = {sites(s).Biases,sites(s).CLeft,sites(s).CRight};%, sitesP_CLeft, sitesP_CRight};
    labels = {'bias','CL sensitivity','CR sensitivity'};

    dotSize=150;
    dotShape='o';
    kimg=imread('\\basket.cortexlab.net\home\stuff\kirkcaldie_brain_BW.PNG');
    
    a=1;
    figure('name',subject{s,1});
    side={'log(\pi_L/\pi_{NG})','log(\pi_R/\pi_{NG})'};
    for d = 1:length(toDisplay)
        for lr=1:2
            subplot(length(toDisplay),2,a);
            imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal');
            set(imX,'alphadata',0.7);
            hold on;
            scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),dotSize,toDisplay{d}(:,lr),dotShape,'filled'); axis equal; colorbar;
            ylim([-6 4]);
            
            pcntl = quantile(toDisplay{d}(:),[0.05 0.95]);
            
            caxis([-1 1]*max(abs(pcntl)));
            %         set(s(lr),'markeredgecolor',[1 1 1]*1,'linewidth',0);
            hold off;
            title(['Laser effect ' side{lr} ' ' labels{d}])
            xlim([-5 5]);
            a=a+1;
        end
    end
    %
    cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
        linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
    colormap(cmap);
end

%% nice plot: LvsR and GOvsNG maps 
for s = 1:size(subject,1)
    toDisplay = {[-diff(sites(s).Biases,[],2) log(sum(exp(sites(s).Biases),2))];
                 [-diff(sites(s).CLeft,[],2) nan(size(sites(s).CLeft,1),1)];
                 [-diff(sites(s).CRight,[],2) nan(size(sites(s).CRight,1),1)]};%, sitesP_CLeft, sitesP_CRight};
    labels = {'bias','CL sensitivity','CR sensitivity'};

    dotSize=150;
    dotShape='o';
    kimg=imread('\\basket.cortexlab.net\home\stuff\kirkcaldie_brain_BW.PNG');
    
    a=1;
    figure('name',subject{s,1});
    side={'log(\pi_L/\pi_{R})','log(\pi_{LuR}/\pi_{NG})'};
    for d = 1:length(toDisplay)
        for lr=1:2
            subplot(length(toDisplay),2,a);
            imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal');
            set(imX,'alphadata',0.7);
            hold on;
            scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),dotSize,toDisplay{d}(:,lr),dotShape,'filled'); axis equal; colorbar;
            ylim([-6 4]);
            
            pcntl = quantile(toDisplay{d}(:),[0.05 0.95]);
            
            if lr==1
                caxis([-1 1]*max(abs(pcntl)));
            else
                caxis(log(2*exp(0))+[-1 1]*max(abs(pcntl)));
            end
            %         set(s(lr),'markeredgecolor',[1 1 1]*1,'linewidth',0);
            hold off;
            title(['Laser effect ' side{lr} ' ' labels{d}])
            xlim([-5 5]);
            a=a+1;
        end
    end
    %
    cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
        linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
    colormap(cmap);
end

%% nice plot: performance effects of the laser (no model)
figure;
kimg=imread('\\basket.cortexlab.net\home\stuff\kirkcaldie_brain_BW.PNG');

for s = 1:size(subject,1)
    subplot(size(subject,1),1,s);
    imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal');
    set(imX,'alphadata',0.7);
    hold on;
    
    PF = pivottable(l(s).data.laserIdx,[],l(s).data.feedbackType==1,'mean');
    PF = PF(2:end)-PF(1);
    scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),150,PF,'o','filled'); axis equal; colorbar;
    ylim([-6 4]);
    
    pcntl = quantile(PF,[0.05 0.95]);
    
    caxis([-1 1]*max(abs(pcntl)));
    hold off;
    title(subject{s,1})
    xlim([-5 5]);
    
    cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
        linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
    colormap(cmap);
end


%% nice plot: reaction time effects and lick effects (no model)
for s = 1:size(subject,1)   
    RTs = pivottable(l(s).data.laserIdx,l(s).data.response,l(s).data.RT,'median');
    RTs = bsxfun(@minus,RTs(2:end,:),RTs(1,:));
    RTs = RTs(:,1:2);
    
    Ls = pivottable(l(s).data.laserIdx,l(s).data.response,l(s).data.lickenergy,'median');
    Ls = bsxfun(@minus,Ls(2:end,:),Ls(1,:));
    Ls = Ls(:,1:2);
    
    toDisplay = {PF,RTs,Ls};%, sitesP_CLeft, sitesP_CRight};
    labels = {'RT^{median} - noLaser RT^{median}','lick^{median} - noLaser lick^{median}'};

    dotSize=150;
    dotShape='o';
    kimg=imread('\\basket.cortexlab.net\home\stuff\kirkcaldie_brain_BW.PNG');
    
    a=1;
    figure('name',subject{s,1});
    side={'left choices','right choices'};
    for d = 1:length(toDisplay)
        for lr=1:2
            subplot(length(toDisplay),2,a);
            imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal');
            set(imX,'alphadata',0.7);
            hold on;
            scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),dotSize,toDisplay{d}(:,lr),dotShape,'filled'); axis equal; colorbar;
            ylim([-6 4]);
            
            pcntl = quantile(toDisplay{d}(:),[0.05 0.95]);
            
            caxis([-1 1]*max(abs(pcntl)));
            %         set(s(lr),'markeredgecolor',[1 1 1]*1,'linewidth',0);
            hold off;
            title([side{lr} ' ' labels{d}])
            xlim([-5 5]);
            a=a+1;
        end
    end
    %
    cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
        linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
    colormap(cmap);
end

%% Statistical test of laser effect across mice
%Testing the site parameters against the null hypothesis that they are
%INDIVIDUALLY zero. Therefore need to do ~50 separate tests and so probably
%need to correct for multiple comparisons.

params = {[sites.Biases],[sites.CLeft],[sites.CRight]};

labels = {'bias','CL sensitivity','CR sensitivity'};
dotSize=150;
dotShape='o';
kimg=imread('B:\stuff\kirkcaldie_brain_BW.PNG');

a=1;
figure('name',subject{s,1});
side={'log(\pi_L/\pi_{NG})','log(\pi_R/\pi_{NG})'};
for d = 1:length(params)
    for lr=1:2
        subplot(length(toDisplay),2,a);
        imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal');
        set(imX,'alphadata',0.7);
        hold on;
        
        p = params{d}(:,lr:2:end)';
        [~,ptest,~]=ttest(p);
        
        scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),dotSize,ptest',dotShape,'filled'); axis equal; colorbar;
        ylim([-6 4]);
        caxis([0 0.05]);
        %         set(s(lr),'markeredgecolor',[1 1 1]*1,'linewidth',0);
        hold off;
        title(['Laser effect ' side{lr} ' ' labels{d}])
        xlim([-5 5]);
        a=a+1;
    end
end
colormap gray

%% Split-half reliability (TODO: PER SUBJECT)
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

%% Plot all psychometric curves (select noLaser session and plot all laser additive curves) (TODO: PER SUBJECT)
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

for siteID = 1:size(l.inactivationSite,1)
    paramsL = [sitesP_Biases(siteID,:);
               sitesP_CLeft(siteID,:);
               sitesP_CRight(siteID,:)];
    paramsL = paramsL + paramsNL;
    posIn = l.inactivationSite(siteID,:)/10 + [0.58 0.465];
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

%% (NO GLM) simple map of change in % choice at C=0 over the laser sites, median over sessions (TODO: PER SUBJECT)

for session = 1:length(expRefs)
    noLaser_responses = l.data.response(l.data.laserIdx==0 & diff(l.data.contrast_cond,[],2)==0 & l.data.sessionID==session);
    noLaser_p = arrayfun(@(r)(sum(noLaser_responses==r)/length(noLaser_responses)),1:3);
    
    Laser_p=[];
    for site = 1:size(l.inactivationSite,1)
        Laser_responses = l.data.response(l.data.laserIdx==site & diff(l.data.contrast_cond,[],2)==0 & l.data.sessionID==session);
        Laser_p(site,:) = arrayfun(@(r)(sum(Laser_responses==r)/length(Laser_responses)),1:3);
    end
    
    pDiff(:,:,session) = bsxfun(@minus,Laser_p,noLaser_p);
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