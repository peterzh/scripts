%% Construct list of expRefs
expt = 'sparse_bilateral_2D';

expRefs = {
    'Spemann',...
    LaserOmnibusCheckSessions('Spemann',expt);
    
    'Murphy',...
    LaserOmnibusCheckSessions('Murphy',expt);
    
        'Morgan',...
        LaserOmnibusCheckSessions('Morgan',expt);
    
    'Whipple',...
    LaserOmnibusCheckSessions('Whipple',expt)
    };

numSubjects = size(expRefs,1);
clear shuffle; disp('done');

%% Load all sessions.
clear l;

c50_fits = cell(1,size(expRefs,1));
cn_fits = cell(1,size(expRefs,1));

for s = 1:numSubjects
    name = expRefs{s,1};
    
    %Combine all sessions
    D=struct;
    %i=1;
    
    for session = 1:length(expRefs{s,2})
        disp([num2str(session) '/' num2str(length(expRefs{s,2}))]);
        %try
        d = laserGLM(expRefs{s,2}{session}).data;
        d = structfun(@(x)(x(6:(end-14),:)),d,'uni',0); %trim first 5 trials and last 15
        if length(d.response) < 100
            error('not enough trials');
        end
        
        e = getrow(d,d.laserIdx==0);
        g=GLM(e).setModel('C50-subset').fit;
        %             figure;g.plotFit;
        c50_fits{s}(session,:) = g.parameterFits(5:6);
        %             c50sub.n(i) = g.parameterFits(5);
        %             c50sub.c50(i) = g.parameterFits(6);
        
        g=GLM(e).setModel('C^N-subset').fit;
        cn_fits{s}(session,:) = g.parameterFits(5);
        
        %         d.sessionID = ones(length(d.response),1)*s;
        d.sessionID = ones(length(d.response),1)*session;
        %i=i+1;
        D = addstruct(D,d);
        %catch
        %   warning(eRefs{session});
        %end
    end
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
X = cell(1,numSubjects);
Y = cell(1,numSubjects);
figure;

for s = 1:numSubjects
    numSessions = max(l(s).data.sessionID);
    numSites = max(l(s).data.laserIdx);
    numTrials = size(l(s).data.response,1);
    
    Y{s} = l(s).data.response;
    
    %construct design matrix X (for loop for ease of understanding not speed!)
    X{s} = sparse(numTrials,3*numSessions + 3*numSites);
    
    contrast_rep_fcn = @(c,n,c50)((c.^n)/(c.^n + c50^n));
    for i = 1:numTrials
        
        thisSession = l(s).data.sessionID(i);
        cfn = @(c)contrast_rep_fcn(c,c50_fits{s}(thisSession,1),c50_fits{s}(thisSession,2));
        X{s}(i,3*thisSession - 2) = 1;
        X{s}(i,3*thisSession - 1) = cfn(l(s).data.contrast_cond(i,1));
        X{s}(i,3*thisSession - 0) = cfn(l(s).data.contrast_cond(i,2));
        
        thisSite = l(s).data.laserIdx(i);
        if thisSite ~= 0
            X{s}(i,3*numSessions + 3*thisSite - 2) = 1;
            X{s}(i,3*numSessions + 3*thisSite - 1) = cfn(l(s).data.contrast_cond(i,1));
            X{s}(i,3*numSessions + 3*thisSite - 0) = cfn(l(s).data.contrast_cond(i,2));
            
        end
    end
    
    opts=struct;
        opts.intr=0; %don't add a global intercept
%     opts.intr=1;
    % opts.alpha=1; %lasso regularisation (L1)
    %     opts.alpha=0; %ridge (L2)
    opts.alpha=0.5; %elasticnet (both)
    
    % penalty = ones(1,size(X,2)); penalty(1:3:end)=0;
    % opts.penalty_factor=penalty; %Try penalising only the contrast terms
        fit=cvglmnet(X{s},Y{s},'multinomial',glmnetSet(opts));
%     fit=glmnet(X{s},Y{s},'multinomial',glmnetSet(opts));
    %     cvglmnetPlot(fit);
        b=cvglmnetCoef(fit,0.001);
%     b=glmnetCoef(fit,0.01);
    
        disp([expRefs{s,1} ' cv lambda min: ' num2str(fit.lambda_min)]);
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
    %     figure('name',expRefs{s,1});
    subplot(numSubjects,2,2*s-1);
    imagesc(X{s});
    hold on;
    pos=0;
    for session = 1:numSessions
        numTr = sum(l(s).data.sessionID==session);
        tx=text(25,pos + 0.5*numTr,[expRefs{s,2}{session} ' n=' num2str(numTr)],'interpreter','none');
        tx.Color=[1 1 1];
        pos = pos + numTr;
    end
    hold off;
    title(expRefs{s,1});
    
    subplot(numSubjects,2,2*s);
    tab = tabulate(l(s).data.laserIdx);
    tab = tab(2:end,2);
    scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),200,tab,'o','filled'); axis equal; colorbar; colormap(gca,'parula')
    axis([-5 5 -4 3]);
    % set(s,'markeredgecolor',[0.8 0.8 0.8]);
    title('Number of trials at each site');
end

%% Plot model parameters over sessions and sites
for s = 1:numSubjects
    name = expRefs{s,1};
    %Plot the session by session values of the parameters to see whether they
    %change much
    figure('name',name);
    subplot(3,2,1); bar(sess(s).Biases); title('Biases for each session');
    subplot(3,2,3); bar(sess(s).CLeft); title('CL scaling for each session');
    subplot(3,2,5); bar(sess(s).CRight); title('CR scaling for each session');
    set(gca,'XTick',1:length(expRefs{s,2}),'XTickLabel',expRefs{s,2},'XTickLabelRotation',90);
    
    %Plot the site by site values of the parameters to see whether they
    %change much
    subplot(3,2,2); bar(sites(s).Biases); title('Biases for each site');
    subplot(3,2,4); bar(sites(s).CLeft); title('CL scaling for each site');
    subplot(3,2,6); bar(sites(s).CRight); title('CR scaling for each site');
end

%% LASER EFFECT MAP: LvsNG and RvsNG maps
for s = 1:numSubjects
    toDisplay = {sites(s).Biases,sites(s).CLeft,sites(s).CRight};%, sitesP_CLeft, sitesP_CRight};
    labels = {'bias','CL sensitivity','CR sensitivity'};
    
    if exist('shuffle','var') %if shuffle analysis done
        null_quarts_95 = {quantile(shuffle(s).Bias,[0.025 0.975],3),...
            quantile(shuffle(s).CLeft,[0.025 0.975],3),...
            quantile(shuffle(s).CRight,[0.025 0.975],3)};
        
        null_quarts_99 = {quantile(shuffle(s).Bias,[0.005 0.995],3),...
            quantile(shuffle(s).CLeft,[0.005 0.995],3),...
            quantile(shuffle(s).CRight,[0.005 0.995],3)};
%     elseif exist('Qshuffle','var')
%         
    end
    
    
    dotSize=150;
    dotShape='o';
    kimg=imread('D:\kirkcaldie_brain_BW.PNG');
    
    a=1;
    figure('name',expRefs{s,1});
    side={'log(\pi_L/\pi_{NG})','log(\pi_R/\pi_{NG})'};
    for d = 1:length(toDisplay)
        for lr=1:2
            subplot(length(toDisplay),2,a);
            imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal');
            set(imX,'alphadata',0.7);
            hold on;
            
            if exist('shuffle','var') %if shuffle analysis done
                sig95 = double(toDisplay{d}(:,lr) < null_quarts_95{d}(:,lr,1) | toDisplay{d}(:,lr) > null_quarts_95{d}(:,lr,2));
                sig99 = double(toDisplay{d}(:,lr) < null_quarts_99{d}(:,lr,1) | toDisplay{d}(:,lr) > null_quarts_99{d}(:,lr,2));
                shuffle(s).sig95(:,lr,d) = sig95;
                shuffle(s).sig99(:,lr,d) = sig99;
                sig95(sig95==0)=0.1;
                
                sig99(sig99==0)=0.1;
                %                 scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),sig*15,'ok','filled');
                scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),sig99*150,toDisplay{d}(:,lr),dotShape,'filled'); axis equal; colorbar;
                scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),sig95*150*0.6,toDisplay{d}(:,lr),dotShape,'filled'); axis equal; colorbar;
                
            else
                scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),dotSize,toDisplay{d}(:,lr),dotShape,'filled'); axis equal; colorbar;
                
            end
            ylim([-6 4]);
            
            %             pcntl = quantile(toDisplay{d}(:),[0.05 0.95]);
            %             caxis([-1 1]*max(abs(pcntl)));
            caxis([-1 1]*max(abs(toDisplay{d}(:))));
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

%% LASER EFFECT MAP: LvsR and GOvsNG maps
for s = 1:numSubjects
    toDisplay = {[-diff(sites(s).Biases,[],2) log(sum(exp(sites(s).Biases),2))];
        [-diff(sites(s).CLeft,[],2) nan(size(sites(s).CLeft,1),1)];
        [-diff(sites(s).CRight,[],2) nan(size(sites(s).CRight,1),1)]};%, sitesP_CLeft, sitesP_CRight};
    labels = {'bias','CL sensitivity','CR sensitivity'};
    
    dotSize=150;
    dotShape='o';
    kimg=imread('\\basket.cortexlab.net\home\stuff\kirkcaldie_brain_BW.PNG');
    
    a=1;
    figure('name',expRefs{s,1});
    side={'log(\pi_L/\pi_{R})','log(\pi_{LuR}/\pi_{NG})'};
    for d = 1:length(toDisplay)
        for lr=1:2
            subplot(length(toDisplay),2,a);
            imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal');
            set(imX,'alphadata',0.7);
            hold on;
            scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),dotSize,toDisplay{d}(:,lr),dotShape,'filled'); axis equal; colorbar;
            ylim([-6 4]);
            
            %             pcntl = quantile(toDisplay{d}(:),[0.05 0.95]);
            %             caxis([-1 1]*max(abs(pcntl)));
            caxis([-1 1]*max(abs(toDisplay{d}(:))));
            %             if lr==1
            %                 caxis([-1 1]*max(abs(pcntl)));
            %             else
            %                 caxis(log(2*exp(0))+[-1 1]*max(abs(pcntl)));
            %             end
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

%% LASER EFFECT MAP: performance effects of the laser (no model)
figure;
kimg=imread('\\basket.cortexlab.net\home\stuff\kirkcaldie_brain_BW.PNG');

for s = 1:numSubjects
    subplot(numSubjects,1,s);
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
    title(expRefs{s,1})
    xlim([-5 5]);
    
    cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
        linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
    colormap(cmap);
end


%% LASER EFFECT MAP: reaction time effects and lick effects (no model)
for s = 1:numSubjects
    PF = pivottable(l(s).data.laserIdx,[],l(s).data.feedbackType==1,'mean');
    PF = PF(2:end)-PF(1);
    PF = [PF nan(length(PF),1)];
    
    RTs = pivottable(l(s).data.laserIdx,l(s).data.response,l(s).data.RT,'median');
    RTs = bsxfun(@minus,RTs(2:end,:),RTs(1,:));
    RTs = RTs(:,1:2);
    
%     Ls = pivottable(l(s).data.laserIdx,l(s).data.response,l(s).data.lickenergy,'median');
%     Ls = bsxfun(@minus,Ls(2:end,:),Ls(1,:));
%     Ls = Ls(:,1:2);
%     
    toDisplay = {RTs,PF};%, sitesP_CLeft, sitesP_CRight};
    labels = {'RT^{median} - noLaser RT^{median}','perf^{median} - noLaser perf^{median}'};
    
    dotSize=150;
    dotShape='o';
    kimg=imread('D:\kirkcaldie_brain_BW.PNG');
    
    a=1;
    figure('name',expRefs{s,1});
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

%% RELIABILITY: Split half reliability
for s = 1:numSubjects
    numTrials = length(Y{s});
    C = cvpartition(numTrials,'KFold',2);
    
    sitesP=[];
    sitesP_Biases_split=[];
    sitesP_CLeft_split=[];
    sitesP_CRight_split=[];
    for split = 1:C.NumTestSets
        Ysplit = Y{s}(C.training(split));
        Xsplit = X{s}(C.training(split),:);
        fit=cvglmnet(Xsplit,Ysplit,'multinomial',glmnetSet(opts));
        b=cvglmnetCoef(fit);
        b=[b{1}-b{3} b{2}-b{3}];
        b(1,:) = [];
        
        numSessions = length(expRefs{s,2});
        
        sitesP = b(3*numSessions+1:end,:);
        sitesP_Biases_split(:,:,split) = sitesP(1:3:end,:);
        sitesP_CLeft_split(:,:,split) = sitesP(2:3:end,:);
        sitesP_CRight_split(:,:,split) = sitesP(3:3:end,:);
    end
    
    %     figure('name',expRefs{s,1});
    %     h(1)=subplot(3,1,1);
    %     Diff = std(sitesP_Biases_split,[],3);
    %     hist(Diff); title('std in biases for all sites');
    %     h(2)=subplot(3,1,2);
    %     Diff = std(sitesP_CLeft_split,[],3);
    %     hist(Diff); title('std in left contrast sens for all sites');
    %     h(3)=subplot(3,1,3);
    %     Diff = std(sitesP_CRight_split,[],3);
    %     hist(Diff); title('std in right contrast sens for all sites');
    %     % linkaxes(h,'x');
    
    figure('name',expRefs{s,1})
    subplot(3,2,1);
    Diff = std(sitesP_Biases_split,[],3);
    sx(1)=scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),300,Diff(:,1),'s','filled'); axis equal; colorbar;
    subplot(3,2,2);
    sx(2)=scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),300,Diff(:,2),'s','filled'); axis equal; colorbar;
    
    subplot(3,2,3);
    Diff = std(sitesP_CLeft_split,[],3);
    sx(3)=scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),300,Diff(:,1),'s','filled'); axis equal; colorbar;
    subplot(3,2,4);
    sx(4)=scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),300,Diff(:,2),'s','filled'); axis equal; colorbar;
    
    subplot(3,2,5);
    Diff = std(sitesP_CRight_split,[],3);
    sx(5)=scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),300,Diff(:,1),'s','filled'); axis equal; colorbar;
    subplot(3,2,6);
    sx(6)=scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),300,Diff(:,2),'s','filled'); axis equal; colorbar;
    
    cmap = colormap(gray);
    colormap(flipud(cmap));
    set(sx,'markeredgecolor',[0 0 0]);
    %     set(get(gcf,'children'),'Xlim',[-5 5],'Ylim',[-6 4])
end

%% RELIABILITY: Shuffle analysis
%Totally randomise the labelling of each trial to a laser site. Under this
%shuffling there should be no true laser effect. Thus model fits will give
%a null distribution
NUM_SHUFFLES = 200;
shuffle = struct;
for s = 1:numSubjects
    numSessions = max(l(s).data.sessionID);
    numSites = max(l(s).data.laserIdx);
    numTrials = size(l(s).data.response,1);
    
    Xsess = sparse(numTrials,3*numSessions);
    
    %CONSTRUCT SESSION PORTION OF DESIGN MATRIX
    contrast_rep_fcn = @(c,n,c50)((c.^n)/(c.^n + c50^n));
    for i = 1:numTrials
        thisSession = l(s).data.sessionID(i);
        cfn = @(c)contrast_rep_fcn(c,c50_fits{s}(thisSession,1),c50_fits{s}(thisSession,2));
        Xsess(i,3*thisSession - 2) = 1;
        Xsess(i,3*thisSession - 1) = cfn(l(s).data.contrast_cond(i,1));
        Xsess(i,3*thisSession - 0) = cfn(l(s).data.contrast_cond(i,2));
    end
    
    for shuff = 1:NUM_SHUFFLES
        disp([num2str(shuff) '/' num2str(NUM_SHUFFLES)]);
        %CONSTRUCT SITE PORTION OF DESIGN MATRIX MANY DIFFERENT TIMES
        Xsites = sparse(numTrials,3*numSites);
        siteIDs = l(s).data.laserIdx(randperm(numTrials));
        for i=1:numTrials
            thisSite = siteIDs(i);
            if thisSite > 0
                Xsites(i,3*thisSite - 2) = 1;
                Xsites(i,3*thisSite - 1) = cfn(l(s).data.contrast_cond(i,1));
                Xsites(i,3*thisSite - 0) = cfn(l(s).data.contrast_cond(i,2));
            end
        end
        
        fit=glmnet([Xsess Xsites],Y{s},'multinomial',glmnetSet(opts));
        
        %                 fit=cvglmnet([Xsess Xsites],Y,'multinomial',glmnetSet(opts));
        % b=cvglmnetCoef(fit,0.01);
        b=glmnetCoef(fit,0.01);
        
        b=[b{1}-b{3} b{2}-b{3}];
        b(1,:) = [];
        
        shuffle(s).sitesP(:,:,shuff) = b(3*numSessions+1:end,:);
        shuffle(s).Bias(:,:,shuff) = shuffle(s).sitesP(1:3:end,:,shuff);
        shuffle(s).CLeft(:,:,shuff) = shuffle(s).sitesP(2:3:end,:,shuff);
        shuffle(s).CRight(:,:,shuff) = shuffle(s).sitesP(3:3:end,:,shuff);
        
    end
    
    PLOTTING_VARS = {shuffle(s).Bias, shuffle(s).CLeft, shuffle(s).CRight};
    PLOTTING_VARS_labels = {'Bias','CLeft','CRight'};
    for plotVar=1:3
        figure('name',[expRefs{s,1} '_' PLOTTING_VARS_labels{plotVar}]);
        kimg=imread('D:\kirkcaldie_brain_BW.PNG');
        imX=image(-5:1:5,4:-1:-6,kimg);
        set(gca,'ydir','normal');
        set(gca,'Position',[0 0 1 1]);
        
        for loc = 1:size(l(s).inactivationSite,1)
            posIn = l(s).inactivationSite(loc,1:2)/10 + [0.58 0.465];
            
            null = squeeze(PLOTTING_VARS{plotVar}(loc,:,:))';
            axes;
            scatter(null(:,1),null(:,2),'b.');
            %             boxplot(null,'whisker',0,'symbol','');
            
            hold on;
            fitted = sites(s).Biases(loc,:);
            scatter(fitted(1),fitted(2),'r+'); plot(0,0,'ko');
            %             plot([1 2],fitted,'mo','markersize',5,'markerfacecolor','m')
            
            hold off; axis auto;
            set(gca,'Position',[posIn(2) posIn(1) 0.07 0.07],'box','off','XColor',[0 0 0 0],'YColor',[0 0 0 0]);
        end
    end
    
end


%% RELIABILITY: Fake inactivation location
%adds another laser site from a certain proportion of NoLaser trials. Since
%we know the laser has no true effect here, the estimated parameters from
%the fit can give an indication of the null distribution under this
%behavioural model.
Qshuffle = struct;
proportion = 0.1;
NUM_SHUFFLES = 300;
for s = 1:numSubjects
    nonL_idx = find(l(s).data.laserIdx==0);
    for shuff = 1:NUM_SHUFFLES
        disp(shuff);
        nL = randsample(nonL_idx,round(proportion*length(nonL_idx)));
        QX = zeros(size(X{s},1),3);
        
        for t = 1:length(nL)
            trial = nL(t);
            c = l(s).data.contrast_cond(trial,:);
            sessionID = l(s).data.sessionID(trial);
            cfn = @(c)contrast_rep_fcn(c,c50_fits{s}(sessionID,1),c50_fits{s}(sessionID,2));
            
            QX(trial,1) = 1;
            QX(trial,2) = cfn(c(1));
            QX(trial,3) = cfn(c(2));
        end
        
        %         QX = rand(size(X{s},1),3);
        %         QX = bsxfun(@minus,QX,mean(QX));
        
        fit=glmnet([X{s} QX],Y{s},'multinomial',glmnetSet(opts));
        %         fit=glmnet([QX],Y{s},'multinomial',glmnetSet(opts));
        
        b=glmnetCoef(fit,0.01);
        b=[b{1}-b{3} b{2}-b{3}];
        b(1,:) = [];
        
        Qshuffle(s).Bias(shuff,:) = b(end-2,:);
        Qshuffle(s).CLeft(shuff,:) = b(end-1,:);
        Qshuffle(s).CRight(shuff,:) = b(end,:);
    end
end

%% Check whether model is fitting well on non-laser trials
for s = 1:numSubjects
    figure('name',expRefs{s,1});
    numSessions = max(l(s).data.sessionID);
    for session = 1:numSessions
        D = getrow(l(s).data,l(s).data.laserIdx==0 & l(s).data.sessionID==session);
        params = [sess(s).Biases(session,1),...
            sess(s).CLeft(session,1),...
            sess(s).CRight(session,1),...
            sess(s).Biases(session,2),...
            sess(s).CLeft(session,2),...
            sess(s).CRight(session,2),...
            c50_fits{s}(session,:)];
        
        D.p_hats = g.calculatePhat(params,D.contrast_cond);
        
        stim_labels = {'C=0','CL=CR','C=L','C=R'};
        stim = cell(1,4);
        stim{1} = sum(D.contrast_cond,2)==0;
        stim{2} = ~stim{1} & (D.contrast_cond(:,1)==D.contrast_cond(:,2));
        stim{3} = diff(D.contrast_cond,[],2)<0;
        stim{4} = diff(D.contrast_cond,[],2)>0;
        
        for stim_type = 1:4
            r = D.response(stim{stim_type});
            tab = [sum(r==1) sum(r==3) sum(r==2)]';
            [p,pci] = binofit(tab,sum(tab));
            
            subplot(numSessions,4,4*session - 4 + stim_type);
            errorbar(1:3,p,p-pci(:,1),pci(:,2)-p,'rs'); %plot actual
            hold on;
            plot(mean(D.p_hats(stim{stim_type},[1 3 2])),'.:','markersize',12);
            
            hold off;
            
            if session==1
                title(stim_labels{stim_type});
            end
        end
    end
    
    ax=get(gcf,'children');
    set(ax,'Xticklabel','');
    set(ax(1),'Xticklabel',{'L','NG','R'},'xtick',1:3);
end

%% Non-laser bootstrapping of parameters
numiter = 100;
for s = 1:numSubjects
    numSessions = length(expRefs{s,2});
    Xb = X{s};
    
    bootstrap = struct;
    for iter = 1:numiter
        for session = 1:numSessions
            sidx = l(s).data.sessionID==session;
            ridx = randsample(find(sidx),sum(sidx),true);
            Xb(sidx,:) = Xb(ridx,:);
        end
        
        %fit
        fit=glmnet(Xb,Y{s},'multinomial',glmnetSet(opts));
        b=glmnetCoef(fit,0.001);
        b=[b{1}-b{3} b{2}-b{3}];
        b(1,:) = [];
        
        sessionP_b = b(1:3*numSessions,:);
        bootstrap.Biases(:,:,iter) = sessionP_b(1:3:end,:);
        bootstrap.CLeft(:,:,iter) = sessionP_b(2:3:end,:);
        bootstrap.CRight(:,:,iter) = sessionP_b(3:3:end,:);
    end
    
    biases = permute(bootstrap.Biases,[3 2 1]);
    sensL = permute(bootstrap.CLeft,[3 2 1]);
    sensR = permute(bootstrap.CRight,[3 2 1]);
    labels = {'b_{LvNG}','sL_{LvNG}','sR_{LvNG}','sR_{RvNG}','sL_{RvNG}','b_{RvNG}'};
    for session = 1:numSessions
        figure('name',[expRefs{s,1} ' session ' num2str(session)]);
        gplotmatrix([biases(:,1,session) sensL(:,1,session) sensR(:,1,session) sensR(:,2,session) sensL(:,2,session) biases(:,2,session)],[],[],'b','.',5,'on','none',labels);
    end
end

%% LASER EFFECT MAP: psychometric curves)
% Only possible for 1D contrast trials since the contrast can be projected
% onto 1 axis.
labels = {'Choose L','Choose R','NG'};
for s = 1:numSubjects
    %     sessionID = 1;
    g = g.setModel('C50');
    %     testC = l(s).data.contrast_cond;
    
    evalCL = linspace(0,max(l(s).data.contrast_cond(:,1)),100);
    evalCR = linspace(0,max(l(s).data.contrast_cond(:,1)),100);
    
    %     paramsNL = [sess(s).Biases(sessionID,:);
    %         sess(s).CLeft(sessionID,:);
    %         sess(s).CRight(sessionID,:)];
    paramsNL = [median(sess(s).Biases,1);
        median(sess(s).CLeft,1);
        median(sess(s).CRight,1)];
    %     shape = c50_fits{s}(sessionID,:);
    shape = median(c50_fits{s},1);
    
    
    %     paramsNL = [-1 -1;3 0; 0 3];
    %     shape = [2 0.1];
    
    prop_NL=nan(length(evalCL),length(evalCR),3);
    for cl = 1:length(evalCL)
        for cr = 1:length(evalCR)
            p = g.calculatePhat([paramsNL(:);shape'],[evalCL(cl) evalCR(cr)]);
            for i=1:3
                prop_NL(cl,cr,i) = p(i);
            end
        end
    end
    
    for plotted_choice = 1:3
        
        figure('name',expRefs{s,1});
        kimg=imread('D:\kirkcaldie_brain_BW.PNG');
        imX=image(-5:1:5,4:-1:-6,kimg); set(gca,'ydir','normal');
        set(gca,'Position',[0 0 1 1]);
        title(labels{plotted_choice});
        text(-4.5,3,sprintf('%s\n%s\n%s',expRefs{s,1},expt,labels{plotted_choice}),'fontsize',20)
        
        for siteID = 1:size(l(s).inactivationSite,1)
            %             pass = 1;
            %             if exist('shuffle','var') %if shuffle analysis done
            %                 if ~any(any(squeeze(shuffle(s).sig95(siteID,:,:))))
            %                     pass = 0;
            %                 end
            %             end
            pass = 1;
            
            if pass==1
                paramsL = [sites(s).Biases(siteID,:);
                    sites(s).CLeft(siteID,:);
                    sites(s).CRight(siteID,:)];
                posIn = l(s).inactivationSite(siteID,1:2)/10 + [0.58 0.465];
                axes;
                
                hold on;
                
                prop_L=nan(length(evalCL),length(evalCR),3);
                for cl = 1:length(evalCL)
                    for cr = 1:length(evalCR)
                        p = g.calculatePhat([paramsL(:)+paramsNL(:);shape'],[evalCL(cl) evalCR(cr)]);
                        for i=1:3
                            prop_L(cl,cr,i) = p(i);
                        end
                    end
                end
                %             imagesc(evalCR,evalCL,sum(prop_L(:,:,1:2),3)-sum(prop_NL(:,:,1:2),3));
                %             set(gca,'YDir','normal'); axis([0 max(evalCR) 0 max(evalCL)]);
                
                num_contours = 6;
                contour(evalCL,evalCR,prop_NL(:,:,plotted_choice),num_contours,'linewidth',1,'linestyle','--'); caxis([0 1]);
                hold on;
                contour(evalCL,evalCR,prop_L(:,:,plotted_choice),num_contours,'linewidth',1,'linestyle','-'); caxis([0 1]);
                
                xlabel('CR');ylabel('CL');
                
                caxis([0 1]);
                %             caxis([-1 1]);
                %             cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
                %                 linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
                %             colormap(cmap);
                
                set(gca,'Position',[posIn(2) posIn(1) 0.07 0.07],'box','off','XColor',[0 0 0 1],'YColor',[0 0 0 1]);
            end
        end
        % set(gcf,'color','white');
        set(get(gcf,'children'),'fontsize',5)
    end
end

%% Laser GLM: overall crossvalidated goodness of fit
cv_gof = [];
for s = 1:numSubjects
    cv = cvpartition(size(X{s},1),'kfold',10);
    
    p_hats = nan(cv.NumObservations,1);
    for iter = 1:cv.NumTestSets
        trainX = X{s}(cv.training(iter),:);
        trainY = Y{s}(cv.training(iter));
        
        testX = X{s}(cv.test(iter),:);
        testY = Y{s}(cv.test(iter));
        
        fit=glmnet(trainX,trainY,'multinomial',glmnetSet(opts));
        p = glmnetPredict(fit,testX,0.001,'response');
        p_hats(cv.test(iter)) = p(sub2ind(size(p), [1:length(testY)]', testY));
    end
    
    tab=tabulate(Y{s}); tab=tab(:,3)/100;
    guess_bpt = sum(tab.*log2(tab));
    loglik_bpt = mean(log2(p_hats));
    
    cv_gof(s) = loglik_bpt-guess_bpt;
end

bar(cv_gof); ylabel('loglik [bits] relative to guessing @~-1.5');
set(gca,'XTickLabel',expRefs(:,1));

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

%% (NO GLM) RT across stimulus conditions on noLaser trials
figure; labels={'left','right'};
for s = 1:numSubjects
    
    con = unique(l(s).data.contrast_cond(:));
    
    for resp=1:2
        rt = [];
        for cl=1:length(con)
            for cr=1:length(con)
                idx = l(s).data.contrast_cond(:,1)==con(cl) & l(s).data.contrast_cond(:,2)==con(cr) & l(s).data.response==resp & l(s).data.laserIdx==0;
                rt(cl,cr) = mean(l(s).data.RT(idx));
            end
        end
        rt(1,1)=nan;
        subplot(numSubjects,2,2*s-1+resp-1);
        imagesc(con,con,rt);
        set(gca,'ydir','normal'); xlabel('CR'); ylabel('CL'); colorbar; axis square
        title([expRefs{s,1} ' chose ' labels{resp}]);
    end
end