%% 
expt = 'sparse_unilateral_2D';

%% Construct + save list of expRefs
expRefs = {
    'Spemann',...
    LaserOmnibusCheckSessions('Spemann',expt);
    
    'Murphy',...
    LaserOmnibusCheckSessions('Murphy',expt);
    
%         'Morgan',...
%         LaserOmnibusCheckSessions('Morgan',expt);
%     
    'Whipple',...
    LaserOmnibusCheckSessions('Whipple',expt)
    };

save(['\\basket.cortexlab.net\home\omnibus_files\' expt '.mat'],'expRefs');
disp('done');

%% Load all sessions.
clear l;
clear Qshuffle;
load(['\\basket.cortexlab.net\home\omnibus_files\' expt '.mat']);

MOREDATA=0;

numSubjects = size(expRefs,1);
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
        d = laserGLM(expRefs{s,2}{session});
        
        if MOREDATA==1
            d = d.addData('lick');
            d = d.addData('pupil');
        end
        
        d = d.data;
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

%% Fix all C50 and N, and ^N parameters constant for each mouse
c50_fits=cellfun(@(a)(repmat(mean(a,1),size(a,1),1)),c50_fits,'uni',0);
cn_fits=cellfun(@(a)(repmat(mean(a,1),size(a,1),1)),cn_fits,'uni',0);

%% Define and fit model
sites = struct('Biases',{},'CLeft',{},'CRight',{},'laserSite',{});
sess = struct('Biases',{},'CLeft',{},'CRight',{});

sites2 = struct('Bias',{},'Sens',{},'laserSite',{});
sess2 = struct('Bias',{},'Sens',{});

X = cell(1,numSubjects);
Y = cell(1,numSubjects);
figure;

fitMode = 1;
for s = 1:numSubjects
    numSessions = max(l(s).data.sessionID);
    numSites = max(l(s).data.laserIdx);
    numTrials = size(l(s).data.response,1);
    
    Y{s} = l(s).data.response;
    
    %construct design matrix X (for loop for ease of understanding not speed!)
    X{s} = sparse(numTrials,3*numSessions + 3*numSites);
    X_LvNG{s} = sparse(numTrials,2*numSessions + 2*numSites);
    X_RvNG{s} = sparse(numTrials,2*numSessions + 2*numSites);

    contrast_rep_fcn = @(c,n,c50)((c.^n)/(c.^n + c50^n));
    for i = 1:numTrials
        
        thisSession = l(s).data.sessionID(i);
        cfn = @(c)contrast_rep_fcn(c,c50_fits{s}(thisSession,1),c50_fits{s}(thisSession,2));
        X{s}(i,3*thisSession - 2) = 1;
        X{s}(i,3*thisSession - 1) = cfn(l(s).data.contrast_cond(i,1));
        X{s}(i,3*thisSession - 0) = cfn(l(s).data.contrast_cond(i,2));
        
        X_LvNG{s}(i,2*thisSession - 1) = 1;
        X_LvNG{s}(i,2*thisSession - 0) = cfn(l(s).data.contrast_cond(i,1));
        
        X_RvNG{s}(i,2*thisSession - 1) = 1;
        X_RvNG{s}(i,2*thisSession - 0) = cfn(l(s).data.contrast_cond(i,2));
        
        thisSite = l(s).data.laserIdx(i);
        if thisSite ~= 0
            X{s}(i,3*numSessions + 3*thisSite - 2) = 1;
            X{s}(i,3*numSessions + 3*thisSite - 1) = cfn(l(s).data.contrast_cond(i,1));
            X{s}(i,3*numSessions + 3*thisSite - 0) = cfn(l(s).data.contrast_cond(i,2));
            
            X_LvNG{s}(i,2*numSessions + 2*thisSite - 1) = 1;
            X_LvNG{s}(i,2*numSessions + 2*thisSite - 0) = cfn(l(s).data.contrast_cond(i,1));
            
            X_RvNG{s}(i,2*numSessions + 2*thisSite - 1) = 1;
            X_RvNG{s}(i,2*numSessions + 2*thisSite - 0) = cfn(l(s).data.contrast_cond(i,2));
        end
    end
       
    opts=struct;
    opts.intr=1; %global intercept
    
    fit=cvglmnet(X{s},Y{s},'multinomial',glmnetSet(opts));
    b=cvglmnetCoef(fit,0);
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
    
    %     
%     %Also try separate LvNG and RvNG fitting
    %Remove irrelevant trials for the LvNG and RvNG matrices (for separate
    %fitting method)
    Y_LvNG{s} = Y{s}(Y{s}~=2); Y_LvNG{s}(Y_LvNG{s}==3)=0; %recode NOGO to be 0
    Y_RvNG{s} = Y{s}(Y{s}~=1); Y_RvNG{s}(Y_RvNG{s}==3)=0; Y_RvNG{s}(Y_RvNG{s}==2)=1;
    X_LvNG{s}(Y{s}==2,:) = [];
    X_RvNG{s}(Y{s}==1,:) = [];

    b2=[];
    fit=glmnet(X_LvNG{s}, Y_LvNG{s}, 'binomial',glmnetSet(opts));
    b2(:,1)=glmnetCoef(fit,0.001);
    fit=glmnet(X_RvNG{s}, Y_RvNG{s}, 'binomial',glmnetSet(opts));
    b2(:,2)=glmnetCoef(fit,0.001);
    b2(1,:) = [];
    
    sessionP2 = b2(1:2*numSessions,:);
    sess2(s).Bias = sessionP2(1:2:end,:);
    sess2(s).Sens = sessionP2(2:2:end,:);
    
    sitesP2 = b2(2*numSessions+1:end,:);
    sites2(s).Bias = sitesP2(1:2:end,:);
    sites2(s).Sens = sitesP2(2:2:end,:);
    sites2(s).laserSite = l(s).inactivationSite;
    
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

%% Two-stage fitting alternative
%Fits nonLaser parameters first and then the laser effects
fitMode = 2;
for s = 1:numSubjects
    twoStageFitOpts = struct;
    twoStageFitOpts.intr = 1;
    
    numSessions = length(expRefs{s,2});
    
    %First fit non-laser portion of data
    nonL_idx = l(s).data.laserIdx==0;
    XnL = X{s}(nonL_idx,:);
    YnL = Y{s}(nonL_idx);
    fit1=cvglmnet(XnL(:,1:(3*numSessions)),YnL,'multinomial',glmnetSet(twoStageFitOpts));
    disp([expRefs{s,1} ' nonLaser fit lambda min: ' num2str(fit1.lambda_min)]);

    %pull out parameter fits
    b=cvglmnetCoef(fit1);
    b=[b{1}-b{3} b{2}-b{3}];
    b(1,:) = [];
    sess(s).Biases = b(1:3:end,:);
    sess(s).CLeft = b(2:3:end,:);
    sess(s).CRight = b(3:3:end,:);
    
    % %Then fit laser portion of data, using nonLaser params as offsets
    XL = X{s}(l(s).data.laserIdx>0,:);
    YL = Y{s}(l(s).data.laserIdx>0);
    twoStageFitOpts.intr = 0;
    twoStageFitOpts.offset = cvglmnetPredict(fit1,XL(:,1:(3*numSessions)),[],'link');
    fit2=cvglmnet(XL(:,((3*numSessions)+1):end),YL,'multinomial',glmnetSet(twoStageFitOpts));
    disp([expRefs{s,1} ' Laser fit lambda min: ' num2str(fit2.lambda_min)]);

    %pull out parameter fits
    b=cvglmnetCoef(fit2);
    b=[b{1}-b{3} b{2}-b{3}];
    b(1,:) = [];
    sites(s).Biases = b(1:3:end,:);
    sites(s).CLeft = b(2:3:end,:);
    sites(s).CRight = b(3:3:end,:);    
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
            
            if exist('Qshuffle','var') && fitMode == 2 %if shuffle analysis done using 2 stage fitting
                ci_95 = Qshuffle(s).CI95(:,lr,d);
                ci_99 = Qshuffle(s).CI99(:,lr,d);
                
                sig95 = double(toDisplay{d}(:,lr) < ci_95(1) | ci_95(2) < toDisplay{d}(:,lr));
                sig99 = double(toDisplay{d}(:,lr) < ci_99(1) | ci_99(2) < toDisplay{d}(:,lr));
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
            set(gca,'xtick','','ytick','');
        set(gca,'box','off','XColor','w','YColor','w');
        end
    end
    %
    cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
        linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
    colormap(cmap);
    
    set(gcf,'Color','w');
end

%% LASER EFFECT MAP: LvsNG and RvsNG maps CONDENSED OVER SUBJECTS (USING NON-CROSS MODEL!)
figure;
for s = 1:numSubjects
    labels = {'bias','CL sensitivity','CR sensitivity'};

    dotSize=150;
    dotShape='o';
    kimg=imread('D:\kirkcaldie_brain_BW.PNG');
    
    subplot(4,numSubjects,s);
    imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal'); set(imX,'alphadata',0.7);
    hold on;
    scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),dotSize,sites2(s).Bias(:,1),dotShape,'filled'); axis equal; colorbar;
    axis([-5 5 -6 4]); set(gca,'box','off','Xtick','','YTick',''); if s==1; set(gca,'YColor','k'); ylabel('Bias LvNG'); end;
    title(expRefs{s,1}); caxis([-1 1]*max(abs(sites(s).Biases(:,1))));
    
    subplot(4,numSubjects,s+numSubjects);
    imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal'); set(imX,'alphadata',0.7);
    hold on;
    scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),dotSize,sites2(s).Bias(:,2),dotShape,'filled'); axis equal; colorbar;
    axis([-5 5 -6 4]); set(gca,'box','off','Xtick','','YTick',''); if s==1; set(gca,'YColor','k'); ylabel('Bias RvNG'); end;
    caxis([-1 1]*max(abs(sites(s).Biases(:,2))));
    
    subplot(4,numSubjects,s+2*numSubjects);
    imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal'); set(imX,'alphadata',0.7);
    hold on;
    scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),dotSize,sites2(s).Sens(:,1),dotShape,'filled'); axis equal; colorbar;
    axis([-5 5 -6 4]); set(gca,'box','off','Xtick','','YTick',''); if s==1; set(gca,'YColor','k'); ylabel('Sens CL (LvNG)'); end;
    caxis([-1 1]*max(abs(sites(s).CLeft(:,1))));
    
        subplot(4,numSubjects,s+3*numSubjects);
    imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal'); set(imX,'alphadata',0.7);
    hold on;
    scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),dotSize,sites2(s).Sens(:,2),dotShape,'filled'); axis equal; colorbar;
    axis([-5 5 -6 4]); set(gca,'box','off','Xtick','','YTick',''); if s==1; set(gca,'YColor','k'); ylabel('Sens CR (RvNG)'); end;
    caxis([-1 1]*max(abs(sites(s).CRight(:,2))));
    
    cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
        linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
    colormap(cmap);
    
    set(gcf,'Color','w');
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
clear PD;
for s = 1:numSubjects
    subplot(1,numSubjects+1,s);
    imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal');
    set(imX,'alphadata',0.7);
    hold on;
    
    PF = pivottable(l(s).data.laserIdx,[],l(s).data.feedbackType==1,'mean');
    scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),150,PF(2:end),'s','filled'); axis equal; colorbar;
    ylim([-6 4]);
    cax = PF(1) + [-1 1]*max(abs(PF(2:end)-PF(1)));
    PD(:,s)=PF(2:end)-PF(1);
%     cax(cax>1)=1;
    caxis(cax);
    hold off;
    title(expRefs{s,1})
    xlim([-5 5]);
    
    cmap = [ linspace(0,1,100)' ones(100,1) linspace(0,1,100)';
        linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
    colormap(flipud(cmap));
    
    set(gca,'xtick','','ytick','');
    set(gca,'box','off','XColor','w','YColor','w');
    
end
% 
% subplot(1,numSubjects+1,s+1);
% imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal');
% set(imX,'alphadata',0.7);
% hold on;
% scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),150,median(PD,2),'s','filled'); axis equal; colorbar;
% ylim([-6 4]);
% caxis([-1 1]*max(abs(median(PD,2))));
% hold off;
% title('Median')
% xlim([-5 5]);

cmap = [ linspace(0,1,100)' ones(100,1) linspace(0,1,100)';
    linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
colormap(flipud(cmap));

set(gca,'xtick','','ytick','');
set(gca,'box','off','XColor','w','YColor','w');

set(gcf,'Color','w');

%% LASER EFFECT MAP: pupil effects of the laser (no model)
figure;
kimg=imread('\\basket.cortexlab.net\home\stuff\kirkcaldie_brain_BW.PNG');
clear PD;
for s = 1:numSubjects
    subplot(1,numSubjects+1,s);
    imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal');
    set(imX,'alphadata',0.7);
    hold on;
    
    PF = pivottable(l(s).data.laserIdx,[],l(s).data.pupilenergy,'nanmean');
    scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),150,PF(2:end),'s','filled'); axis equal; colorbar;
    ylim([-6 4]);
    cax = PF(1) + [-1 1]*max(abs(PF(2:end)-PF(1)));
    PD(:,s)=PF(2:end)-PF(1);
%     cax(cax>1)=1;
    caxis(cax);
    hold off;
    title(expRefs{s,1})
    xlim([-5 5]);
    
    cmap = [ linspace(0,1,100)' ones(100,1) linspace(0,1,100)';
        linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
    colormap(flipud(cmap));
    
    set(gca,'xtick','','ytick','');
    set(gca,'box','off','XColor','w','YColor','w');
    
end
% 
% subplot(1,numSubjects+1,s+1);
% imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal');
% set(imX,'alphadata',0.7);
% hold on;
% scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),150,median(PD,2),'s','filled'); axis equal; colorbar;
% ylim([-6 4]);
% caxis([-1 1]*max(abs(median(PD,2))));
% hold off;
% title('Median')
% xlim([-5 5]);

cmap = [ linspace(0,1,100)' ones(100,1) linspace(0,1,100)';
    linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
colormap(flipud(cmap));

set(gca,'xtick','','ytick','');
set(gca,'box','off','XColor','w','YColor','w');

set(gcf,'Color','w');

%% LASER EFFECT MAP: reaction time effects (no model)
figure;
for s = 1:numSubjects
%     PF = pivottable(l(s).data.laserIdx,[],l(s).data.feedbackType==1,'mean');
%     PF = PF(2:end)-PF(1);
%     PF = [PF nan(length(PF),1)];
%     
    RTs = pivottable(l(s).data.laserIdx,l(s).data.response,l(s).data.RT,'median');
%     RTs = bsxfun(@minus,RTs(2:end,:),RTs(1,:));
    RTs = RTs(:,1:2);
    
%     Ls = pivottable(l(s).data.laserIdx,l(s).data.response,l(s).data.lickenergy,'median');
%     Ls = bsxfun(@minus,Ls(2:end,:),Ls(1,:));
%     Ls = Ls(:,1:2);
%     
    dotSize=150;
    dotShape='s';
    kimg=imread('D:\kirkcaldie_brain_BW.PNG');
    
    a=1;
    
    side={'left choices','right choices'};
    for lr=1:2
        subplot(2,numSubjects,s + (lr-1)*numSubjects);
        imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal');
        set(imX,'alphadata',0.7);
        hold on;
        scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),dotSize,RTs(2:end,lr),dotShape,'filled'); axis equal; colorbar;
        ylim([-6 4]);
        
        cax = RTs(1,lr) + [-1 1]*max(abs(RTs(2:end,lr)-RTs(1,lr)));
        caxis(cax);
        
        %         set(s(lr),'markeredgecolor',[1 1 1]*1,'linewidth',0);
        hold off;
        title(expRefs{s,1});
        xlim([-5 5]);
        a=a+1;
        set(gca,'xtick','','ytick','');
        set(gca,'box','off','XColor','w','YColor','w');
    end
    %
    cmap = [ linspace(0,1,100)' ones(100,1) linspace(0,1,100)';
        linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
    colormap(flipud(cmap));
end
set(gcf,'Color','w');

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

%% RELIABILITY: Fake inactivation location with two-stage fitting
%adds another laser site from a certain proportion of NoLaser trials. Since
%we know the laser has no true effect here, the estimated parameters from
%the fit can give an indication of the null distribution under this
%behavioural model.
Qshuffle = struct;
% proportion = 0.1;
NUM_SHUFFLES = 200;
figure;
for s = 1:numSubjects
    tab = tabulate(l(s).data.laserIdx);
    proportion = mean(tab(2:end,3))/100;
    
    numSessions = length(expRefs{s,2});
    nonL_idx = l(s).data.laserIdx==0;
    for shuff = 1:NUM_SHUFFLES
        disp([num2str(shuff) '/' num2str(NUM_SHUFFLES)]);
        nL = randsample(find(nonL_idx),round(proportion*sum(nonL_idx)));
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
               
%         %one stage fitting
%         fit=glmnet([X{s} QX],Y{s},'multinomial',glmnetSet(opts));
%         b=glmnetCoef(fit,0);
%         b=[b{1}-b{3} b{2}-b{3}];
%         b(1,:) = [];
        %
        %two stage fitting
        twoStageFitOpts = struct;
        twoStageFitOpts.intr = 1;

        %First fit non-laser portion of data
        XnL = X{s}(l(s).data.laserIdx==0 & QX(:,1)==0,:);
        YnL = Y{s}(l(s).data.laserIdx==0 & QX(:,1)==0);
%         fit=glmnet(XnL(:,1:(3*numSessions)),YnL,'multinomial',glmnetSet(twoStageFitOpts));
                fit=cvglmnet(XnL(:,1:(3*numSessions)),YnL,'multinomial',glmnetSet(twoStageFitOpts));

        %         % %Then fit laser portion of data, using nonLaser params as offsets
        XL = [X{s} QX];
        XL = XL(l(s).data.laserIdx>0 | QX(:,1)==1,:);
        YL = Y{s}(l(s).data.laserIdx>0 | QX(:,1)==1);
        twoStageFitOpts.intr = 0;
%         twoStageFitOpts.offset = glmnetPredict(fit,XL(:,1:(3*numSessions)),0,'link');
%         fit2=glmnet(XL(:,((3*numSessions)+1):end),YL,'multinomial',glmnetSet(twoStageFitOpts));
                        twoStageFitOpts.offset = cvglmnetPredict(fit,XL(:,1:(3*numSessions)),[],'link');
                fit2=cvglmnet(XL(:,((3*numSessions)+1):end),YL,'multinomial',glmnetSet(twoStageFitOpts));
        %
        %         %pull out parameter fits
        %         b=cvglmnetCoef(fit2,0);
        b=glmnetCoef(fit2,[]);
        b=[b{1}-b{3} b{2}-b{3}];
        b(1,:) = [];


        Qshuffle(s).Bias(shuff,:) = b(end-2,:);
        Qshuffle(s).CLeft(shuff,:) = b(end-1,:);
        Qshuffle(s).CRight(shuff,:) = b(end,:);
    end
    
    Qshuffle(s).Bias95CI = quantile(Qshuffle(s).Bias,[0.025 0.975]);
    Qshuffle(s).Bias99CI = quantile(Qshuffle(s).Bias,[0.005 0.995]);
    
    Qshuffle(s).CLeft95CI = quantile(Qshuffle(s).CLeft,[0.025 0.975]);
    Qshuffle(s).CLeft99CI = quantile(Qshuffle(s).CLeft,[0.005 0.995]);
    
    Qshuffle(s).CRight95CI = quantile(Qshuffle(s).CRight,[0.025 0.975]);
    Qshuffle(s).CRight99CI = quantile(Qshuffle(s).CRight,[0.005 0.995]);

    Qshuffle(s).CI95(:,:,1) = Qshuffle(s).Bias95CI;
    Qshuffle(s).CI95(:,:,2) = Qshuffle(s).CLeft95CI;
    Qshuffle(s).CI95(:,:,3) = Qshuffle(s).CRight95CI;
        
    Qshuffle(s).CI99(:,:,1) = Qshuffle(s).Bias99CI;
    Qshuffle(s).CI99(:,:,2) = Qshuffle(s).CLeft99CI;
    Qshuffle(s).CI99(:,:,3) = Qshuffle(s).CRight99CI;

    
    subplot(numSubjects,1,s);
    hist(Qshuffle(s).Bias);
    title(expRefs{s,1});
    
end

%% Check whether model is fitting well on non-laser trials %NEED TO ADD GLOBAL INTERCEPT
g = g.setModel('C50');
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
    set(ax,'Xticklabel','','box','off');
    set(ax(1),'Xticklabel',{'L','NG','R'},'xtick',1:3);
end

%% Check whether model is fitting well on laser trials %NEED TO ADD GLOBAL INTERCEPT
g = g.setModel('C50');
figure;
for s = 1:numSubjects
    numSessions = max(l(s).data.sessionID);
    actual = nan(size(l(s).inactivationSite,1)+1,4,3,numSessions);
    pred = nan(size(l(s).inactivationSite,1)+1,4,3,numSessions);
%     figure('name',expRefs{s,1});
    numSessions = max(l(s).data.sessionID);
    for session = 1:numSessions
        D = getrow(l(s).data,l(s).data.sessionID==session);
        tested_sites = unique(D.laserIdx);
        
        stim_labels = {'C=0','CL=CR','C=L','C=R'};
        stim = cell(1,4);
        stim{1} = sum(D.contrast_cond,2)==0;
        stim{2} = ~stim{1} & (D.contrast_cond(:,1)==D.contrast_cond(:,2));
        stim{3} = diff(D.contrast_cond,[],2)<0;
        stim{4} = diff(D.contrast_cond,[],2)>0;
        stim_idx = sum(bsxfun(@times,1:4,cell2mat(stim)),2);
        
        params = [sess(s).Biases(session,1),...
            sess(s).CLeft(session,1),...
            sess(s).CRight(session,1),...
            sess(s).Biases(session,2),...
            sess(s).CLeft(session,2),...
            sess(s).CRight(session,2),...
            c50_fits{s}(session,:)];
        
        for locidx = 1:length(tested_sites)
            loc = tested_sites(locidx);
            if loc==0
                paramsL = params;
            else
                paramsL = params + [sites(s).Biases(loc,1),...
                                    sites(s).CLeft(loc,1),...
                                    sites(s).CRight(loc,1),...
                                    sites(s).Biases(loc,2),...
                                    sites(s).CLeft(loc,2),...
                                    sites(s).CRight(loc,2),0,0];
            end
            
            for stim_type = 1:4
                pred(loc+1,stim_type,:,session) = mean(g.calculatePhat(paramsL, D.contrast_cond(stim{stim_type},:)));
               
                r = D.response(D.laserIdx==loc & stim_idx == stim_type);
                actual(loc+1,stim_type,:,session) = sum([r==1 r==2 r==3])/length(r);
                
            end
            
            
        end
        
%         %TODO: when combining 'actual' and 'pred' data across sessions there aren't an equal number of sites tested
%         actual(:,:,1,session) = pivottable(D.laserIdx,stim_idx,D.response==1,'mean');
%         actual(:,:,2,session) = pivottable(D.laserIdx,stim_idx,D.response==2,'mean');
%         actual(:,:,3,session) = pivottable(D.laserIdx,stim_idx,D.response==3,'mean');
        

    end
    
    dotSize = 100;
    coords = [4 3;l(s).inactivationSite(:,1) l(s).inactivationSite(:,2)];
    for stim_type = 1:4
        subplot(numSubjects,4,4*s - 4 + stim_type);
        scatter(coords(:,2),coords(:,1),dotSize,nanmean(actual(:,stim_type,1,:),4),'s','filled'); caxis([0 1]);
        hold on;
        scatter(coords(:,2)+10,coords(:,1),dotSize,nanmean(actual(:,stim_type,3,:),4),'s','filled');caxis([0 1]);
        scatter(coords(:,2)+20,coords(:,1),dotSize,nanmean(actual(:,stim_type,2,:),4),'s','filled');caxis([0 1]);
        
        scatter(coords(:,2),coords(:,1)-12,dotSize,nanmean(pred(:,stim_type,1,:),4),'s','filled'); caxis([0 1]);
        scatter(coords(:,2)+10,coords(:,1)-12,dotSize,nanmean(pred(:,stim_type,3,:),4),'s','filled'); caxis([0 1]);
        scatter(coords(:,2)+20,coords(:,1)-12,dotSize,nanmean(pred(:,stim_type,2,:),4),'s','filled'); caxis([0 1]);
        
        hold off
%         xlim([0 14]); 
        %             axis square;
        
        if s == 1
            title(stim_labels{stim_type});
        end
        
        if stim_type == 1
            ylabel(expRefs{s,1});
        end
    end
    
    ax=get(gcf,'children');
    set(ax,'Xticklabel','','Yticklabel','');
    set(ax(1),'XTickLabel',{'pL','pNG','pR'},'Xtick',[2 7 12],'YTickLabel',{'pred','actual'},'ytick',[-12 0])
%     set(ax(1),'Xticklabel',{'L','NG','R'},'xtick',1:3);
end

%% Bootstrapping of non-laser parameters to see if they tradeoff
%also plot model output to see whether parameter tradeoff causes variation
%in the model output
numiter = 100;
stim_labels = {'C=0','CL=CR','C=L','C=R'};
for s = 1:numSubjects
    D = getrow(l(s).data,l(s).data.laserIdx==0);
    XnL = X{s}(l(s).data.laserIdx==0,1:(3*numSessions));
    YnL = Y{s}(l(s).data.laserIdx==0);
    
    numSessions = max(D.sessionID);
   
    bootstrap = struct;
    
    model_output = nan(numiter,3,4,numSessions);
    for iter = 1:numiter
        disp(iter);
        
        Xb = XnL;
        Yb = YnL;
        for session = 1:numSessions
            sidx = D.sessionID==session;
            ridx = randsample(find(sidx),sum(sidx),true);
            Xb(sidx,:) = XnL(ridx,:);
            Yb(sidx) = YnL(ridx);
        end
%         imagesc(Xb);
        
        %fit
        fit=glmnet(Xb,Yb,'multinomial',glmnetSet(opts));
        b=glmnetCoef(fit,0);
        b=[b{1}-b{3} b{2}-b{3}];
        b(1,:) = [];
        
        bootstrap.Biases(:,:,iter) = b(1:3:end,:);
        bootstrap.CLeft(:,:,iter) = b(2:3:end,:);
        bootstrap.CRight(:,:,iter) = b(3:3:end,:);
        
        %model output
        D.y_hat = glmnetPredict(fit,XnL,0,'response');
        for session = 1:numSessions
            E = getrow(D,D.sessionID==session);
            stim = cell(1,4);
            stim{1} = sum(E.contrast_cond,2)==0;
            stim{2} = ~stim{1} & (E.contrast_cond(:,1)==E.contrast_cond(:,2));
            stim{3} = diff(E.contrast_cond,[],2)<0;
            stim{4} = diff(E.contrast_cond,[],2)>0;
            
            for stim_type = 1:4
                r = E.response(stim{stim_type});
                model_output(iter,:,stim_type,session)=mean(E.y_hat(stim{stim_type},:));
            end
        end
    end
    
    biases = permute(bootstrap.Biases,[3 2 1]);
    sensL = permute(bootstrap.CLeft,[3 2 1]);
    sensR = permute(bootstrap.CRight,[3 2 1]);

        %reparameterise the parameters at lambda=0
    splus = (sensL + sensR)/2;
    sminus = (sensL - sensR)/2;
    czero = -biases./splus;
    
    labels_old = {'b_{LvNG}','sL_{LvNG}','sR_{LvNG}','sR_{RvNG}','sL_{RvNG}','b_{RvNG}'};
    labels_new = {'czero_{LvNG}','splus_{LvNG}','sminus_{LvNG}','sminus_{RvNG}','splus_{RvNG}','czero_{RvNG}'};
    
    figure('name',expRefs{s,1});
    for session = 1:numSessions
        old_params = [biases(:,1,session) sensL(:,1,session) sensR(:,1,session) sensR(:,2,session) sensL(:,2,session) biases(:,2,session)];
        new_params = [czero(:,1,session) splus(:,1,session) sminus(:,1,session) sminus(:,2,session) splus(:,2,session) czero(:,2,session)];
        
        subplot(2,numSessions,session);
        imagesc(corrcoef(old_params));
        if session == 1
            set(gca,'XTickLabel',labels_old,'Xtick',1:6,'YTickLabel',labels_old,'Ytick',1:6,'XTickLabelRotation',90,'fontsize',9);
            title('old parameterisation');
        else
            set(gca,'XTickLabel','','YTickLabel','');
        end
        axis square; caxis([-1 1]);
        subplot(2,numSessions,session + numSessions);
        imagesc(corrcoef(new_params));
       if session == 1
            set(gca,'XTickLabel',labels_new,'Xtick',1:6,'YTickLabel',labels_new,'Ytick',1:6,'XTickLabelRotation',90,'fontsize',9);
            title('new parameterisation');
        else
            set(gca,'XTickLabel','','YTickLabel','');
        end
        axis square;caxis([-1 1]);
        
cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
        linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
    colormap(cmap);
%         figure('name',[expRefs{s,1} ' session ' num2str(session) ' CONVENTIONAL PARAMETERS']);
%         gplotmatrix([biases(:,1,session) sensL(:,1,session) sensR(:,1,session) sensR(:,2,session) sensL(:,2,session) biases(:,2,session)],[],[],'b','.',5,'on','none',labels_old);
%         
%         figure('name',[expRefs{s,1} ' session ' num2str(session) ' NEW PARAMETERS']);
%         gplotmatrix([czero(:,1,session) splus(:,1,session) sminus(:,1,session) sminus(:,2,session) splus(:,2,session) czero(:,2,session)],[],[],'b','.',5,'on','none',labels_new);
%         drawnow;
    end
    
    f=figure('name',expRefs{s,1});
    for session = 1:numSessions 
        for stim_type = 1:4
            subplot(numSessions,4,4*session - 4 + stim_type);
            plot(1:3,model_output(:,:,stim_type,session),'b.'); ylim([0 1]); xlim([0.5 3.5]);
%             boxplot(model_output(:,[1 3 2],stim_type,session),'plotstyle','compact'); ylim([0 1]);
            if session == 1
                title(stim_labels{stim_type});
            end
        end
    end
    ax=get(f,'children');
    set(ax,'Xticklabel','','Yticklabel','','box','off');
    set(ax(1),'XTickLabel',{'pL','pR','pNG'});
    ylabel(ax(1),'probability');
end

%% LASER EFFECT MAP: psychometric curves) %NEED TO ADD GLOBAL INTERCEPT
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
    numSessions = length(expRefs{s,2});
    cv = cvpartition(size(X{s},1),'kfold',5);
    
    p_hats = nan(cv.NumObservations,3);
    for iter = 1:cv.NumTestSets
        disp(iter);
        trainX = X{s}(cv.training(iter),:);
        trainY = Y{s}(cv.training(iter));
        
        testX = X{s}(cv.test(iter),:);
        testY = Y{s}(cv.test(iter));
        
        %one stage fitting without a global intercept
        opts.intr=0;
        fit=cvglmnet(trainX,trainY,'multinomial',glmnetSet(opts));
        p = cvglmnetPredict(fit,testX,'lambda_min','response');
        p_hats(cv.test(iter),1) = p(sub2ind(size(p), [1:length(testY)]', testY));
        
        %one stage fitting with a global intercept
        opts.intr=1;
        fit = cvglmnet(trainX,trainY,'multinomial',glmnetSet(opts));
        p = cvglmnetPredict(fit,testX,'lambda_min','response');
        p_hats(cv.test(iter),2) = p(sub2ind(size(p), [1:length(testY)]', testY));
        
        %two stage fitting with a nonLaser global intercept only
        twoStageFitOpts = struct;
        twoStageFitOpts.intr = 1;
        nonL_idx = l(s).data.laserIdx(cv.training(iter))==0;
        fit1=cvglmnet(trainX(nonL_idx,1:(3*numSessions)),trainY(nonL_idx),'multinomial',glmnetSet(twoStageFitOpts));
        
        L_idx = l(s).data.laserIdx(cv.training(iter))>0;
        twoStageFitOpts.offset = cvglmnetPredict(fit1,trainX(L_idx,1:(3*numSessions)),[],'link');
        twoStageFitOpts.intr = 0;
        fit2=cvglmnet(trainX(L_idx,((3*numSessions)+1):end),trainY(L_idx),'multinomial',glmnetSet(twoStageFitOpts));
        
        nL_idx_test = l(s).data.laserIdx(cv.test(iter))==0;
        L_idx_test = l(s).data.laserIdx(cv.test(iter))>0;
        
        offset = cvglmnetPredict(fit1,testX(:,1:(3*numSessions)),[],'link');
        p = cvglmnetPredict(fit2,testX(:,((3*numSessions)+1):end),[],'response',[],offset);
        p_hats(cv.test(iter),3) = p(sub2ind(size(p), [1:length(testY)]', testY));
       
    end
    
    tab=tabulate(Y{s}); tab=tab(:,3)/100;
    guess_bpt = sum(tab.*log2(tab));
    
    cv_gof(s,:) = mean(log2(p_hats))-guess_bpt;
end 
disp('done');

bar(cv_gof); ylabel('loglik [bits] relative to guessing @~-1.5');
set(gca,'XTickLabel',expRefs(:,1),'xtick',1:numSubjects);

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

%% MISC: Plot a table of contrast stimulus settings
p = load(dat.expFilePath(expRefs{end,2}{end},'parameters','m'));
p = p.parameters;
counts = p.numRepeats(p.rewardOnStimulus(2,:)==1.1);
c = p.visCueContrast(:,p.rewardOnStimulus(2,:)==1.1);

uC = unique(c(:));
for cl = 1:length(uC)
    for cr = 1:length(uC)
        T(cl,cr)=sum(counts(c(1,:) == uC(cl) & c(2,:) == uC(cr)));
    end
end

imagesc(uC,uC,T); set(gca,'ydir','normal','box','off'); xlabel('CR'); ylabel('CL'); 
set(gcf,'Color','w');

%% (NO GLM) Test IIA assumption in data
% figure;
for s = 1:numSubjects
    numSessions = max(l(s).data.sessionID);
    
    figure('name',expRefs{s,1});
    pratios=[];
    %     for se = 1:numSessions
    uC = unique(l(s).data.contrast_cond(:));
    %         D = getrow(l(s).data,l(s).data.sessionID==se & l(s).data.laserIdx==0);
    D = getrow(l(s).data,l(s).data.laserIdx==0);
    
    for relC = 1:length(uC)
        %             ax(se)=subplot(numSessions,length(uC),length(uC)*(se-1) + relC );
        ax(relC)=subplot(length(uC),1,relC);
        
        for c=1:length(uC)
            %test pL/pNG over different right contrasts
            r=D.response(D.contrast_cond(:,2)==uC(c) & D.contrast_cond(:,1)==uC(relC));
            pratios(c,1,se) = sum(r==1)/sum(r==3);
            
            %test pR/pNG over different left contrasts
            r=D.response(D.contrast_cond(:,1)==uC(c) & D.contrast_cond(:,2)==uC(relC));
            pratios(c,2,se) = sum(r==2)/sum(r==3);
        end
        plot(uC,log(pratios(:,:,se)));  xlim([0 .6]);
        set(gca,'box','off');
        
        if relC == length(uC)
            xlabel('Contrast of irrelevant alternative');
            ylabel('log odds'); legend('ln[pL/pNG]','ln[pR/pNG]');
        end
        
        if relC ~= length(uC)
            set(gca,'XTickLabel','','XColor','w');
        end
        
        title(['Relevant contrast: ' num2str(uC(relC))]);
    end
    
    %         if se == 1
    %             title(expRefs{s,1});
    %         end
    %     end
    set(gcf,'Color','w');
    %     set(ax,'XColor','w');
    %     linkaxes(ax,'x');
end

%% Crossval
figure;
for s = 1:numSubjects
   numSessions =  max(l(s).data.sessionID);
   subplot(1,numSubjects,s);
   
   ll = [];
   for se = 1:numSessions
       D = getrow(l(s).data,l(s).data.laserIdx==12 & l(s).data.sessionID==se);
       g = GLM(D);
       
       g.regularise = @(b)0.001*sum(b.^2);
       
       ll(se,1) = mean(log2(g.setModel('C50-subset').fitCV(5).p_hat)) - g.guess_bpt;
       ll(se,2) = mean(log2(g.setModel('C50').fitCV(5).p_hat)) - g.guess_bpt;
   end
   
   bar(ll','stacked'); set(gca,'box','off','XTickLabel',{'C50-subset','C50'},'XTick',[1 2]);
 
   title(expRefs{s,1});
   
   if s==1
       ylabel('log_2 likelihood relative to guessing [sum over sessions]');
   end
   
   set(gcf,'Color','w'); colormap gray;
end