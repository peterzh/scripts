%% 
expt = 'sparse_unilateral_2D';
bilateral_mirror = 0;
MOREDATA=0;


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


numSubjects = size(expRefs,1);
c50_fits = cell(1,size(expRefs,1));
cn_fits = cell(1,size(expRefs,1));

for s = 1:numSubjects
    name = expRefs{s,1};
    
    %Combine all sessions
    D=struct;
    %i=1;
    
    figure('name',expRefs{s,1});
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
        
                
        cVal = unique(d.contrast_cond(:));
        prop=nan(length(cVal),length(cVal),3);
        
        for cl = 1:length(cVal)
            for cr = 1:length(cVal)
                r = d.response(d.contrast_cond(:,1) == cVal(cl) & d.contrast_cond(:,2) == cVal(cr) & d.laserIdx==0);
                for i=1:3
                    prop(cl,cr,i) = sum(r==i)/length(r);
                end
            end
        end
        
        titles = {'p( Left | c)','p( Right | c)','p( NoGo | c)'};
        for i=1:3
            subplot(length(expRefs{s,2}),3,3*session - 3 + i);
            
            imagesc(cVal,cVal,prop(:,:,i),[0 1]);
            set(gca,'YDir','normal','box','off');
            
            %                         xlabel('Contrast right');
            %                         ylabel('Contrast left');
%             title(titles{i});
            axis square;
            set(gca,'XTick','','YTick','','XColor','w','YColor','w');
            if i > 1
                set(gca,'XTick','','ytick','');
            end
            
%             if i == 1
%                 %                                 xlabel('Contrast right');
%                 ylabel('Contrast left');
%             end
        end
    end
    
    set(gcf,'color','w');
    % figure; subplot(1,2,1); hist(c50sub.n); xlabel('n');
    % subplot(1,2,2); hist(c50sub.c50); xlabel('c50');
    D = rmfield(D,'laserIdx');
    D = getrow(D,D.repeatNum==1);
    l(s) = laserGLM(D);
    % save('C:\Users\Peter\Desktop\scripts\LaserOmnibus_data.mat', 'l');
end

% Fix all C50 and N, and ^N parameters constant for each mouse
c50_fits=cellfun(@(a)(repmat(mean(a,1),size(a,1),1)),c50_fits,'uni',0);
cn_fits=cellfun(@(a)(repmat(mean(a,1),size(a,1),1)),cn_fits,'uni',0);

%% Define and fit model
sites = struct('Biases',{},'CLeft',{},'CRight',{},'laserSite',{});
sess = struct('Biases',{},'CLeft',{},'CRight',{});

% sites2 = struct('Bias',{},'Sens',{},'laserSite',{});
% sess2 = struct('Bias',{},'Sens',{});

X = cell(1,numSubjects); X2=cell(1,numSubjects);
Y = cell(1,numSubjects); Y2=cell(1,numSubjects);
figure;

fitMode = 1;
for s = 1:numSubjects
    numSessions = max(l(s).data.sessionID);
    numSites = max(l(s).data.laserIdx);
    numTrials = size(l(s).data.response,1);
    
    Y{s} = l(s).data.response;
    Y2{s} = Y{s}; Y2{s}(Y2{s}==3)=0; Y2{s}=Y2{s}+1;
    
    %construct design matrix X (for loop for ease of understanding not speed!)
    X{s} = sparse(numTrials,3*numSessions + 3*numSites);
%     X2{s} = sparse(numTrials,3*numSessions + numSites);
    X2{s} = sparse(numTrials,3 + numSites);
%     X_LvNG{s} = sparse(numTrials,2*numSessions + 2*numSites);
%     X_RvNG{s} = sparse(numTrials,2*numSessions + 2*numSites);

    contrast_rep_fcn = @(c,n,c50)((c.^n)/(c.^n + c50^n));
    for i = 1:numTrials
        
        thisSession = l(s).data.sessionID(i);
        cfn = @(c)contrast_rep_fcn(c,c50_fits{s}(thisSession,1),c50_fits{s}(thisSession,2));
        X{s}(i,3*thisSession - 2) = 1;
        X{s}(i,3*thisSession - 1) = cfn(l(s).data.contrast_cond(i,1));
        X{s}(i,3*thisSession - 0) = cfn(l(s).data.contrast_cond(i,2));
        
%         X2{s}(i,3*thisSession - 2) = 1;
%         X2{s}(i,3*thisSession - 1) = cfn(l(s).data.contrast_cond(i,1));
%         X2{s}(i,3*thisSession - 0) = cfn(l(s).data.contrast_cond(i,2));
%         
%         if thisSession == 1
%             X2{s}(i,3*thisSession - 2) = -1;
%         end
        
        X2{s}(i,1) = 1;
        X2{s}(i,2) = cfn(l(s).data.contrast_cond(i,1));
        X2{s}(i,3) = cfn(l(s).data.contrast_cond(i,2));

        thisSite = l(s).data.laserIdx(i);
        if thisSite ~= 0
            X{s}(i,3*numSessions + 3*thisSite - 2) = 1;
            X{s}(i,3*numSessions + 3*thisSite - 1) = cfn(l(s).data.contrast_cond(i,1));
            X{s}(i,3*numSessions + 3*thisSite - 0) = cfn(l(s).data.contrast_cond(i,2));
            
%             X2{s}(i,3*numSessions + thisSite) = 1;
            X2{s}(i,3 + thisSite) = 1;
%             X2{s}(i,3 + 3*thisSite - 1) = cfn(l(s).data.contrast_cond(i,1));
%             X2{s}(i,3 + 3*thisSite - 0) = cfn(l(s).data.contrast_cond(i,2));

        end
    end
       
    opts=struct;
    opts.intr=1; %global intercept
    opts.alpha=1; %lasso L1 penalty
    
    fit=cvglmnet(X{s},Y{s},'multinomial',glmnetSet(opts));
    b=cvglmnetCoef(fit,[]);
    disp([expRefs{s,1} ' cv lambda 1se: ' num2str(fit.lambda_min)]);
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

%% Fit hierarchical logit model instead using mnrfit (no regularisation and reduced model)
% Overparameterisation is a bigger problem when there is no regularisation
% and so I reduce the model here to fit a hierarchical model
mode = 2;
figure('name',['Mode: ' num2str(mode)]);
laserEffect = [];
for s = 1:numSubjects
    disp(s);
    numSessions = max(l(s).data.sessionID);
    numSites = size(l(s).inactivationSite,1);
    
    switch(mode)
        case 1 %Combine all sessions + all laser effects as biasing
            xh = X2{s};
            [b,dev,stat]=mnrfit(xh(:,2:end),Y2{s},'model','hierarchical');
            laserEffect(:,:,s) = b(end-numSites+1:end,:);
            
            maxEffect = max(max(abs(laserEffect(:,:,s))));
            sig = 50*(0.1 + 2*(stat.p(end-numSites+1:end,:)<0.05) + 2*(stat.p(end-numSites+1:end,:)<0.01));
%             sig = 50*abs(stat.t(end-numSites+1:end,:));
            
        case 2 %Have seperate non-laser session fits + all laser effects as biasing
            xh = [X{s}(:,1:3*numSessions) X{s}(:,3*numSessions+1 :3: end)];
            [b,dev,stat]=mnrfit(xh(:,2:end),Y2{s},'model','hierarchical');
            laserEffect(:,:,s) = b(end-numSites+1:end,:);
            
            maxEffect = max(max(abs(laserEffect(:,:,s))));
            sig = 50*(0.1 + 2*(stat.p(end-numSites+1:end,:)<0.05) + 2*(stat.p(end-numSites+1:end,:)<0.01));
%             sig = 50*abs(stat.t(end-numSites+1:end,:));
            
        case 3 %Separate non-laser session fits + laser biasing + laser sensitivity effects
            xh = X{s};
            [b,dev,stat]=mnrfit(xh(:,2:end),Y2{s},'model','hierarchical');
            laserEffect(:,:,s,1) = b(3*numSessions + 1:3:end,:);
%             laserEffect(:,:,s,2) = b(3*numSessions + 2:3:end,:);
%             laserEffect(:,:,s,3) = b(3*numSessions + 3:3:end,:);
            sig = 50*(0.1 + 2*(stat.p(3*numSessions + 1:3:end,:)<0.05) + 2*(stat.p(3*numSessions + 1:3:end,:)<0.01));
%             sig = 50*abs(stat.t(3*numSessions + 1:3:end,:));


        otherwise
            error('Not implemented');
    end
    
    ap = l(s).inactivationSite(:,1);
    ml = l(s).inactivationSite(:,2);
    out = laserEffect;
    if bilateral_mirror == 1
        ap = [ap;ap];
        ml = [ml;-ml];
        out = [out;out];
        sig=[sig;sig];
    end
    
    
    subplot(numSubjects,3,3*s-2);
    scatter(ml,ap,sig(:,1),out(:,1),'s','filled'); axis equal;
    caxis([-1 1]*maxEffect); ylabel(expRefs{s,1}); if s == 1; title('Laser biasing pNG/pG'); end;
    
    subplot(numSubjects,3,3*s-1);
    scatter(ml,ap,sig(:,2),out(:,2),'s','filled'); axis equal;
    caxis([-1 1]*maxEffect);  if s == 1; title('Laser biasing pL/pR'); end;
    
    subplot(numSubjects,3,3*s); 
    imagesc(stat.coeffcorr); caxis([-1 1]); axis square; if s == 1; title('Correlation among params'); end;

end
cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
        linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
    colormap(cmap);
    set(get(gcf,'children'),'box','off','xtick','','ytick','','xcolor','w','ycolor','k');

figure;
avg = mean(laserEffect,3);
ap = l(s).inactivationSite(:,1);
ml = l(s).inactivationSite(:,2);
out = avg;
if bilateral_mirror == 1
    ap = [ap;ap];
    ml = [ml;-ml];
    out = [out;out];
end
subplot(1,2,1);
scatter(ml,ap,200,out(:,1),'s','filled');  caxis([-1 1]*5); axis equal;
subplot(1,2,2);
scatter(ml,ap,200,out(:,2),'s','filled');  caxis([-1 1]*5); axis equal;
    colormap(cmap); set(get(gcf,'children'),'box','off','xtick','','ytick','','xcolor','w','ycolor','w');
    set(gcf,'color','w'); colorbar;



%% Two-stage fitting alternative
%Fits nonLaser parameters first and then the laser effects
fitMode = 2;
lambda1se = nan(numSubjects,2);
for s = 1:numSubjects
    twoStageFitOpts = struct;
    twoStageFitOpts.intr = 1;
    twoStageFitOpts.alpha = 0; %I want ridge
    
    numSessions = length(expRefs{s,2});
    
    %First fit non-laser portion of data
    nonL_idx = l(s).data.laserIdx==0;
    XnL = X{s}(nonL_idx,:);
    YnL = Y{s}(nonL_idx);
    fit1=cvglmnet(XnL(:,1:(3*numSessions)),YnL,'multinomial',glmnetSet(twoStageFitOpts));
    disp([expRefs{s,1} ' nonLaser fit lambda 1se: ' num2str(fit1.lambda_1se)]);
    
    %pull out parameter fits
    b=cvglmnetCoef(fit1,[]);
    b=[b{1}-b{3} b{2}-b{3}];
    sess(s).GlobalIntercept = b(1,:);
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
    disp([expRefs{s,1} ' Laser fit lambda 1se: ' num2str(fit2.lambda_1se)]);

    %pull out parameter fits
    b=cvglmnetCoef(fit2,[]);
    b=[b{1}-b{3} b{2}-b{3}];
    b(1,:) = [];
    sites(s).Biases = b(1:3:end,:);
    sites(s).CLeft = b(2:3:end,:);
    sites(s).CRight = b(3:3:end,:);    
    
    lambda1se(s,:) = [fit1.lambda_1se fit2.lambda_1se];
end

%% SIGNIFICANCE: Fake inactivation location with two-stage fitting
%adds another laser site from a certain proportion of NoLaser trials. Since
%we know the laser has no true effect here, the estimated parameters from
%the fit can give an indication of the null distribution under this
%behavioural model.
Qshuffle = struct;
% proportion = 0.1;
NUM_SHUFFLES = 250;
figure;
for s = 1:numSubjects
    
    FIRST_LAMBDA = lambda1se(s,1);
    SEC_LAMBDA = lambda1se(s,2);

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
        twoStageFitOpts.alpha = 0; %I want ridge
        
        %First fit non-laser portion of data
        XnL = X{s}(l(s).data.laserIdx==0 & QX(:,1)==0,:);
        YnL = Y{s}(l(s).data.laserIdx==0 & QX(:,1)==0);
        %         fit=glmnet(XnL(:,1:(3*numSessions)),YnL,'multinomial',glmnetSet(twoStageFitOpts));
        fit1=glmnet(XnL(:,1:(3*numSessions)),YnL,'multinomial',glmnetSet(twoStageFitOpts));
        
        %         % %Then fit laser portion of data, using nonLaser params as offsets
        XL = [X{s} QX];
        XL = XL(l(s).data.laserIdx>0 | QX(:,1)==1,:);
        YL = Y{s}(l(s).data.laserIdx>0 | QX(:,1)==1);
        twoStageFitOpts.intr = 0;
        %         twoStageFitOpts.offset = glmnetPredict(fit,XL(:,1:(3*numSessions)),0,'link');
        %         fit2=glmnet(XL(:,((3*numSessions)+1):end),YL,'multinomial',glmnetSet(twoStageFitOpts));
        twoStageFitOpts.offset = glmnetPredict(fit1,XL(:,1:(3*numSessions)),FIRST_LAMBDA,'link');
        fit2=glmnet(XL(:,((3*numSessions)+1):end),YL,'multinomial',glmnetSet(twoStageFitOpts));
        %
        %         %pull out parameter fits
        %         b=cvglmnetCoef(fit2,0);
        b=glmnetCoef(fit2,SEC_LAMBDA);
        b=[b{1}-b{3} b{2}-b{3}];
        b(1,:) = [];

        Qshuffle(s).Bias(shuff,:) = b(end-2,:);
        Qshuffle(s).CLeft(shuff,:) = b(end-1,:);
        Qshuffle(s).CRight(shuff,:) = b(end,:);
% 
%         %one-stage fitting LvNG and RvNG separately DOESNT WORK
%         Xn = [X{s} QX];
%         Yn = Y{s};
%         
%         Yn_LvNG = Yn(Yn~=2); Yn_LvNG(Yn_LvNG==3)=0;
%         Yn_RvNG = Yn(Yn~=1); Yn_RvNG(Yn_RvNG==3)=0; Yn_RvNG(Yn_RvNG==2)=1;
%         Xn_LvNG = Xn; Xn_LvNG(Yn==2,:) = []; Xn_LvNG(:,3:3:end) = [];
%         Xn_RvNG = Xn; Xn_RvNG(Yn==1,:) = []; Xn_RvNG(:,2:3:end) = [];
% 
%         b2=[];
%         fit=glmnet(Xn_LvNG, Yn_LvNG, 'binomial',glmnetSet(opts));
%         b2(:,1)=glmnetCoef(fit,0.001);
%         fit=glmnet(Xn_RvNG, Yn_RvNG, 'binomial',glmnetSet(opts));
%         b2(:,2)=glmnetCoef(fit,0.001);
%         b2(1,:) = [];
% 
%         Qshuffle(s).Bias(shuff,:) = b2(end-1,:);
%         Qshuffle(s).Sens(shuff,:) = b2(end,:);
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
    
    
    
    subplot(numSubjects,3,3*s-2);
    hist(Qshuffle(s).Bias); title(expRefs{s,1}); ylabel('Bias');
    
    subplot(numSubjects,3,3*s-1);
    hist(Qshuffle(s).CLeft); ylabel('CLeft');
    
    subplot(numSubjects,3,3*s);
    hist(Qshuffle(s).CRight); ylabel('CRight');
    
end

% parameter correlations in the fake site
c = [];
figure;
for s = 1:numSubjects
    ax(1)=subplot(4,numSubjects,s);
    hist(Qshuffle(s).Bias); title(expRefs{s,1}); if s == 1; xlabel('Bias'); end
    
        ax(2)=subplot(4,numSubjects,numSubjects+s);
    hist(Qshuffle(s).CLeft);  if s == 1; xlabel('CLeft'); end
    
        ax(3)=subplot(4,numSubjects,2*numSubjects+s);
    hist(Qshuffle(s).CRight); if s == 1; xlabel('CRight'); end
    
    set(ax,'ytick','','ycolor','w');

    subplot(4,numSubjects,3*numSubjects+s);
    
    c(:,:,s) = corrcoef([Qshuffle(s).Bias(:,1) Qshuffle(s).CLeft(:,1) Qshuffle(s).CRight(:,1) Qshuffle(s).CRight(:,2) Qshuffle(s).CLeft(:,2) Qshuffle(s).Bias(:,2)]);
    imagesc(c(:,:,s)); caxis([-1 1]); colormap(cmap); axis square; 
    
        
    labels = {'b_{LvNG}','sL_{LvNG}','sL_{RvNG}','sR_{RvNG}','sL_{RvNG}','b_{RvNG}'};
    if s == 1
        set(gca,'XTickLabel',labels,'Xtick',1:6,'YTickLabel',labels,'Ytick',1:6,'XTickLabelRotation',90,'fontsize',9);
    else
        set(gca,'Xtick','','Ytick','');
    end
end
set(gcf,'color','w');
set(get(gcf,'children'),'box','off');

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
    dotShape='s';
    kimg=imread('D:\kirkcaldie_brain_BW.PNG');
    
    a=1;
    figure('name',expRefs{s,1},'units','normalized','position',[.1 .1 .7 .7]);
    
    side={'log(\pi_L/\pi_{NG})','log(\pi_R/\pi_{NG})'};
    for d = 1:length(toDisplay)
        for lr=1:2
            subplot(length(toDisplay),2,a);
%             imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal');
%             set(imX,'alphadata',0.7);
%             hold on;

%             
            if exist('Qshuffle','var') && fitMode == 2 %if shuffle analysis done using 2 stage fitting
                ci_95 = Qshuffle(s).CI95(:,lr,d);
                ci_99 = Qshuffle(s).CI99(:,lr,d);
                
                sig95 = double(toDisplay{d}(:,lr) < ci_95(1) | ci_95(2) < toDisplay{d}(:,lr));
                sig99 = double(toDisplay{d}(:,lr) < ci_99(1) | ci_99(2) < toDisplay{d}(:,lr));
                
                sig = sum([sig95 sig99],2);
                sig(sig==0)=0.2; 
                
                ap = l(s).inactivationSite(:,1);
                ml = l(s).inactivationSite(:,2);
                out = toDisplay{d}(:,lr);
                if bilateral_mirror == 1
                    ap = [ap;ap];
                    ml = [ml;-ml];
                    out = [out;out];
                    sig=[sig;sig];
                end
                

                s1=scatter(ml,ap,sig*80,out,dotShape,'filled'); axis equal; colorbar;
                s1.MarkerEdgeColor=[1 1 1]*0.9;
        
        
%                 %                 scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),sig*15,'ok','filled');
%                 s1=scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),sig99*150,toDisplay{d}(:,lr),dotShape,'filled'); axis equal; colorbar;
%                 hold on; %s1.MarkerEdgeColor=[1 1 1]*0.8;
%                 s1=scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),sig95*150*0.6,toDisplay{d}(:,lr),dotShape,'filled'); axis equal; colorbar;
%                 hold off; %s1.MarkerEdgeColor=[1 1 1]*0.8;
            else
                
                ap = l(s).inactivationSite(:,1);
                ml = l(s).inactivationSite(:,2);
                out = toDisplay{d}(:,lr);
                if bilateral_mirror == 1
                    ap = [ap;ap];
                    ml = [ml;-ml];
                    out = [out;out];
                end
                scatter(ml,ap,dotSize,out,dotShape,'filled'); axis equal; colorbar;
                
            end
            
            
%             s1 = scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),dotSize,toDisplay{d}(:,lr),dotShape,'filled'); axis equal; colorbar;
            ylim([-6 4]);
            
            %             pcntl = quantile(toDisplay{d}(:),[0.05 0.95]);
            %             caxis([-1 1]*max(abs(pcntl)));
            if d == 1
                caxis([-1 1]);
            else
                caxis([-1 1]*2);
            end
            
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

% %% LASER EFFECT MAP: LvsNG and RvsNG maps CONDENSED OVER SUBJECTS 
% figure;
% for s = 1:numSubjects
%     labels = {'bias','CL sensitivity','CR sensitivity'};
% 
%     dotSize=150;
%     dotShape='o';
%     kimg=imread('D:\kirkcaldie_brain_BW.PNG');
%     
%     subplot(4,numSubjects,s);
%     imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal'); set(imX,'alphadata',0.7);
%     hold on;
%     scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),dotSize,sites2(s).Bias(:,1),dotShape,'filled'); axis equal; colorbar;
%     axis([-5 5 -6 4]); set(gca,'box','off','Xtick','','YTick',''); if s==1; set(gca,'YColor','k'); ylabel('Bias LvNG'); end;
%     title([expRefs{s,1} ' n=' num2str(size(X{s},1))]); caxis([-1 1]*max(abs(sites2(s).Bias(:,1))));
%     
%     subplot(4,numSubjects,s+numSubjects);
%     imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal'); set(imX,'alphadata',0.7);
%     hold on;
%     scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),dotSize,sites2(s).Bias(:,2),dotShape,'filled'); axis equal; colorbar;
%     axis([-5 5 -6 4]); set(gca,'box','off','Xtick','','YTick',''); if s==1; set(gca,'YColor','k'); ylabel('Bias RvNG'); end;
%     caxis([-1 1]*max(abs(sites2(s).Bias(:,2))));
%     
%     subplot(4,numSubjects,s+2*numSubjects);
%     imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal'); set(imX,'alphadata',0.7);
%     hold on;
%     scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),dotSize,sites2(s).Sens(:,1),dotShape,'filled'); axis equal; colorbar;
%     axis([-5 5 -6 4]); set(gca,'box','off','Xtick','','YTick',''); if s==1; set(gca,'YColor','k'); ylabel('Sens CL (LvNG)'); end;
%     caxis([-1 1]*max(abs(sites2(s).Sens(:,1))));
%     
%         subplot(4,numSubjects,s+3*numSubjects);
%     imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal'); set(imX,'alphadata',0.7);
%     hold on;
%     scatter(l(s).inactivationSite(:,2),l(s).inactivationSite(:,1),dotSize,sites2(s).Sens(:,2),dotShape,'filled'); axis equal; colorbar;
%     axis([-5 5 -6 4]); set(gca,'box','off','Xtick','','YTick',''); if s==1; set(gca,'YColor','k'); ylabel('Sens CR (RvNG)'); end;
%     caxis([-1 1]*max(abs(sites2(s).Sens(:,2))));
%     
%     cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
%         linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
%     colormap(cmap);
%     
% end
%     set(gcf,'Color','w');
% 
% set(get(gcf,'children'),'box','off','fontsize',20)


%% LASER EFFECT MAP: LvsNG and RvsNG maps MEAN OVER SUBJECTS 
figure; %average map

BI = []; CLS = []; CRS = [];
for s = 1:numSubjects
    BI(:,:,s) = sites(s).Biases;
    CLS(:,:,s) = sites(s).CLeft;
    CRS(:,:,s) = sites(s).CRight;
end
BI = mean(BI,3);
CLS = mean(CLS,3);
CRS = mean(CRS,3);

dotSize=150;
dotShape='s';
kimg=imread('D:\kirkcaldie_brain_BW_outline.PNG');

ap = l(s).inactivationSite(:,1);
ml = l(s).inactivationSite(:,2);
out = [BI CLS CRS];
if bilateral_mirror == 1
    ap = [ap;ap];
    ml = [ml;-ml];
    out = [out;out];
end

subplot(3,2,1);
% imX=image(-5:1:5,4:-1:-6,kimg); 
% axis square; set(gca,'ydir','normal'); set(imX,'alphadata',0.7);
% hold on;
s1=scatter(ml,ap,dotSize,out(:,1),dotShape,'filled'); axis equal; colorbar;
axis([-5 5 -6 4]); set(gca,'box','off','Xtick','','YTick','','XColor','w','YColor','w'); title('Bias_L'); 
caxis([-1 1]*1);

subplot(3,2,2);
% imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal'); set(imX,'alphadata',0.7);
% hold on;
s1=scatter(ml,ap,dotSize,out(:,2),dotShape,'filled'); axis equal; colorbar;
axis([-5 5 -6 4]); set(gca,'box','off','Xtick','','YTick','','XColor','w','YColor','w'); title('Bias_R');
caxis([-1 1]*1);

subplot(3,2,3);
% imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal'); set(imX,'alphadata',0.7);
% hold on;
s1=scatter(ml,ap,dotSize,out(:,3),dotShape,'filled'); axis equal; colorbar;
axis([-5 5 -6 4]); set(gca,'box','off','Xtick','','YTick','','XColor','w','YColor','w'); title('Sens_L (LvNG)');
caxis([-1 1]*2);

subplot(3,2,4);
% imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal'); set(imX,'alphadata',0.7);
% hold on;
s1=scatter(ml,ap,dotSize,out(:,4),dotShape,'filled'); axis equal; colorbar;
axis([-5 5 -6 4]); set(gca,'box','off','Xtick','','YTick','','XColor','w','YColor','w'); title('Sens_L (RvNG)');
caxis([-1 1]*2);

subplot(3,2,5);
% imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal'); set(imX,'alphadata',0.7);
% hold on;
s1=scatter(ml,ap,dotSize,out(:,5),dotShape,'filled'); axis equal; colorbar;
axis([-5 5 -6 4]); set(gca,'box','off','Xtick','','YTick','','XColor','w','YColor','w'); title('Sens_R (LvNG)');
caxis([-1 1]*2);

subplot(3,2,6);
% imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal'); set(imX,'alphadata',0.7);
% hold on;
s1=scatter(ml,ap,dotSize,out(:,6),dotShape,'filled'); axis equal; colorbar;
axis([-5 5 -6 4]); set(gca,'box','off','Xtick','','YTick','','XColor','w','YColor','w'); title('Sens_R (RvNG)');
caxis([-1 1]*2);%s1.MarkerEdgeColor=[1 1 1]*0.8;

cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
    linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
colormap(cmap);
    set(gcf,'Color','w');

set(get(gcf,'children'),'box','off','fontsize',21)

%% LASER EFFECT MAP: performance effects of the laser (no model)
figure;
dotShape = 's';
kimg=imread('\\basket.cortexlab.net\home\stuff\kirkcaldie_brain_BW.PNG');
clear PD;
for s = 1:numSubjects
    subplot(1,numSubjects,s);
    imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal');
    set(imX,'alphadata',0.7);
    hold on;
    
    PF = 100*pivottable(l(s).data.laserIdx,[],l(s).data.feedbackType==1,'mean');
    
    % NULL TEST
    tab = tabulate(l(s).data.laserIdx);
    proportion = mean(tab(2:end,3)/100);
    NUM_SHUFFLES = 1000;
    nullP = nan(NUM_SHUFFLES,1);
    nonL_idx = l(s).data.laserIdx==0;
    
    for shuff = 1:NUM_SHUFFLES
        disp(shuff);
        L_idx = l(s).data.laserIdx;
        nL = randsample(find(nonL_idx),round(proportion*sum(nonL_idx)));
        L_idx(nL) = max(l(s).data.laserIdx)+1;
        
        shuf_p = 100*pivottable(L_idx,[],l(s).data.feedbackType==1,'mean');
        nullP(shuff) = shuf_p(end);
    end
    
    nullP_ci95 = quantile(nullP,[0.025 0.975]);
    nullP_ci99 = quantile(nullP,[0.005 0.995]);
   
    
    ap = l(s).inactivationSite(:,1);
    ml = l(s).inactivationSite(:,2);
    out = PF(2:end);
    if bilateral_mirror == 1
        ap = [ap;ap];
        ml = [ml;-ml];
        out = [out;out];
    end
    
    %sig testing
    sig95 = double((out < nullP_ci95(1)) | (nullP_ci95(2) < out));
    sig99 = double((out < nullP_ci99(1)) | (nullP_ci99(2) < out));
    sig = sum([sig95 sig99],2);
    sig(sig==0)=0.2;
    
    s1=scatter(ml,ap,sig*80,out,dotShape,'filled'); axis equal; colorbar;
    %         hold on;
    s1.MarkerEdgeColor=[1 1 1]*0.9;
        
%     s1=scatter(ml,ap,150,out,'s','filled'); axis equal; colorbar;
    ylim([-6 4]); %s1.MarkerEdgeColor=[1 1 1]*0.8;
    cax = PF(1) + [-1 1]*max(abs(PF(2:end)-PF(1)));
    try
        PD(:,s)=PF(2:end)-PF(1);
    catch
    end
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

figure;
% imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal');
% set(imX,'alphadata',0.7);
% hold on;

        ap = l(s).inactivationSite(:,1);
        ml = l(s).inactivationSite(:,2);
        out = mean(PD,2);
        if bilateral_mirror == 1
            ap = [ap;ap];
            ml = [ml;-ml];
            out = [out;out];
        end

s1=scatter(ml,ap,150,out,'s','filled'); axis equal; colorbar;
ylim([-6 4]); %s1.MarkerEdgeColor=[1 1 1]*0.8;
cax = [-1 1]*max(abs(mean(PD,2)));
caxis(cax);
hold off;
xlim([-5 5]);
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
figure; RTs_all=[];
for s = 1:numSubjects
%     PF = pivottable(l(s).data.laserIdx,[],l(s).data.feedbackType==1,'mean');
%     PF = PF(2:end)-PF(1);
%     PF = [PF nan(length(PF),1)];
    numSessions = max(l(s).data.sessionID);
    rt_z = nan(length(l(s).data.sessionID),1);
    for sessionID = 1:numSessions
        for choice = 1:2
            idx = l(s).data.sessionID==sessionID & l(s).data.response==choice;
            rt_z(idx) = zscore(l(s).data.RT(idx));
        end
    end
    RTs = pivottable(l(s).data.laserIdx,l(s).data.response,rt_z,'mean');
%     RTs = bsxfun(@minus,RTs(2:end,:),RTs(1,:));
    RTs = RTs(:,1:2);
    
%     Ls = pivottable(l(s).data.laserIdx,l(s).data.response,l(s).data.lickenergy,'median');
%     Ls = bsxfun(@minus,Ls(2:end,:),Ls(1,:));
%     Ls = Ls(:,1:2);
%     


% NULL TEST
    NUM_SHUFFLES = 1000;
    nullz = nan(NUM_SHUFFLES,2);

    for lr = 1:2
        tab = tabulate(l(s).data.laserIdx(l(s).data.response==lr));
        proportion = mean(tab(2:end,3)/100);
        nonL_idx = l(s).data.laserIdx==0 & l(s).data.response==lr;

        for shuff = 1:NUM_SHUFFLES
            disp(shuff);
            L_idx = l(s).data.laserIdx;
            nL = randsample(find(nonL_idx),round(proportion*sum(nonL_idx)));
            L_idx(nL) = max(l(s).data.laserIdx)+1;
            
            shuf_rtz = pivottable(L_idx,l(s).data.response,rt_z,'mean');
            nullz(shuff,lr) = shuf_rtz(end,lr);
        end
    end
    nullz_ci95 = quantile(nullz,[0.025 0.975]);
    nullz_ci99 = quantile(nullz,[0.005 0.995]);
    
    dotSize=150;
    dotShape='s';
    kimg=imread('D:\kirkcaldie_brain_BW.PNG');
    
    a=1;
    
    side={'left choices','right choices'};
    for lr=1:2
        subplot(2,numSubjects,s + (lr-1)*numSubjects);
%         imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal');
%         set(imX,'alphadata',0.7);
%         hold on;
            
        ap = l(s).inactivationSite(:,1);
        ml = l(s).inactivationSite(:,2);
        out = RTs(2:end,lr);
        if bilateral_mirror == 1
            ap = [ap;ap];
            ml = [ml;-ml];
            out = [out;out];
        end
        
        %sig testing
        sig95 = double((out < nullz_ci95(1,lr)) | (nullz_ci95(2,lr) < out));
        sig99 = double((out < nullz_ci99(1,lr)) | (nullz_ci99(2,lr) < out));
        sig = sum([sig95 sig99],2);
        sig(sig==0)=0.2; 

        s1=scatter(ml,ap,sig*80,out,dotShape,'filled'); axis equal; colorbar;
%         hold on;  
        s1.MarkerEdgeColor=[1 1 1]*0.9;
%         s2=scatter(ml,ap,sig95*150*0.4,out,dotShape,'filled'); axis equal; colorbar;
%         hold off; s2.MarkerEdgeColor=[1 1 1]*0.9;
                
                
%         s1=scatter(ml,ap,dotSize,out,dotShape,'filled'); axis equal; colorbar;
        ylim([-6 4]);% s1.MarkerEdgeColor=[1 1 1]*0.8;
%         
%         cax = RTs(1,lr) + [-1 1]*max(abs(RTs(2:end,lr)-RTs(1,lr)));
%         caxis(cax);
        caxis(gca,[-1 1]*2);
        
        %         set(s(lr),'markeredgecolor',[1 1 1]*1,'linewidth',0);
        hold off;
        title(expRefs{s,1});
        xlim([-5 5]);
        a=a+1;
        set(gca,'xtick','','ytick','');
        set(gca,'box','off','XColor','w','YColor','w');
        
        if s == 1
            if lr == 1
                ylabel('RT to left choices');
            else
                ylabel('RT to right choices');
            end
        end
    end
    %
    cmap = [ linspace(0,1,100)' ones(100,1) linspace(0,1,100)';
        linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
    colormap(flipud(cmap));
    
    try
    RTs_all(:,:,s)=RTs;
    catch
    end
end
set(gcf,'Color','w');

figure;
r = mean(RTs_all,3);
for i = 1:2
    subplot(2,1,i);
%     imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal');
%     set(imX,'alphadata',0.7);
%     hold on;
    ap = l(s).inactivationSite(:,1);
    ml = l(s).inactivationSite(:,2);
    out = r(2:end,i);
    if bilateral_mirror == 1
        ap = [ap;ap];
        ml = [ml;-ml];
        out = [out;out];
    end

    s1=scatter(ml,ap,dotSize,out,dotShape,'filled'); axis equal; colorbar;
    ylim([-6 4]);        xlim([-5 5]); %s1.MarkerEdgeColor=[1 1 1]*0.8;
    caxis(gca,[-1 1]*2);
%     caxis(cax);
    hold off;
    set(gca,'xtick','','ytick','');
    set(gca,'box','off','XColor','w','YColor','w');
end
colormap(flipud(cmap)); set(gcf,'color','w');

%% LASER EFFECT MAP:  map of change in % choice over the laser sites (no model)
pDiff = []; pDiff_all = [];
        stim_labels = {'C=0','High CL','High CR','High CL and CR'};
        NUM_SHUFFLES = 1000;

nullP = nan(NUM_SHUFFLES,3,numSubjects,4);

for s = 1:numSubjects
    maxC = max(l(s).data.contrast_cond(:));
    for type = 1:4
        if type==1
            STIM_IDX = diff(l(s).data.contrast_cond,[],2)==0 & sum(l(s).data.contrast_cond,2)==0; %cl=cr=0
        elseif type==2
            STIM_IDX = l(s).data.contrast_cond(:,1)==maxC & l(s).data.contrast_cond(:,2)==0; %cl=high
        elseif type==3
            STIM_IDX = l(s).data.contrast_cond(:,2)==maxC & l(s).data.contrast_cond(:,1)==0; %cr=high
        elseif type==4
            STIM_IDX = l(s).data.contrast_cond(:,2)==maxC & l(s).data.contrast_cond(:,1)==maxC; %cr=cr=high
        end
        
        
        
        pDiff = [];
        for session = 1:length(expRefs{s,2})
            idx = l(s).data.laserIdx==0 & l(s).data.sessionID==session & STIM_IDX;
            noLaser_responses = l(s).data.response(idx);
            noLaser_p = 100*arrayfun(@(r)(sum(noLaser_responses==r)/length(noLaser_responses)),1:3);
            
            Laser_p=[];
            for site = 1:size(l(s).inactivationSite,1)
                idx = l(s).data.laserIdx==site & l(s).data.sessionID==session & STIM_IDX;
                Laser_responses = l(s).data.response(idx);
                Laser_p(site,:) = 100*arrayfun(@(r)(sum(Laser_responses==r)/length(Laser_responses)),1:3);
            end
            
            pDiff(:,:,session) = bsxfun(@minus,Laser_p,noLaser_p);
        end
        pDiff_all(:,:,s,type) = nanmean(pDiff,3);
        
%         
%         %NULL
%         tab = tabulate(l(s).data.laserIdx(l(s).data.response==lr));
%         proportion = mean(tab(2:end,3)/100);
%         nonL_idx = l(s).data.laserIdx==0 & l(s).data.response==lr;
%         for shuff = 1:NUM_SHUFFLES
%             tab = tabulate(l(s).data.laserIdx(l(s).data.response==lr));
%             proportion = mean(tab(2:end,3)/100);
%             nonL_idx = l(s).data.laserIdx==0 & l(s).data.response==lr;
% 
%             disp(shuff);
%             L_idx = l(s).data.laserIdx;
%             nL = randsample(find(nonL_idx),round(proportion*sum(nonL_idx)));
%             L_idx(nL) = max(l(s).data.laserIdx)+1;
%             
%             shuf_rtz = pivottable(L_idx,l(s).data.response,rt_z,'mean');
%             nullz(shuff,lr) = shuf_rtz(end,lr);
%         
%         end
    end
    
    
    figure('name',expRefs{s,1});
    titles = {'\Delta pL','\Delta pR','\Delta pNG'};
    for type = 1:4
        for lr=1:3
            pd = pDiff_all(:,:,s,type);
            subplot(4,3,3*type - 3 + lr);
            
            ap = l(s).inactivationSite(:,1);
            ml = l(s).inactivationSite(:,2);
            out = pd(:,lr);
            if bilateral_mirror == 1
                ap = [ap;ap];
                ml = [ml;-ml];
                out = [out;out];
            end
            
            scatter(ml,ap,300,out,'s','filled'); axis equal; colorbar;
            ylim([-6 4]);
            xlim([-4 4]);
            caxis([-1 1]*40);
            set(gca,'xtick','','ytick','','xcolor','w','ycolor','w');
            if type == 1
                title(titles{lr});
            end
        end
    end
    set(gcf,'color','w');
cmap = [ linspace(0,1,100)' ones(100,1) linspace(0,1,100)';
        linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
    colormap(flipud(cmap));
end

figure('name','average over mice');
pda = mean(pDiff_all,3);
for type=1:4
    for lr=1:3
        subplot(4,3,3*type - 3 + lr)
        
        ap = l(1).inactivationSite(:,1);
        ml = l(1).inactivationSite(:,2);
        out = pda(:,lr,1,type);
        if bilateral_mirror == 1
            ap = [ap;ap];
            ml = [ml;-ml];
            out = [out;out];
        end
        
        scatter(ml,ap,200,out,'s','filled'); axis equal; %colorbar;
        ylim([-6 4]);
        xlim([-4 4]);
        caxis([-1 1]*40);
        
        if type == 1
        title(titles{lr});
        end
        
        set(gca,'box','off','xtick','','ytick','','XColor','w','YColor','w');
        
%         if lr==1
%             ylabel(stim_labels{type});
%         end
    end
end
set(gcf,'color','w');

cmap = [ linspace(0,1,100)' ones(100,1) linspace(0,1,100)';
        linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
    colormap(flipud(cmap));

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

%% Check whether model is fitting well on non-laser trials %NEED TO ADD GLOBAL INTERCEPT
g = g.setModel('C50');
for s = 1:numSubjects
    figure('name',expRefs{s,1});
    numSessions = max(l(s).data.sessionID);
    for session = 1:numSessions
        D = getrow(l(s).data,l(s).data.laserIdx==0 & l(s).data.sessionID==session);
        params = [sess(s).Biases(session,1) + sess(s).GlobalIntercept(1),...
            sess(s).CLeft(session,1),...
            sess(s).CRight(session,1),...
            sess(s).Biases(session,2) + sess(s).GlobalIntercept(2),...
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
        
        params = [sess(s).Biases(session,1) + sess(s).GlobalIntercept(1),...
            sess(s).CLeft(session,1),...
            sess(s).CRight(session,1),...
            sess(s).Biases(session,2) + sess(s).GlobalIntercept(2),...
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
figure;
for s = 1:numSubjects
    D = getrow(l(s).data,l(s).data.laserIdx==0);
    numSessions = max(D.sessionID);
    
    XnL = X{s}(l(s).data.laserIdx==0,1:(3*numSessions));
    YnL = Y{s}(l(s).data.laserIdx==0);

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
        
%         %fit
        fit=glmnet(Xb,Yb,'multinomial',glmnetSet(opts));
        b=glmnetCoef(fit,0.01);
        b=[b{1}-b{3} b{2}-b{3}];
        b(1,:) = [];

        bootstrap.Biases(:,:,iter) = b(1:3:end,:);
        bootstrap.CLeft(:,:,iter) = b(2:3:end,:);
        bootstrap.CRight(:,:,iter) = b(3:3:end,:);
%

%         %model output
%         D.y_hat = glmnetPredict(fit,XnL,0,'response');
%         for session = 1:numSessions
%             E = getrow(D,D.sessionID==session);
%             stim = cell(1,4);
%             stim{1} = sum(E.contrast_cond,2)==0;
%             stim{2} = ~stim{1} & (E.contrast_cond(:,1)==E.contrast_cond(:,2));
%             stim{3} = diff(E.contrast_cond,[],2)<0;
%             stim{4} = diff(E.contrast_cond,[],2)>0;
%             
%             for stim_type = 1:4
%                 r = E.response(stim{stim_type});
%                 model_output(iter,:,stim_type,session)=mean(E.y_hat(stim{stim_type},:));
%             end
%         end
    end
    
    
    biases = permute(bootstrap.Biases,[3 2 1]);
    sensL = permute(bootstrap.CLeft,[3 2 1]);
    sensR = permute(bootstrap.CRight,[3 2 1]);
    
    labels_old = {'b_{LvNG}','s_{LvNG}','s_{RvNG}','b_{RvNG}'};
%     
    figure;
    for session = 1:numSessions
        old_params = [biases(:,1,session) sensL(:,1,session) sensR(:,1,session) sensR(:,2,session) sensL(:,2,session) biases(:,2,session)];

        subplot(numSubjects,numSessions,(s-1)*numSessions + session);
        imagesc(corrcoef(old_params));
        if session == 1
            set(gca,'XTickLabel',labels_old,'Xtick',1:6,'YTickLabel',labels_old,'Ytick',1:6,'XTickLabelRotation',90,'fontsize',9);
            %             title('old parameterisation');
            ylabel(expRefs{s,1});
        else
            set(gca,'XTickLabel','','YTickLabel','');
        end
        axis square; caxis([-1 1]);

        cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
            linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
        colormap(cmap);

    end
    
%     for session = 1:numSessions
%         figure('name',[expRefs{s,1} ' session ' num2str(session) ' CONVENTIONAL PARAMETERS']);
%         gplotmatrix([biases(:,1,session) sensit(:,1,session) sensit(:,2,session) biases(:,2,session)],[],[],'b','.',5,'on','none',labels_old);
%         
% %         figure('name',[expRefs{s,1} ' session ' num2str(session) ' NEW PARAMETERS']);
% %         gplotmatrix([czero(:,1,session) splus(:,1,session) sminus(:,1,session) sminus(:,2,session) splus(:,2,session) czero(:,2,session)],[],[],'b','.',5,'on','none',labels_new);
% %         drawnow;
%     end
%     
%     
%     biases = permute(bootstrap.Biases,[3 2 1]);
%     sensL = permute(bootstrap.CLeft,[3 2 1]);
%     sensR = permute(bootstrap.CRight,[3 2 1]);
% 
%         %reparameterise the parameters at lambda=0
%     splus = (sensL + sensR)/2;
%     sminus = (sensL - sensR)/2;
%     czero = -biases./splus;
%     
%     labels_old = {'b_{LvNG}','sL_{LvNG}','sR_{LvNG}','sR_{RvNG}','sL_{RvNG}','b_{RvNG}'};
%     labels_new = {'czero_{LvNG}','splus_{LvNG}','sminus_{LvNG}','sminus_{RvNG}','splus_{RvNG}','czero_{RvNG}'};
%     
%     figure('name',expRefs{s,1});
%     for session = 1:numSessions
%         old_params = [biases(:,1,session) sensL(:,1,session) sensR(:,1,session) sensR(:,2,session) sensL(:,2,session) biases(:,2,session)];
%         new_params = [czero(:,1,session) splus(:,1,session) sminus(:,1,session) sminus(:,2,session) splus(:,2,session) czero(:,2,session)];
%         
%         subplot(2,numSessions,session);
%         imagesc(corrcoef(old_params));
%         if session == 1
%             set(gca,'XTickLabel',labels_old,'Xtick',1:6,'YTickLabel',labels_old,'Ytick',1:6,'XTickLabelRotation',90,'fontsize',9);
%             title('old parameterisation');
%         else
%             set(gca,'XTickLabel','','YTickLabel','');
%         end
%         axis square; caxis([-1 1]);
%         subplot(2,numSessions,session + numSessions);
%         imagesc(corrcoef(new_params));
%        if session == 1
%             set(gca,'XTickLabel',labels_new,'Xtick',1:6,'YTickLabel',labels_new,'Ytick',1:6,'XTickLabelRotation',90,'fontsize',9);
%             title('new parameterisation');
%         else
%             set(gca,'XTickLabel','','YTickLabel','');
%         end
%         axis square;caxis([-1 1]);
%         
% cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
%         linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
%     colormap(cmap);
% %         figure('name',[expRefs{s,1} ' session ' num2str(session) ' CONVENTIONAL PARAMETERS']);
% %         gplotmatrix([biases(:,1,session) sensL(:,1,session) sensR(:,1,session) sensR(:,2,session) sensL(:,2,session) biases(:,2,session)],[],[],'b','.',5,'on','none',labels_old);
% %         
% %         figure('name',[expRefs{s,1} ' session ' num2str(session) ' NEW PARAMETERS']);
% %         gplotmatrix([czero(:,1,session) splus(:,1,session) sminus(:,1,session) sminus(:,2,session) splus(:,2,session) czero(:,2,session)],[],[],'b','.',5,'on','none',labels_new);
% %         drawnow;
%     end
    
%     f=figure('name',expRefs{s,1});
%     for session = 1:numSessions 
%         for stim_type = 1:4
%             subplot(numSessions,4,4*session - 4 + stim_type);
%             plot(1:3,model_output(:,:,stim_type,session),'b.'); ylim([0 1]); xlim([0.5 3.5]);
% %             boxplot(model_output(:,[1 3 2],stim_type,session),'plotstyle','compact'); ylim([0 1]);
%             if session == 1
%                 title(stim_labels{stim_type});
%             end
%         end
%     end
%     ax=get(f,'children');
%     set(ax,'Xticklabel','','Yticklabel','','box','off');
%     set(ax(1),'XTickLabel',{'pL','pR','pNG'});
%     ylabel(ax(1),'probability');
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
model_labels = {'Orig MNR (One stage with intercept) L=cv',...
                'Orig MNR (Two stage) L=cv',...
                'Hier MNR (Combined sessions + Laser biasing L=0',...
                'Hier MNR (Separate sessions + Laser biasing L=0'};
for s = 1:numSubjects
    numSessions = max(l(s).data.sessionID);
    cv = cvpartition(size(X{s},1),'kfold',5);
    
    p_hats = nan(cv.NumObservations,4);
    for iter = 1:cv.NumTestSets
        disp(iter);
        disp(iter);
        trainX = X{s}(cv.training(iter),:);
        trainY = Y{s}(cv.training(iter));
        
        testX = X{s}(cv.test(iter),:);
        testY = Y{s}(cv.test(iter));
        
%         %one stage fitting without a global intercept
%         opts.intr=0;
%         fit=cvglmnet(trainX,trainY,'multinomial',glmnetSet(opts));
%         p = cvglmnetPredict(fit,testX,'lambda_min','response');
%         p_hats(cv.test(iter),1) = p(sub2ind(size(p), [1:length(testY)]', testY));
%         
        %one stage fitting with a global intercept: no laser terms
        opts.intr=1;
        fit = cvglmnet(trainX(:,1:numSessions*3),trainY,'multinomial',glmnetSet(opts));
        p = cvglmnetPredict(fit,testX(:,1:numSessions*3),'lambda_min','response');
        p_hats(cv.test(iter),1) = p(sub2ind(size(p), [1:length(testY)]', testY));
%         
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
        p_hats(cv.test(iter),2) = p(sub2ind(size(p), [1:length(testY)]', testY));
        
        %Hier MNR with combined sessions + laser biasing only
        trainX2 = X2{s}(cv.training(iter),:);
        trainY2 = Y2{s}(cv.training(iter));
        testX2 = X2{s}(cv.test(iter),:);
        testY2 = Y2{s}(cv.test(iter));
        
        [b,dev,stat]=mnrfit(trainX2(:,2:end),trainY2,'model','hierarchical');
        p = mnrval(b,testX2(:,2:end));
        p_hats(cv.test(iter),3) = p(sub2ind(size(p), [1:length(testY2)]', testY2));

        %Hier MNR with separate sessions + laser biasing only
        trainX2 = [X{s}(cv.training(iter),1:3*numSessions) X{s}(cv.training(iter),3*numSessions+1 :3: end)];
        trainY2 = Y2{s}(cv.training(iter));
        testX2 = [X{s}(cv.test(iter),1:3*numSessions) X{s}(cv.test(iter),3*numSessions+1 :3: end)];
        testY2 = Y2{s}(cv.test(iter));
        
        [b,dev,stat]=mnrfit(trainX2(:,2:end),trainY2,'model','hierarchical');
        p = mnrval(b,testX2(:,2:end));
        p_hats(cv.test(iter),4) = p(sub2ind(size(p), [1:length(testY2)]', testY2));
    end
    
    tab=tabulate(Y{s}); tab=tab(:,3)/100;
    guess_bpt = sum(tab.*log2(tab));
    
    cv_gof(s,:) = mean(log2(p_hats))-guess_bpt;
end 
disp('done');

bar(cv_gof); ylabel('loglik [bits] relative to guessing @~-1.5'); legend(model_labels);
set(gca,'XTickLabel',expRefs(:,1),'xtick',1:numSubjects);

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

%% Crossval C50 vs C50-subset
figure;
for s = 1:numSubjects
   numSessions =  max(l(s).data.sessionID);
   subplot(1,numSubjects,s);
   
   ll = [];
   for se = 1:numSessions
       D = getrow(l(s).data,l(s).data.laserIdx==0 & l(s).data.sessionID==se);
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