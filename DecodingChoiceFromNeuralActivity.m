%% 2AUC data


cfn = @(c,n,c50)((c.^n)./(c.^n + c50^n));
lambda = 'lambda_1se'; 
lambda = 0.05; 
numFolds = 10;
numT = 81; %number of time bins to sample neural activity from within an epoch
interval = [-1 1];
for session = 1:8
    [spikeTimes,behav,behavT,o]=getSpikingData('nick laveran',session);

    guess_bpt = [];
    tab = tabulate(behav.response);
    tab = tab(:,3)/100;
    guess_bpt(1)=sum(tab.*log2(tab));
    
    tab = tabulate(behav.response==3);
    tab = cell2mat(tab(:,3))/100;
    guess_bpt(2)=sum(tab.*log2(tab));
    
    tab = tabulate(behav.response(behav.response<3));
    tab = tab(:,3)/100;
    guess_bpt(3)=sum(tab.*log2(tab));

    g = GLM(behav);
    
    %Fit GLM to find c50/n parameters
    g = g.setModel('C50-subset').fit;
    
    %Construct design matrix
    X = sparse(cfn(behav.contrast_cond,g.parameterFits(5),g.parameterFits(6)));
    Y = behav.response;
    
    %Fit MNR model with global intercept
    clear opts;
    opts.intr = 1;
%     opts.alpha = 0.5;
    fit = cvglmnet(X,Y,'multinomial',glmnetSet(opts));
    
    %Calculate offset which will be supplied to the model when adding
    %neural activity
    offset_MNR = cvglmnetPredict(fit,X,[],'link');
    
    %Calculate baseline non-neural model performance
    cv = cvpartition(size(X,1),'kfold',numFolds);
    phat_baseline = nan(cv.NumObservations,3);
    for fold = 1:cv.NumTestSets
        trainX = X(cv.training(fold),:);
        trainY = Y(cv.training(fold));
        testX = X(cv.test(fold),:);
        testY = Y(cv.test(fold));
        
        fitcv = cvglmnet(trainX,trainY,'multinomial',glmnetSet(opts));
        p = cvglmnetPredict(fitcv,testX,lambda,'response');
        
        %Orig MNR
        phat_baseline(cv.test(fold),1) = p(sub2ind(size(p), [1:length(testY)]', testY));
        
        %Go vs NoGo decoder
        pGO = 1-p(:,3);
        phat_baseline(cv.test(fold),2) = (testY<3).*pGO + (testY==3).*(1-pGO);
        
        %L vs R decoder
        pL_G = p(:,1)./pGO;
        pR_G = p(:,2)./pGO;
        obs = (testY==1).*pL_G + (testY==2).*pR_G; obs(obs==0)=nan;
        phat_baseline(cv.test(fold),3) = obs;
    end
    baseline_bpt = nanmean(log2(phat_baseline));
    
    figure('name',o{1});
    %Second: Add the activity from all neurons at a given time to the model
    for e = 1:length(behavT.timestamps) %for each epoch
        t_zero = behavT.(behavT.timestamps{e});
        tsteps = linspace(interval(1),interval(2),numT);
            
        ax(e)=subplot(1,length(behavT.timestamps),e);
        bpt=[];
        for t = 1:length(tsteps) %for each timestep in that epoch
            firingRate = nan(length(t_zero),length(spikeTimes)); %num trials x num neurone
            queryT = t_zero + tsteps(t);
            
            for neuron = 1:length(spikeTimes) %get activity for all neurons
                for trial = 1:length(t_zero)
                    firingRate(trial,neuron) =sum(((queryT(trial)-0.05) < spikeTimes{neuron}) & (spikeTimes{neuron} < (queryT(trial)+0.05))) / 0.1;
                end
            end
            
            %Centre firing Rates
            firingRate = bsxfun(@minus,firingRate,mean(firingRate,1));
            
            %Fit model
            clear opts;
            opts.intr=1;
%             opts.alpha = 0.5;
            cv = cvpartition(size(firingRate,1),'kfold',numFolds);
            phat = nan(cv.NumObservations,3);
            for fold = 1:cv.NumTestSets
                trainX = firingRate(cv.training(fold),:);
                trainY = Y(cv.training(fold));
                testX = firingRate(cv.test(fold),:);
                testY = Y(cv.test(fold));
                
                %ORIGINAL MNR
                opts.offset = offset_MNR(cv.training(fold),:);
                fit = glmnet(trainX,trainY,'multinomial',glmnetSet(opts));
                p = glmnetPredict(fit,testX,lambda,'response',[],offset_MNR(cv.test(fold),:));
                phat(cv.test(fold),1) = p(sub2ind(size(p), [1:length(testY)]', testY));
                
                %Go vs NoGo decoder
                pGO = 1-p(:,3);
                phat(cv.test(fold),2) = (testY<3).*pGO + (testY==3).*(1-pGO);
                
                %L VS R decoder
                pL_G = p(:,1)./pGO;
                pR_G = p(:,2)./pGO;
                obs = (testY==1).*pL_G + (testY==2).*pR_G; obs(obs==0)=nan;
                phat(cv.test(fold),3) = obs;
                
            end
            
%             bpt(t,:) = nanmean(log2(phat)) - guess_bpt;
            bpt(t,:) = nanmean(log2(phat)) - baseline_bpt;
            
            plot(ax(e),tsteps(1:length(bpt)),bpt,'LineWidth',2); drawnow; 
        end
        set(gca,'box','off');
        if e == 1
%             ylabel('log_2 likelihood relative to guessing [bits]');
            ylabel('log_2 likelihood relative to model baseline [bits]');

            legend('Full choice decoder','G vs NG decoder','L v R decoder'); legend boxoff;
            title(o{1});
        else
            set(gca,'ytick','','ycolor','w'); ylabel('');
        end
        xlabel(behavT.timestamps{e}); ylim([-0.2 0.8]); 
    end
    linkaxes(ax,'y');
    set(gcf,'color','w'); 
    savefig(gcf,['\\basket.cortexlab.net\home\figures\GLM+NeuralActivity\glmnet_decoders\2AUC_' o{1}])
    
end

%% 2AFC data
cfn = @(c,n)(c.^n);
lambda = 'lambda_1se'; 
lambda = 0.05; 
numFolds = 10;
numT = 101; %number of time bins to sample neural activity from within an epoch
interval = [-1 1];

for session = 1:8
    [spikeTimes,behav,behavT,o]=getSpikingData('mush v1',session);

    guess_bpt = [];
    tab = tabulate(behav.response);
    tab = tab(:,3)/100;
    guess_bpt=sum(tab.*log2(tab));
    
    g = GLM(behav);
    
    %Fit GLM to find c50/n parameters
    g = g.setModel('C^N-subset-2AFC').fit;
    
    %Construct design matrix
    X = sparse(cfn(behav.contrast_cond,g.parameterFits(4)));
    Y = behav.response;
    
    %Fit MNR model with global intercept
    clear opts;
    opts.intr = 1;
    fit = cvglmnet(X,Y,'binomial',glmnetSet(opts));
    
    %Calculate offset which will be supplied to the model when adding
    %neural activity
    offset_MNR = cvglmnetPredict(fit,X,[],'link');
    
    %Calculate baseline non-neural model performance
    cv = cvpartition(size(X,1),'kfold',numFolds);
    phat_baseline = nan(cv.NumObservations,1);
    for fold = 1:cv.NumTestSets
        trainX = X(cv.training(fold),:);
        trainY = Y(cv.training(fold));
        testX = X(cv.test(fold),:);
        testY = Y(cv.test(fold));
        
        fitcv = cvglmnet(trainX,trainY,'binomial',glmnetSet(opts));
        pR = cvglmnetPredict(fitcv,testX,lambda,'response');
        %Orig MNR
        phat_baseline(cv.test(fold),1) = (testY==1).*(1-pR) + (testY==2).*pR;
    end
    baseline_bpt = nanmean(log2(phat_baseline));
    
    figure('name',o{1});
    %Second: Add the activity from all neurons at a given time to the model
    for e = 1:length(behavT.timestamps) %for each epoch
        t_zero = behavT.(behavT.timestamps{e});
        tsteps = linspace(interval(1),interval(2),numT);
  
        ax(e)=subplot(1,length(behavT.timestamps),e);
        bpt=[];
        for t = 1:length(tsteps) %for each timestep in that epoch
            firingRate = nan(length(t_zero),length(spikeTimes)); %num trials x num neurone
            queryT = t_zero + tsteps(t);
            
            for neuron = 1:length(spikeTimes) %get activity for all neurons
                for trial = 1:length(t_zero)
                    firingRate(trial,neuron) =sum(((queryT(trial)-0.05) < spikeTimes{neuron}) & (spikeTimes{neuron} < (queryT(trial)+0.05))) / 0.1;
                end
            end
%             
%             if tsteps(t)>0 %Test if adding choice to firing rate is
%             detected in my decoder. (It is)
%                 firingRate(:,1) = Y;
%             end
            
            %Centre firing Rates
            firingRate = bsxfun(@minus,firingRate,mean(firingRate,1));
            
            %Fit model
            clear opts;
            opts.intr=1; 
            cv = cvpartition(size(firingRate,1),'kfold',numFolds);
            phat = nan(cv.NumObservations,1);
            for fold = 1:cv.NumTestSets
                trainX = firingRate(cv.training(fold),:);
                trainY = Y(cv.training(fold));
                testX = firingRate(cv.test(fold),:);
                testY = Y(cv.test(fold));
                
                %ORIGINAL MNR
                opts.offset = offset_MNR(cv.training(fold),:);
                fit = glmnet(trainX,trainY,'binomial',glmnetSet(opts));
                pR = glmnetPredict(fit,testX,lambda,'response',[],offset_MNR(cv.test(fold),:));
                phat(cv.test(fold),1) = (testY==1).*(1-pR) + (testY==2).*pR;
            end
            
%             bpt(t,:) = nanmean(log2(phat)) - guess_bpt;
            bpt(t,:) = nanmean(log2(phat)) - baseline_bpt;
            
            plot(ax(e),tsteps(1:length(bpt)),bpt,'LineWidth',2); drawnow; 
        end
        
        set(gca,'box','off');
        if e == 1
%             ylabel('log_2 likelihood relative to guessing [bits]');
            ylabel('log_2 likelihood relative to model baseline [bits]');

            legend('Full choice decoder'); legend boxoff;
            title(o{1});
        else
            set(gca,'ytick','','ycolor','w'); ylabel('');
        end
        xlabel(behavT.timestamps{e}); ylim([-0.2 0.6]); 
    end
    set(gcf,'color','w'); linkaxes(ax,'y');
    savefig(gcf,['\\basket.cortexlab.net\home\figures\GLM+NeuralActivity\glmnet_decoders\2AFC_' o{1}])
    
end
