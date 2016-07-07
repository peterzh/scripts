dirs = {'\\basket.cortexlab.net\data\nick\M140528_NS1\20141130\cw_1','M140528_NS1 1';
        '\\basket.cortexlab.net\data\nick\M140528_NS1\20141201\cw_1','M140528_NS1 2';
        '\\basket.cortexlab.net\data\nick\M140528_NS1\20141202',     'M140528_NS1 3';
        '\\basket.cortexlab.net\data\nick\M150218_NS1LAV\20150529\cw_1','Laveran Left PPC 1';
        '\\basket.cortexlab.net\data\nick\M150218_NS1LAV\20150530\cw_1','Laveran Left Cingulate 1';
        '\\basket.cortexlab.net\data\nick\M150218_NS1LAV\20150601\cw_1','Laveran Right Cingulate 1';
        '\\basket.cortexlab.net\data\nick\M150218_NS1LAV\20150602\cw_1','Laveran Left V1 1';
        '\\basket.cortexlab.net\data\nick\M150218_NS1LAV\20150528\cw_1','Laveran Left Cingulate 2'};
        };
    
dirs = {};
    
dirs = {};

    
cfn = @(c,n,c50)((c.^n)./(c.^n + c50^n));
lambda = 'lambda_1se'; 
lambda = 0.05; 
numFolds = 5;
numT = 51; %number of time bins to sample neural activity from within an epoch

for session = 4:length(dirs)
    n=neurGLM(dirs{session,1},dirs{session,2});
    neurons = n.spikes.shanks.singleUnitInds;
    
    %First: Construct behav model
    %Put into GLM object
    D = struct;
    D.contrast_cond = [n.cwLabels.contrastLeft n.cwLabels.contrastRight];
    D.response = n.cwLabels.responseMade;
    D.repeatNum = n.cwLabels.repeatNum;
    D.feedbackType = n.cwLabels.feedbackType;
%     D.RT = n.cwLabels.reactionTime;
    
    guess_bpt = [];
    tab = tabulate(D.response);
    tab = tab(:,3)/100;
    guess_bpt(1)=sum(tab.*log2(tab));
    
    tab = tabulate(D.response==3);
    tab = cell2mat(tab(:,3))/100;
    guess_bpt(2)=sum(tab.*log2(tab));
    
    tab = tabulate(D.response(D.response<3));
    tab = tab(:,3)/100;
    guess_bpt(3)=sum(tab.*log2(tab));
     
    g = GLM(D);
    
    %Fit GLM to find c50/n parameters
    g = g.setModel('C50-subset').fit;
    
    %Construct design matrix
    X = sparse(cfn(D.contrast_cond,g.parameterFits(5),g.parameterFits(6)));
    Y = D.response;
    
    %Fit MNR model with global intercept
    clear opts;
    opts.intr = 1;
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
    
%     %Also run a hierarchical model on the behav data
%     fit1 = cvglmnet(X,Y==3,'binomial',glmnetSet(opts));%NG vs GO fits
%     offset_NGvG = cvglmnetPredict(fit1,X,[],'link');
%     fit2 = cvglmnet(X(Y==1 | Y==2,:),Y(Y==1 | Y==2),'binomial',glmnetSet(opts)); %L v R fits
%     offset_LvR = cvglmnetPredict(fit2,X(Y==1 | Y==2,:),[],'link');

    figure('name',dirs{session,2});
    %Second: Add the activity from all neurons at a given time to the model
    for e = 1:size(n.epochs,1) %for each epoch
        t_zero = n.cwEvents.(n.epochs{e,1}); %Time of the epoch beginning
        interval = n.epochs{e,2};
        tsteps = linspace(interval(1),interval(2),numT);
            
        ax(e)=subplot(1,size(n.epochs,1),e);
        bpt=[];
        for t = 1:length(tsteps) %for each timestep in that epoch
            firingRate = nan(length(t_zero),length(neurons)); %num trials x num neurone
            queryT = t_zero + tsteps(t);
            
            for neuron = 1:length(neurons) %get activity for all neurons
                spikeTimes = n.spikes.shanks.units(neurons(neuron)).spikeTimes;
                for trial = 1:length(t_zero)
                    firingRate(trial,neuron) =sum(((queryT(trial)-0.05) < spikeTimes) & (spikeTimes < (queryT(trial)+0.05))) / 0.1;
                end
            end
            
            %Centre firing Rates
            firingRate = bsxfun(@minus,firingRate,mean(firingRate,1));
            
            %Fit model
            clear opts;
            opts.intr=1;
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
                
%                 %NG vs Go decoder
%                 opts.offset = offset_NGvG(cv.training(fold),:);
%                 fit1 = glmnet(trainX,trainY==3,'binomial',glmnetSet(opts));
%                 pNG = glmnetPredict(fit1,testX,lambda,'response',[],offset_NGvG(cv.test(fold),:)); 
%                 p = [1-pNG pNG];
%                 phat(cv.test(fold),2) = p(sub2ind(size(p), [1:length(testY)]', (testY==3)+1));
% 
% %                 keyboard;
% %                 %L vs R decoder
%                 trainIdx = cv.training(fold); 
%                 opts.offset = offset_LvR(trainIdx(Y==1 | Y==2));
%                 fit2 = glmnet(trainX(trainY==1 | trainY==2,:),trainY(trainY==1 | trainY==2),'binomial',glmnetSet(opts));
%                 
%                 testIdx = cv.test(fold);
%                 pR_G = glmnetPredict(fit2,testX(testY==1 | testY==2,:),lambda,'response',[],offset_LvR(testIdx(Y==1 | Y==2))); 
%                 p = [1-pR_G pR_G];
%                 phat(testIdx & Y<3,3) = p(sub2ind(size(p), [1:length(testY(testY==1 | testY==2))]', testY(testY==1 | testY==2)));
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
            title(dirs{session,2});
        else
            set(gca,'ytick','','ycolor','w'); ylabel('');
        end
        xlabel(n.epochs{e,1}); ylim([-0.2 0.8]); linkaxes(ax,'y');
    end
    set(gcf,'color','w'); 
    savefig(gcf,['\\basket.cortexlab.net\home\figures\GLM+NeuralActivity\glmnet_decoders\' dirs{session,2}])
    
end