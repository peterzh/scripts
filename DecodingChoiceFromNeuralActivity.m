%%  

cfn = @(c,n,c50)((c.^n)./(c.^n + c50^n));
numFolds = 5;
numT = 81; %number of time bins to sample neural activity from within an epoch
% numT = 3;
interval = [-1 1];

for session = 2:(17*2)
    [spikeTimes,behav,behavT,o]=getSpikingData('nick neuropixels',session);

    cols = [         0    0.4470    0.7410
        0.8500    0.3250    0.0980
        0.4940    0.1840    0.5560
        0.4660    0.6740    0.1880
        0.3010    0.7450    0.9330
        0.6350    0.0780    0.1840];

    figure('name',o{1});
    for e = 1:length(behavT.timestamps) %for each epoch
        subplot(4,length(behavT.timestamps),e); hold on;
        
        h=plot(zeros(length(behavT.(behavT.timestamps{e})),1),1:length(behavT.(behavT.timestamps{e})),'.');
        h.Color = cols(e,:);
        h.YData=h.YData+75;
        ylim([0 max(h.YData)]);
        
        otherIdx = 1:length(behavT.timestamps);
        otherIdx(e) = [];
        
        for ot = 1:length(otherIdx)
            
            dt = behavT.(behavT.timestamps{otherIdx(ot)}) - behavT.(behavT.timestamps{e});
%             dt = sort(dt);
            h = plot(dt,1:length(dt),'.');
            h.Color = cols(otherIdx(ot),:);
            h.YData=h.YData+75;
            
            h = histogram(dt,100);
            h.FaceColor = cols(otherIdx(ot),:);
            h.EdgeColor = cols(otherIdx(ot),:);
%             h.EdgeColor = [0 0 0];
        end

        xlabel(behavT.timestamps{e},'Color',cols(e,:));
        xlim([-1 1]*1);
        
        set(gca, 'ytick','','ycolor','w');
    end
    
    
    guess_bpt = [];
    tab = tabulate(behav.response);
    tab = tab(:,3)/100;
    guess_bpt(1)=sum(tab.*log2(tab));
    
    if max(behav.response)==3
        
        tab = tabulate(behav.response==3);
        tab = cell2mat(tab(:,3))/100;
        guess_bpt(2)=sum(tab.*log2(tab));
        
        tab = tabulate(behav.response(behav.response<3));
        tab = tab(:,3)/100;
        guess_bpt(3)=sum(tab.*log2(tab));
        g = GLM(behav).setModel('C50-subset').fit;
        
        n = g.parameterFits(5);
        c50 = g.parameterFits(6);
    else
        g = GLM(behav).setModel('C50-subset-2AFC').fit;
        
        n = g.parameterFits(4);
        c50 = g.parameterFits(5);
    end

    %Construct design matrix
    X = sparse(cfn(behav.contrast_cond,n,c50));
    Y = behav.response;
    
    csvwrite('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_Y.dat',Y); %Write input data to a file
    csvwrite('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_X.dat',full(X));
    csvwrite('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_Xtest.dat',full(X));
    tic;
    system('"C:\Program Files\R\R-3.4.0\bin\R" CMD BATCH C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_interface.R');
    toc;
    offset_MNR_NC50 = csvread('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_link.dat');
    
    
    %Calculate baseline behavioural model log likelihood (crossval)
    cv = cvpartition(Y,'kfold',numFolds);
    phat_baseline = nan(cv.NumObservations,3);
    for fold = 1:cv.NumTestSets
        trainX = X(cv.training(fold),:);
        trainY = Y(cv.training(fold));
        testX = X(cv.test(fold),:);
        testY = Y(cv.test(fold));
        
        csvwrite('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_Y.dat',trainY); %Write input data to a file
        csvwrite('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_X.dat',full(trainX));
        csvwrite('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_Xtest.dat',full(testX));
        system('"C:\Program Files\R\R-3.4.0\bin\R" CMD BATCH C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_interface.R');
        p = csvread('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_phat.dat');

        %Orig MNR
        phat_baseline(cv.test(fold),1) = p(:,1).*(testY==1) + p(:,2).*(testY==2) + p(:,3).*(testY==3);
        
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
    
%     X = sparse(behav.contrast_cond);
%     Y = behav.response;
%     csvwrite('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_Y.dat',Y); %Write input data to a file
%     csvwrite('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_X.dat',full(X));
%     csvwrite('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_Xtest.dat',full(X));
%     system('"C:\Program Files\R\R-3.4.0\bin\R" CMD BATCH C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_interface.R');
%     offset_MNR_C = csvread('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_link.dat');
%     
    for e = 1:length(behavT.timestamps) %for each epoch
        t_zero = behavT.(behavT.timestamps{e});
        tsteps = linspace(interval(1),interval(2),numT);
        
        ax1(e)=subplot(3,length(behavT.timestamps),e+length(behavT.timestamps));
%         ax2(e)=subplot(4,length(behavT.timestamps),e+2*length(behavT.timestamps));
        ax3(e)=subplot(3,length(behavT.timestamps),e+2*length(behavT.timestamps));
        
        bptN=[];
        bptNB_C=[];
        bptNB_NC50=[];
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
            cv = cvpartition(size(firingRate,1),'kfold',numFolds);
            phatN = nan(cv.N,3);
%             phatNB_C = nan(cv.N,3);
            phatNB_NC50 = nan(cv.N,3);
            for fold = 1:cv.NumTestSets
                disp([num2str(fold) '/' num2str(cv.NumTestSets)]);
                trainX = firingRate(cv.training(fold),:);
                trainY = Y(cv.training(fold));
                testX = firingRate(cv.test(fold),:);
                testY = Y(cv.test(fold));
                
                %If trainY contains fewer than 2 entries for L/R/NG, then
                %skip this fold
                if any(sum(trainY==[1 2 3],1)<2)
                    warning('Skipping fold due to too few trials');
                else
                    
                    %NEURAL MODEL ONLY
                    csvwrite('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_Y.dat',trainY);
                    csvwrite('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_X.dat',full(trainX));
                    csvwrite('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_Xtest.dat',full(testX));
                    flag=system('"C:\Program Files\R\R-3.4.0\bin\R" CMD BATCH C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_interface.R');
                    p = csvread('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_phat.dat');
                    
                    
                    
                    if max(behav.response)==2
                        p = [1-p p];
                        phatN(cv.test(fold),1) = p(:,1).*(testY==1) + p(:,2).*(testY==2);
                        
                    else
                        phatN(cv.test(fold),1) = p(:,1).*(testY==1) + p(:,2).*(testY==2) + p(:,3).*(testY==3);
                        
                        %Go vs NoGo decoder
                        pGO = 1-p(:,3);
                        phatN(cv.test(fold),2) = (testY<3).*pGO + (testY==3).*(1-pGO);
                        
                        %L VS R decoder
                        pL_G = p(:,1)./pGO;
                        pR_G = p(:,2)./pGO;
                        obs = (testY==1).*pL_G + (testY==2).*pR_G; obs(obs==0)=nan;
                        phatN(cv.test(fold),3) = obs;
                    end
                    
                    %                 %NEURAL + C BEHAV MODEL
                    %
                    %                 trainoffset = offset_MNR_C(cv.training(fold),:);
                    %                 testoffset = offset_MNR_C(cv.test(fold),:);
                    %
                    %                 csvwrite('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_offset.dat',trainoffset); %Write input data to a file
                    %                 csvwrite('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_offsettest.dat',testoffset); %Write input data to a file
                    %                 csvwrite('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_Y.dat',trainY);
                    %                 csvwrite('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_X.dat',full(trainX));
                    %                 csvwrite('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_Xtest.dat',full(testX));
                    %                 system('"C:\Program Files\R\R-3.4.0\bin\R" CMD BATCH C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_interface_withoffset.R')
                    %                 p = csvread('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_phat.dat');
                    %
                    %                 if max(behav.response)==2
                    %                     p = [1-p p];
                    %                     phatNB_C(cv.test(fold),1) = p(:,1).*(testY==1) + p(:,2).*(testY==2);
                    %                 else
                    %                     phatNB_C(cv.test(fold),1) = p(:,1).*(testY==1) + p(:,2).*(testY==2) + p(:,3).*(testY==3);
                    %
                    %                     %Go vs NoGo decoder
                    %                     pGO = 1-p(:,3);
                    %                     phatNB_C(cv.test(fold),2) = (testY<3).*pGO + (testY==3).*(1-pGO);
                    %
                    %                     %L VS R decoder
                    %                     pL_G = p(:,1)./pGO;
                    %                     pR_G = p(:,2)./pGO;
                    %                     obs = (testY==1).*pL_G + (testY==2).*pR_G; obs(obs==0)=nan;
                    %                     phatNB_C(cv.test(fold),3) = obs;
                    %                 end
                    %
                    
                    %NEURAL + NC50 BEHAV MODEL
                    
                    trainoffset = offset_MNR_NC50(cv.training(fold),:);
                    testoffset = offset_MNR_NC50(cv.test(fold),:);
                    
                    csvwrite('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_offset.dat',trainoffset); %Write input data to a file
                    csvwrite('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_offsettest.dat',testoffset); %Write input data to a file
                    csvwrite('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_Y.dat',trainY);
                    csvwrite('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_X.dat',full(trainX));
                    csvwrite('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_Xtest.dat',full(testX));
                    system('"C:\Program Files\R\R-3.4.0\bin\R" CMD BATCH C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_interface_withoffset.R')
                    p = csvread('C:\Users\Peter\Desktop\GLMNET_DATA\glmnet_phat.dat');
                    
                    if max(behav.response)==2
                        p = [1-p p];
                        phatNB_NC50(cv.test(fold),1) = p(:,1).*(testY==1) + p(:,2).*(testY==2);
                    else
                        phatNB_NC50(cv.test(fold),1) = p(:,1).*(testY==1) + p(:,2).*(testY==2) + p(:,3).*(testY==3);
                        
                        %Go vs NoGo decoder
                        pGO = 1-p(:,3);
                        phatNB_NC50(cv.test(fold),2) = (testY<3).*pGO + (testY==3).*(1-pGO);
                        
                        %L VS R decoder
                        pL_G = p(:,1)./pGO;
                        pR_G = p(:,2)./pGO;
                        obs = (testY==1).*pL_G + (testY==2).*pR_G; obs(obs==0)=nan;
                        phatNB_NC50(cv.test(fold),3) = obs;
                    end
                end
            end
            
%             keyboard;
            
            bptN(t,:) = nansum(log2(phatN),1)/size(phatN,1) - guess_bpt;
            bptNB_NC50(t,:) = nansum(log2(phatNB_NC50),1)/size(phatNB_NC50,1) - baseline_bpt;
                        
            plot(ax1(e),tsteps(1:length(bptN)),bptN,'LineWidth',2); 
%             plot(ax2(e),tsteps(1:length(bptNB_C)),bptNB_C,'LineWidth',2); 
            plot(ax3(e),tsteps(1:length(bptNB_NC50)),bptNB_NC50,'LineWidth',2); 
            drawnow;
        end
%         set([ax1(e) ax2(e) ax3(e)],'box','off');
        set([ax1(e) ax3(e)],'box','off');
        
        if e == 1
            ylabel(ax1(e),'NEUR ONLY log lik - GUESS');
%             ylabel(ax2(e),'NEUR + C log lik');
            ylabel(ax3(e),'NEUR + NC50 log lik - NC50 BASELINE');
            title(ax1(e),o{1});
            
            if max(behav.response)==2
                legend('L v R decoder');
            else
                legend('Full choice decoder','G vs NG decoder','L v R decoder');
            end
            legend boxoff;
        else
%             set([ax1(e) ax2(e) ax3(e)],'ytick','','ycolor','w'); ylabel('');
            set([ax1(e) ax3(e)],'ytick','','ycolor','w'); ylabel('');
        end
        xlabel(ax3(e),behavT.timestamps{e},'Color',cols(e,:));
%         ylim([-0.2 0.8]);
    end
    
    linkaxes(ax1,'y');
%     linkaxes(ax2,'y');
    linkaxes(ax3,'y');
%     linkaxes(get(gcf,'children'),'x');
    set(gcf,'color','w'); 
    
    if max(behav.response)==3
        dir = '2AUC';
    else
        dir = '2AFC';
    end
    
    savefig(gcf,['\\basket.cortexlab.net\home\figures\GLM+NeuralActivity\glmnet_decoders\' dir '\' o{1}]);
    
    clear spikeTimes behav behavT;
end
