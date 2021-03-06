%% Load widefield data
tic;
session =1;
[W,behav,behavT,otherInfo]=getWidefieldData('nick widefield',session);
type = 'widefield';
toc;


%% Generate figures for all neuropixels sessions
% for session = 35:48
for session = 2:48
    choiceDecodingFromNeuropixels(session);
    close all;
end

%% Extract traces from the figures and compile together to compare
figDir = '\\basket.cortexlab.net\home\decoding_files';
files = dir(figDir); files = {files.name}; files(1:2) = [];
isfig = cellfun(@(path) contains(path,'fig'), files);
files = files(isfig);

time = nan(length(files),50);
full = nan(size(time));
gng = nan(size(time));
lvr = nan(size(time));

region = cell(length(files),1);
for f = 1:length(files)
    fig = open(fullfile(figDir,files{f}));
    region{f} = fig.Name;
    
    axs = get(fig,'children');
    plt = get(axs(2),'children');
    
    time(f,:) = plt(1).XData;
    full(f,:) = plt(3).YData;
    gng(f,:) = plt(2).YData;
    lvr(f,:) = plt(1).YData;
    
    close(fig);
end

%Sort by maximum
[magnitude,peakidx] = max(lvr,[],2); 
[~,idx]=sort(magnitude);
% [~,idx]=sort(peakidx);
time = time(idx,:);
full = full(idx,:);
gng = gng(idx,:);
lvr = lvr(idx,:);
region = region(idx);


% figure;
% % subplot(2,1,1);
% plot(time',full');
figure;
imagesc(time(1,:),1:length(files),full);
caxis([min([full(:);gng(:);lvr(:)]) max([full(:);gng(:);lvr(:)])]);
set(gca,'YTickLabel',region,'ytick',1:length(region));
title('Full decoding'); grid on;

figure;
% subplot(1,3,2);
imagesc(time(1,:),1:length(files),gng);
caxis([min([full(:);gng(:);lvr(:)]) max([full(:);gng(:);lvr(:)])]);
set(gca,'YTickLabel',region,'ytick',1:length(region));

title('GO v NOGO decoding');  grid on;


figure;
imagesc(time(1,:),1:length(files),lvr);
caxis([min([full(:);gng(:);lvr(:)]) max([full(:);gng(:);lvr(:)])]);
set(gca,'YTickLabel',region,'ytick',1:length(region));

title('L v R decoding');   grid on;
colorbar;



%% OLD CODE
% %% Fit behav-only model to get Offsets
% tic;
% cfn = @(c,n,c50)((c.^n)./(c.^n + c50^n));
% guess_bpt = [];
% tab = tabulate(behav.response);
% tab = tab(:,3)/100;
% guess_bpt(1)=sum(tab.*log2(tab));
% 
% if max(behav.response)==3
%     
%     tab = tabulate(behav.response==3);
%     tab = cell2mat(tab(:,3))/100;
%     guess_bpt(2)=sum(tab.*log2(tab));
%     
%     tab = tabulate(behav.response(behav.response<3));
%     tab = tab(:,3)/100;
%     guess_bpt(3)=sum(tab.*log2(tab));
%     g = GLM(behav).setModel('C50-subset').fit;
%     
%     n = g.parameterFits(5);
%     c50 = g.parameterFits(6);
% else
%     g = GLM(behav).setModel('C50-subset-2AFC').fit;
%     
%     n = g.parameterFits(4);
%     c50 = g.parameterFits(5);
% end
% 
% %Construct design matrix
% X = sparse(cfn(behav.stimulus,n,c50));
% Y = behav.response;
% 
% csvwrite('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_Y.dat',Y); %Write input data to a file
% csvwrite('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_X.dat',full(X));
% csvwrite('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_Xtest.dat',full(X));
% system('"C:\Program Files\R\R-3.4.0\bin\R" CMD BATCH C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_interface.R');
% offset_MNR_NC50 = csvread('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_link.dat');
% toc;
% 
% %% Cross-validation: behav-only model baseline
% tic;
% numFolds = 5;
% cv = cvpartition(Y,'kfold',numFolds);
% phat_baseline = nan(cv.NumObservations,3);
% for fold = 1:cv.NumTestSets
%     trainX = X(cv.training(fold),:);
%     trainY = Y(cv.training(fold));
%     testX = X(cv.test(fold),:);
%     testY = Y(cv.test(fold));
%     
%     csvwrite('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_Y.dat',trainY); %Write input data to a file
%     csvwrite('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_X.dat',full(trainX));
%     csvwrite('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_Xtest.dat',full(testX));
%     system('"C:\Program Files\R\R-3.4.0\bin\R" CMD BATCH C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_interface.R');
%     p = csvread('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_phat.dat');
%     
%     %Orig MNR
%     phat_baseline(cv.test(fold),1) = p(:,1).*(testY==1) + p(:,2).*(testY==2) + p(:,3).*(testY==3);
%     
%     %Go vs NoGo decoder
%     pGO = 1-p(:,3);
%     phat_baseline(cv.test(fold),2) = (testY<3).*pGO + (testY==3).*(1-pGO);
%     
%     %L vs R decoder
%     pL_G = p(:,1)./pGO;
%     pR_G = p(:,2)./pGO;
%     obs = (testY==1).*pL_G + (testY==2).*pR_G; obs(obs==0)=nan;
%     phat_baseline(cv.test(fold),3) = obs;
% end
% baseline_bpt = nanmean(log2(phat_baseline));
% toc;

% %% Go through each epoch, extract neural data
% tic;
% numT = 150; %number of time bins to sample neural activity from within an epoch
% interval = [-1 1];
% tsteps = linspace(interval(1),interval(2),numT);
% 
% switch(type)
%     case 'spikes'
% %         numChannels = length(spikeTimes); %number of neurons
%         numChannels = numRegions; %number of population regions
%     case 'widefield'
%         numChannels = size(W.V,1); %number of components
% end
% numTrials = length(behav.response);
% numEpochs = length(behavT.timestamps);
% neuralActivity = nan(numTrials,numChannels,numT,numEpochs);
% 
% for e = 1:numEpochs %for each epoch
%     t_zero = behavT.(behavT.timestamps{e});
%     
%     for t = 1:length(tsteps) %for each timestep in that epoch
%         queryT = t_zero + tsteps(t);
%         
%         fprintf('%d/%d \n',t,length(tsteps));
%         
%         for channel = 1:numChannels %get activity for each channel
%             for trial = 1:numTrials
%                 
%                 switch(type)
%                     case 'spikes'
%                         neuralActivity(trial,channel,t,e) = sum(((queryT(trial)-0.05) < spikeTimesPop{channel}) & (spikeTimesPop{channel} < (queryT(trial)+0.05))) / 0.1;
%                     case 'widefield'
%                         idx = find(W.t>queryT(trial),1,'first');
%                         neuralActivity(trial,channel,t,e) = W.V(channel,idx); 
%                 end
%             end
%         end
%         
%     end
% end
% 
% %Remove mean activity across trials
% neuralActivity = bsxfun(@minus,neuralActivity,mean(neuralActivity,1));
% toc;
% 
% %% Per-epoch, per timepoint, predict choice from neural activity
% tic;
% bpt = nan(numEpochs,length(tsteps));
% for e = 1:numEpochs %for each epoch
%     disp(['Epoch ' num2str(e) '/' num2str(numEpochs)]);
%     for t = 1:length(tsteps)
%         disp(['Time ' num2str(t) '/' num2str(length(tsteps))]);
%         X = neuralActivity(:,:,t,e);
%         Y = behav.response;
%         
%         cv = cvpartition(Y,'kfold',numFolds);
%         phatNB_NC50 = nan(length(Y),1);
%         for fold = 1:numFolds
%             
%             trainX = X(cv.training(fold),:);
%             trainY = Y(cv.training(fold));
%             testX = X(cv.test(fold),:);
%             testY = Y(cv.test(fold));
%             trainoffset = offset_MNR_NC50(cv.training(fold),:);
%             testoffset = offset_MNR_NC50(cv.test(fold),:);
%             
%             csvwrite('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_offset.dat',trainoffset); %Write input data to a file
%             csvwrite('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_offsettest.dat',testoffset); %Write input data to a file
%             csvwrite('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_Y.dat',trainY);
%             csvwrite('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_X.dat',trainX);
%             csvwrite('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_Xtest.dat',testX);
%             system('"C:\Program Files\R\R-3.4.0\bin\R" CMD BATCH C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_interface_withoffset.R')
%             p = csvread('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_phat.dat');
%             
%             if max(behav.response)==2
%                 p = [1-p p];
%                 phatNB_NC50(cv.test(fold),1) = p(:,1).*(testY==1) + p(:,2).*(testY==2);
%             else
%                 phatNB_NC50(cv.test(fold),1) = p(:,1).*(testY==1) + p(:,2).*(testY==2) + p(:,3).*(testY==3);
%                 
% %                 %Go vs NoGo decoder
% %                 pGO = 1-p(:,3);
% %                 phatNB_NC50(cv.test(fold),2) = (testY<3).*pGO + (testY==3).*(1-pGO);
% %                 
% %                 %L VS R decoder
% %                 pL_G = p(:,1)./pGO;
% %                 pR_G = p(:,2)./pGO;
% %                 obs = (testY==1).*pL_G + (testY==2).*pR_G; obs(obs==0)=nan;
% %                 phatNB_NC50(cv.test(fold),3) = obs;
%             end
%         end
%         
%         bpt(e,t) = mean(log2(phatNB_NC50)) - baseline_bpt(1);
% 
%     end
% end
% toc;
% 
% %% Plot everything
% cols = [         0    0.4470    0.7410
%         0.8500    0.3250    0.0980
%         0.4940    0.1840    0.5560
%         0.4660    0.6740    0.1880
%         0.3010    0.7450    0.9330
%         0.6350    0.0780    0.1840];
% 
% f1=figure('name',otherInfo{1});
% for e = 1:length(behavT.timestamps) %for each epoch
%     subplot(2,length(behavT.timestamps),e); hold on;
%     
%     h=plot(zeros(length(behavT.(behavT.timestamps{e})),1),1:length(behavT.(behavT.timestamps{e})),'.');
%     h.Color = cols(e,:);
% %     h.YData=h.YData+75;
%     ylim([0 max(h.YData)]);
%     
%     otherIdx = 1:length(behavT.timestamps);
%     otherIdx(e) = [];
%     
%     for ot = 1:length(otherIdx)
%         
%         dt = behavT.(behavT.timestamps{otherIdx(ot)}) - behavT.(behavT.timestamps{e});
%         %             dt = sort(dt);
%         h = plot(dt,1:length(dt),'.');
%         h.Color = cols(otherIdx(ot),:);
%         h.YData=h.YData+75;
%         
%         h = histogram(dt,100);
%         h.FaceColor = cols(otherIdx(ot),:);
%         h.EdgeColor = cols(otherIdx(ot),:);
%         %             h.EdgeColor = [0 0 0];
%     end
%     
%     xlabel(behavT.timestamps{e},'Color',cols(e,:));
%     xlim([-1 1]*1);
%     
%     set(gca, 'ytick','','ycolor','w');
% end
% 
% h=[];
% for e = 1:numEpochs
%     h(e)=subplot(2,numEpochs,numEpochs+e);
%     plot(tsteps, bpt(e,:), 'LineWidth', 2);
%     
%     xlabel(behavT.timestamps{e},'Color',cols(e,:));
%     
%     if e == 1
%         ylabel('NEUR + NC50 log lik - NC50 BASELINE');
%         
%         if max(behav.response)==2
%             legend('L v R decoder');
%         else
%             legend('Full choice decoder');
%         end
%         legend boxoff;
%     else
%         set(gca,'ytick','','ycolor','w'); ylabel('');
%     end
%     
% end
% linkaxes(h,'xy');
% set(h,'box','off');
% 
% set(gcf,'color','w');
% 
% if max(behav.response)==3
%     dir = '2AUC';
% else
%     dir = '2AFC';
% end
% 
% makepretty;
% 
% savefig(gcf,['\\basket.cortexlab.net\home\figures\GLM+NeuralActivity\glmnet_decoders\' dir '\' otherInfo{1}]);


% Fit complete model to get full parameter set
X =  neuralActivity(:,:,1,1);
csvwrite('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_offset.dat',offset_MNR_NC50); %Write input data to a file
csvwrite('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_offsettest.dat',offset_MNR_NC50); %Write input data to a file
csvwrite('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_Y.dat',Y);
csvwrite('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_X.dat',X);
csvwrite('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_Xtest.dat',X);
system('"C:\Program Files\R\R-3.4.0\bin\R" CMD BATCH C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_interface_withoffset.R')
p = csvread('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_phat.dat');

neurParams = csvread('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_params.dat');
figure;
plot(aax,neurParams(:,1),neurParams(:,2),'ko'); title(aax,num2str(tsteps(t)));
drawnow;

%% Fit  temporal component of widefield to predict choice
tic;
figure('color','w','name',otherInfo{1});
for e = 1:numEpochs %for each epoch
    disp(['Epoch ' num2str(e) '/' num2str(numEpochs)]);
    subplot(2,numEpochs,e);
    
    numComponents = size(W.U,3);
    beta = nan(numComponents,2,length(tsteps));
    for t = 1:length(tsteps)
        X = [behav.stimulus neuralActivity(:,:,t,e)];
        Y = behav.response;
        csvwrite('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_Y.dat',Y);
        csvwrite('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_X.dat',X);
        csvwrite('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_Xtest.dat',X);
        system('"C:\Program Files\R\R-3.4.0\bin\R" CMD BATCH C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_interface.R')
        
        neurParams = csvread('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_params.dat');

%         b = mnrfit(neuralActivity(:,:,t,e),behav.response);
        beta(:,:,t) = neurParams(4:end,:);
        
        imagesc([svdFrameReconstruct(W.U,beta(:,1,t)) svdFrameReconstruct(W.U,beta(:,2,t))]);
        title(num2str(tsteps(t)));
%         caxis([-1 1]*1e-3);
        drawnow;
%         pause(0.05);
       
    end
    
    subplot(2,numEpochs,e+numEpochs);    
    title(behavT.timestamps{e});

    plot(tsteps,squeeze(beta(:,1,:)),'-',tsteps,squeeze(beta(:,2,:)),'--');
end
toc;




%% Create widefield movie
MOVIE = struct('cdata',[],'colormap',[]);

figure('color','w','name',otherInfo{1});
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
for e = 1:numEpochs
    subplot(2,numEpochs,e);   hold on;
    title(behavT.timestamps{e});

    plot(tsteps,squeeze(beta(:,1,:)),'-',tsteps,squeeze(beta(:,2,:)),'--');
    set(gca,'box','off');
    l(e) = line([0 0],get(gca,'ylim'));
    l(e).Color=[0 0 0];
end

for t = 1:length(tsteps)
    for e = 1:numEpochs

        subplot(2,numEpochs,e+numEpochs); 
        img = [svdFrameReconstruct(W.U,beta(:,1,t)) svdFrameReconstruct(W.U,beta(:,2,t))];
        imagesc([svdFrameReconstruct(W.U,beta(:,1,t)) svdFrameReconstruct(W.U,beta(:,2,t))]);
        set(gca,'xtick','','ytick','','box','off','xcolor','w','ycolor','w');
        caxis([-1 1]*1e-6); axis equal;  
        
        l(e).XData=[1 1]*tsteps(t);
%         title(num2str(tsteps(t)));
    end
    drawnow;
    MOVIE(t) = getframe(gcf);
end

%Write file
v = VideoWriter(['\\basket.cortexlab.net\home\figures\GLM+Widefield\' otherInfo{1} '.avi']);
open(v);
writeVideo(v,MOVIE);
close(v);





















