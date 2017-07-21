
for s = 1:100
    
    session =s;
    [W,behav,behavT,otherInfo]=getWidefieldData('nick widefield',session);
    type = 'widefield';
    
    % Go through each epoch, extract neural data
    
    numT = 150; %number of time bins to sample neural activity from within an epoch
    interval = [-1 1];
    tsteps = linspace(interval(1),interval(2),numT);
    
    switch(type)
        case 'spikes'
            numChannels = length(spikeTimes); %number of neurons
        case 'widefield'
            numChannels = size(W.V,1); %number of components
    end
    numTrials = length(behav.response);
    numEpochs = length(behavT.timestamps);
    neuralActivity = nan(numTrials,numChannels,numT,numEpochs);
    
    for e = 1:numEpochs %for each epoch
        t_zero = behavT.(behavT.timestamps{e});
        
        for t = 1:length(tsteps) %for each timestep in that epoch
            queryT = t_zero + tsteps(t);
            
            for channel = 1:numChannels %get activity for each channel
                for trial = 1:numTrials
                    
                    switch(type)
                        case 'spikes'
                            neuralActivity(trial,channel,t,e) = sum(((queryT(trial)-0.05) < spikeTimes{channel}) & (spikeTimes{channel} < (queryT(trial)+0.05))) / 0.1;
                        case 'widefield'
                            idx = find(W.t>queryT(trial),1,'first');
                            neuralActivity(trial,channel,t,e) = W.V(channel,idx);
                    end
                end
            end
            
        end
    end
    
    %Remove mean activity across trials
    neuralActivity = bsxfun(@minus,neuralActivity,mean(neuralActivity,1));
    
    
    
    
    tic;
    figure('color','w','name',otherInfo{1});
    for e = 1:numEpochs %for each epoch
        disp(['Epoch ' num2str(e) '/' num2str(numEpochs)]);
        subplot(2,numEpochs,e);
        
        numComponents = size(W.U,3);
        beta = nan(numComponents,2,length(tsteps));
        for t = 1:length(tsteps)
%             X = [behav.stimulus neuralActivity(:,:,t,e)];
            X = [ neuralActivity(:,:,t,e) ];
            Y = behav.response;
            csvwrite('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_Y.dat',Y);
            csvwrite('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_X.dat',X);
            csvwrite('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_Xtest.dat',X);
            system('"C:\Program Files\R\R-3.4.0\bin\R" CMD BATCH C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_interface.R')
            
            neurParams = csvread('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_params.dat');
            
            %         b = mnrfit(neuralActivity(:,:,t,e),behav.response);
            beta(:,:,t) = neurParams(2:end,:);
            
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
    v = VideoWriter(['\\basket.cortexlab.net\home\figures\GLM+Widefield\' otherInfo{1} '_neurOnly.avi']);
    open(v);
    writeVideo(v,MOVIE);
    close(v);
    
end