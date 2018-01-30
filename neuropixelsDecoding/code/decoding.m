preprocFiles = dir('../preproc/*.mat');

numT = 150;
epoch_dt = linspace(-1,+1,numT);

cols = [         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

for sess = 1:length(preprocFiles)
    load( fullfile(preprocFiles(sess).folder, preprocFiles(sess).name) );
    
    %Fit GLM to get parameters
    D = struct;
    D.stimulus = behav.stimulus;
    D.response = double(behav.response);
    g = GLM(D).setModel('C50-subset').fit;
    D.offset_ZL = g.ZL(g.parameterFits, g.Zinput(g.data));
    D.offset_ZR = g.ZR(g.parameterFits, g.Zinput(g.data));
    
    %Cross-validate C50-model to get baseline log likelihood
    g = GLM(D).setModel('C50-subset').fitCV(10);
    pL = g.p_hat(:,1);    pR = g.p_hat(:,2);    pNG = g.p_hat(:,3);
    bpt_baseline = struct;
    likelihood_Full = pL.*(g.data.response==1) + pR.*(g.data.response==2) + pNG.*(g.data.response==3);
    bpt_baseline.full = mean( log2( likelihood_Full ) );
    likelihood_GOvNOGO = (1 - pNG).*(g.data.response<3) + pNG.*(g.data.response==3);
    bpt_baseline.GOvNOGO = mean( log2( likelihood_GOvNOGO ) );
    
    likelihood_LvR = (pL./(1-pNG)).*(g.data.response==1) + (pR./(1-pNG)).*(g.data.response==2) + 10*(g.data.response==3);
    likelihood_LvR(likelihood_LvR==10) = [];
    bpt_baseline.LvR = mean( log2( likelihood_LvR ) );
    
    fig = figure('name',preprocFiles(sess).name,'units','normalized','outerposition',[0 0 1 1],'color','w');
    for e = 1:length(epoches)
        subplot(4,length(epoches), e); hold on; %Plot behavioural timestamps
        tstamps = behav.(epoches{e});
        hx=plot(zeros(size(tstamps)),1:length(tstamps),'.');
        hx.Color = cols(e,:);
        hx.MarkerSize=20;
        ylim([0 max(hx.YData)]);
        otherIdx = 1:length(epoches);
        otherIdx(e) = [];
        
        for ot = 1:length(otherIdx)
            dt = behav.(epoches{otherIdx(ot)}) - tstamps;
            hx = plot(dt,1:length(dt),'.');
            hx.Color = cols(otherIdx(ot),:);
            hx.YData=hx.YData;
            hx.MarkerSize=20;
        end
        xlabel(epoches{e},'Color',cols(e,:));
        xlim([-1 1]*1);
        set(gca, 'ytick','','ycolor','w');
        drawnow;
    end
        
    %     brainRegions = intersect(unique(population.region),...
    %                     {'Secondary motor area','Primary visual area'});
    brainRegions = unique(population.region);
        
    %Preallocate memory for bits per trial
    bpt = struct;
    bpt.full = nan(numT, length(brainRegions), length(epoches));
    bpt.GOvNOGO = nan(numT, length(brainRegions), length(epoches));
    bpt.LvR = nan(numT, length(brainRegions), length(epoches));

    numberOfRuns = length(epoches)*numT;
    tic; i = 1;
    
    traceAxs = gobjects(3,length(epoches));
    for ep = 1:length(epoches) %For each epoch
        fprintf('Epoch: %d/%d\n',ep,length(epoches));
        t_zero = behav.(epoches{ep});
        for t = 1:numT %for each timestep in that epoch
            fprintf('\tTime: %d/%d',t,numT);
            queryT = t_zero + epoch_dt(t);
            
            for region = 1:length(brainRegions) %For each region
                D.neur = interp1(population.Properties.UserData.fr_times, population.firingRates(region,:)', queryT);
                
                if any(isnan(D.neur))
                    keyboard;
                end
                
                g = GLM(D).setModel('neur').fitCV(10);
                pL = g.p_hat(:,1);    pR = g.p_hat(:,2);    pNG = g.p_hat(:,3);
                
                likelihood_Full = pL.*(g.data.response==1) + pR.*(g.data.response==2) + pNG.*(g.data.response==3);
                bpt.full(t,region,ep) = mean( log2( likelihood_Full ) );
                
                likelihood_GOvNOGO = (1 - pNG).*(g.data.response<3) + pNG.*(g.data.response==3);
                bpt.GOvNOGO(t,region,ep) = mean( log2( likelihood_GOvNOGO ) );
                
                likelihood_LvR = (pL./(1-pNG)).*(g.data.response==1) + (pR./(1-pNG)).*(g.data.response==2) + 10*(g.data.response==3);
                likelihood_LvR(likelihood_LvR==10) = [];
                bpt.LvR(t,region,ep) = mean( log2( likelihood_LvR ) );
            end
            
            %Also do some plotting
            traceAxs(1,ep) = subplot(4,length(epoches), ep + length(epoches) );
            plot(epoch_dt,bpt.full(:,:,ep) - bpt_baseline.full); drawnow;
            
            traceAxs(2,ep) = subplot(4,length(epoches), ep + 2*length(epoches) );
            plot(epoch_dt,bpt.GOvNOGO(:,:,ep) - bpt_baseline.GOvNOGO); drawnow;
            
            traceAxs(3,ep) = subplot(4,length(epoches), ep + 3*length(epoches) );
            plot(epoch_dt,bpt.LvR(:,:,ep) - bpt_baseline.LvR); drawnow;
            
            timePerRun = toc/i;
            timeLeft = timePerRun*( numberOfRuns-i );
            timeLeft = datestr(timeLeft/(60*60*24),'HH:MM:SS');
            fprintf('\t Time Left: %s\n',timeLeft);
            i = i+1;
            
        end
        
    end
    
    axs = get(gcf,'children');
    set(axs,'box','off');
    linkaxes(axs,'x');
    
    linkaxes(traceAxs,'y');
    
    ylabel(traceAxs(3,1),'L v R decoding');
    ylabel(traceAxs(2,1),'GO v NOGO decoding');
    ylabel(traceAxs(1,1),'Full decoding');
    
    
    legend(traceAxs(end),brainRegions);
    legend(traceAxs(end),'boxoff');
    set(fig,'PaperOrientation','landscape');

    
    print(fig,['../figures/' preprocFiles(sess).name(1:end-4)],'-dpdf','-bestfit');
    savefig(fig,['../figures/' preprocFiles(sess).name(1:end-4) '.fig']);
    
    save(['../decoding/' preprocFiles(sess).name ],'bpt_baseline','brainRegions','epoches','epoch_dt','bpt')
    
end
