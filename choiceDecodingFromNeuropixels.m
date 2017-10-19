function choiceDecodingFromNeuropixels(session)
numT = 50; %number of time bins to sample neural activity from within an epoch
interval = [-0.3 0.3];
tsteps = linspace(interval(1),interval(2),numT);

[spikeTimes,behav,behavT,otherInfo]=getSpikingData('neuropixels',session);

%Trim away spike Times beyond the time period of the experiment
spikeTimes = cellfun(@(st) st( (behavT.stimOn(1)-2) < st & st < behavT.feedbackTime(end)+2), spikeTimes, 'uni', 0);

%Overwrite relevant timestamps
behavT.timestamps = {'stimOn','firstMove'};

%Group spikes into cortical regions
numRegions = size(otherInfo{2},1);
spikeTimesPop = cell(numRegions,1);
cluDepths = otherInfo{1};

ratePop = cell(numRegions,1);
t_bin = 0.1;
rate_t = (behavT.stimOn(1)-2):t_bin:(behavT.feedbackTime(end)+2);
for region = 1:numRegions
    UB = otherInfo{2}.upperBorder(region);
    LB = otherInfo{2}.lowerBorder(region);
    idx = LB < cluDepths & cluDepths < UB;
    spikeTimesPop{region} = sort(cat(1,spikeTimes{idx}));
    ratePop{region} = hist(spikeTimesPop{region}, rate_t)/t_bin;
end
clear spikeTimes; clear spikeTimesPop;


% Behav-only model
if max(behav.response)==2
    warning('not coded for 2AFC');
    return;
end
% 
% gBase = GLM(behav).setModel('Offset').fitCV(10); ph = gBase.p_hat;
% phat_full = ph(:,1).*(behav.response==1) + ph(:,2).*(behav.response==2) + ph(:,3).*(behav.response==3);
% phat_GvNG = (1-ph(:,3)).*(behav.response<3) + (ph(:,3)).*(behav.response==3);
% phat_LvR = (ph(:,1)./(1-ph(:,3))).*(behav.response==1) + (ph(:,2)./(1-ph(:,3))).*(behav.response==2); phat_LvR(phat_LvR==0) = [];
% 
% baseline_Base(1) =  mean( log2( phat_full ) );
% baseline_Base(2) =  mean( log2( phat_GvNG ) );
% baseline_Base(3) =  mean( log2( phat_LvR ) );


gC50 = GLM(behav).setModel('C50-subset').fitCV(10); ph = gC50.p_hat;
phat_full = ph(:,1).*(behav.response==1) + ph(:,2).*(behav.response==2) + ph(:,3).*(behav.response==3);
phat_GvNG = (1-ph(:,3)).*(behav.response<3) + (ph(:,3)).*(behav.response==3);
phat_LvR = (ph(:,1)./(1-ph(:,3))).*(behav.response==1) + (ph(:,2)./(1-ph(:,3))).*(behav.response==2); phat_LvR(phat_LvR==0) = [];

baseline_C50(1) =  mean( log2( phat_full ) );
baseline_C50(2) =  mean( log2( phat_GvNG ) );
baseline_C50(3) =  mean( log2( phat_LvR ) );

gC50Neur = GLM(behav).setModel('C50-subset-neur');
gNeurOnly = GLM(behav).setModel('neur');

%Add regularisation for neural regressor
gC50Neur.regularise = @(p) 0.001*sum(p(7:8).^2);

thisSubj = otherInfo{3};
thisDate = otherInfo{4};
regionNames = otherInfo{2}.name;
numEpochs = length(behavT.timestamps);

cols = [         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

for region = 1:numRegions
    figure('color','w','name',[thisSubj ' ' thisDate ' ' regionNames{region}]);
    
    hxFull={}; hxGNG={};
    for epoch = 1:numEpochs %for each epoch
        subplot(2,numEpochs,epoch); hold on; %Plot behavioural timestamps
        hx=plot(zeros(length(behavT.(behavT.timestamps{epoch})),1),1:length(behavT.(behavT.timestamps{epoch})),'.');
        hx.Color = cols(epoch,:);
        hx.MarkerSize=20;
        ylim([0 max(hx.YData)]);
        otherIdx = 1:numEpochs;
        otherIdx(epoch) = [];
        
        for ot = 1:length(otherIdx)
            dt = behavT.(behavT.timestamps{otherIdx(ot)}) - behavT.(behavT.timestamps{epoch});
            hx = plot(dt,1:length(dt),'.');
            hx.Color = cols(otherIdx(ot),:);
            hx.YData=hx.YData;
            hx.MarkerSize=20;
            
            hx = histogram(dt,100);
            hx.FaceColor = cols(otherIdx(ot),:);
            hx.EdgeColor = cols(otherIdx(ot),:);
            %             h.EdgeColor = [0 0 0];
        end

        xlabel(behavT.timestamps{epoch},'Color',cols(epoch,:));
        xlim([-1 1]*1);
        
        set(gca, 'ytick','','ycolor','w');
        
        
        hxFull{epoch,1}=subplot(2,numEpochs,epoch + numEpochs); grid on;
%         hxGNG{epoch,1}=subplot(4,numEpochs,epoch + 2*numEpochs); grid on;
%         hxLR{epoch,1}=subplot(4,numEpochs,epoch + 3*numEpochs); grid on;
        t_zero = behavT.(behavT.timestamps{epoch});
        
        bptC50 = [];
        bptNeurOnly=[];
        for t = 1:numT %for each timestep in that epoch
            queryT = t_zero + tsteps(t);
            
            %             tic;
            %             wr = WithinRanges(spikeTimesPop{region}, queryT + [-0.05 0.05], [1:length(queryT)]' );
            %             toc;
            %             estFiringRate = sum(wr,1)'/0.1;
            estFiringRate = interp1(rate_t,ratePop{region},queryT);
            estFiringRate(isnan(estFiringRate))=0;
            
            gC50Neur.data.neur = estFiringRate;
            gC50Neur = gC50Neur.fitCV(10); ph = gC50Neur.p_hat;
            phat_full = ph(:,1).*(behav.response==1) + ph(:,2).*(behav.response==2) + ph(:,3).*(behav.response==3);
            phat_GvNG = (1-ph(:,3)).*(behav.response<3) + (ph(:,3)).*(behav.response==3);
            phat_LvR = (ph(:,1)./(1-ph(:,3))).*(behav.response==1) + (ph(:,2)./(1-ph(:,3))).*(behav.response==2); phat_LvR(phat_LvR==0) = [];
            
            bptC50(t,1) =  mean( log2( phat_full ) );
            bptC50(t,2) =  mean( log2( phat_GvNG ) );
            bptC50(t,3) =  mean( log2( phat_LvR ) );
            
            
%             gNeurOnly.data.neur = estFiringRate;
%             gNeurOnly = gNeurOnly.fitCV(10); ph = gNeurOnly.p_hat;
%             phat_full = ph(:,1).*(behav.response==1) + ph(:,2).*(behav.response==2) + ph(:,3).*(behav.response==3);
%             phat_GvNG = (1-ph(:,3)).*(behav.response<3) + (ph(:,3)).*(behav.response==3);
%             phat_LvR = (ph(:,1)./(1-ph(:,3))).*(behav.response==1) + (ph(:,2)./(1-ph(:,3))).*(behav.response==2); phat_LvR(phat_LvR==0) = [];
%             
%             bptNeurOnly(t,1) =  mean( log2( phat_full ) );
%             bptNeurOnly(t,2) =  mean( log2( phat_GvNG ) );
%             bptNeurOnly(t,3) =  mean( log2( phat_LvR ) );
%             
            
      plot(hxFull{epoch},tsteps(1:t), bptC50(:,1) - baseline_C50(1), 'k-',...
                tsteps(1:t), bptC50(:,2) - baseline_C50(2), 'r-',...
                tsteps(1:t), bptC50(:,3) - baseline_C50(3), 'b-',...
                'LineWidth',1);
%             
%             plot(hxFull{epoch},tsteps(1:t), bptC50(:,1) - baseline_C50(1), 'k-',...
%                 tsteps(1:t), bptNeurOnly(:,1) - baseline_Base(1), 'k:',...
%                 'LineWidth',1);
            
%             plot(hxGNG{epoch},tsteps(1:t), bptC50(:,2) - baseline_C50(2), 'r-',...
%                 tsteps(1:t), bptNeurOnly(:,2) - baseline_Base(2), 'r:',...
%                 'LineWidth',1);
%             
%             plot(hxLR{epoch},tsteps(1:t), bptC50(:,3) - baseline_C50(3), 'b-',...
%                 tsteps(1:t), bptNeurOnly(:,3) - baseline_Base(3), 'b:',...
%                 'LineWidth',1);
            
            
            %             plot(hxFull{epoch},tsteps(1:t), bpt(:,1) - baseline_bpt(1), 'k-',...
            %                 tsteps(1:t), bpt(:,2) - baseline_bpt(2), 'r-',...
            %                 tsteps(1:t), bpt(:,3) - baseline_bpt(3), 'b-',...
            %                 'LineWidth',2);
            
            
            
            %             plot(hxGNG{epoch},tsteps(1:t), bptNeur(:,1) - g.guess_bpt, 'k-',...
            %                 tsteps(1:t), bptNeur(:,2) - baseline_bpt(2), 'r-',...
            %                 tsteps(1:t), bptNeur(:,3) - baseline_bpt(3), 'b-',...
            %                 'LineWidth',2);
            %              plot(hxGNG{epoch},tsteps(1:t), bptNeur(:,1) - g.guess_bpt, 'k-','LineWidth',2);
            
            drawnow;
            
        end
        %         title(behavT.timestamps{epoch});
        
        %         if epoch == 1
        %             ylabel('loglik relative to m_0');
        %         end
        
    end
    linkaxes(get(gcf,'children'),'x');
%     linkaxes([hxFull{:} hxGNG{:} hxLR{:}],'y');
    linkaxes([hxFull{:}],'y');
    
    set(get(gcf,'children'),'box','off');
    
    
%     
%     set([hxLR{2:end} hxGNG{2:end} hxFull{2:end} ],'ytick','','ycolor','w');
%     set([hxFull{:} hxGNG{:}],'xtick','','xcolor','w');
    
        set([hxFull{2:end} ],'ytick','','ycolor','w');
%     set([hxFull{:} ],'xtick','','xcolor','w');
    
    
    
    hleg = legend(hxFull{end},'full choice','go v nogo','L v R');
%     htitle = get(hleg,'Title'); set(htitle,'String','Full choice decoder');
%     htitle.Color = 'k'; hleg.Box='off';
    
%     hleg = legend(hxGNG{end},'neur + stim','neur only');
%     htitle = get(hleg,'Title'); set(htitle,'String','GO vs NOGO decoder');
%     htitle.Color = 'r'; hleg.Box='off';
%     
%     hleg = legend(hxLR{end},'neur + stim','neur only');
%     htitle = get(hleg,'Title'); set(htitle,'String','L vs R decoder');
%     htitle.Color = 'b'; hleg.Box='off';
%     
    
    filename = get(gcf,'name');
    savefig(gcf,['\\basket.cortexlab.net\home\decoding_files\' filename '.fig'])
end
end