preprocFiles = dir('../preproc/*.mat');
decodingFiles = dir('../decoding/*.mat');

ROIs = {'Primary visual area','Secondary motor area'}; %Get averages for these areas only
relbpt.full = nan(80, length(ROIs), 4, length(preprocFiles));
relbpt.gonogo = nan(80, length(ROIs), 4, length(preprocFiles));
relbpt.lr = nan(80, length(ROIs), 4, length(preprocFiles));
                
for sess = 1:length(preprocFiles)
    load( fullfile(decodingFiles(sess).folder,decodingFiles(sess).name) );
    load( fullfile(preprocFiles(sess).folder,preprocFiles(sess).name) );
    epoches = {'stimulusOnTime', 'firstMoveTime'};

    
    f = figure('color','w','name',decodingFiles(sess).name(1:end-4));
    %
    cols = [         0    0.4470    0.7410
        0.8500    0.3250    0.0980
        0.4940    0.1840    0.5560
        0.4660    0.6740    0.1880
        0.3010    0.7450    0.9330
        0.6350    0.0780    0.1840];
    
    dax = [];
    i=1;
    for e = 1:length(epoches)
        subplot(4,length(epoches),e); %Plot behavioural histogram
        hold on;
        
        hx=plot(zeros(length( behav.(epoches{e}) ),1),1:length( behav.(epoches{e}) ),'.');
        hx.Color = cols(e,:);
        hx.MarkerSize=20;
        ylim([0 max(hx.YData)]);
        otherIdx = 1:length(epoches);
        otherIdx(e) = [];
        
        for ot = 1:length(otherIdx)
            dt = behav.(epoches{otherIdx(ot)}) - behav.(epoches{e});
            hx = plot(dt,1:length(dt),'.');
            hx.Color = cols(otherIdx(ot),:);
            hx.YData=hx.YData;
            hx.MarkerSize=20;
            
            %             hx = histogram(dt,100);
            %             hx.FaceColor = cols(otherIdx(ot),:);
            %             hx.EdgeColor = cols(otherIdx(ot),:);
            %             h.EdgeColor = [0 0 0];
        end
        
        xlabel(epoches{e},'Color',cols(e,:));
        xlim([-1 1]*1);
        
        set(gca, 'ytick','','ycolor','w');
        
        dax(i) = subplot(4,length(epoches),e + length(epoches)); i=i+1; %Plot decoding
        hx=plot(epoch_dt, squeeze(bpt.full(:,:,e)) - bpt_baseline.full );
        
        if length(hx) >= 8
            set(hx(8:end),'LineStyle','--');
        end
        
        set(gca,'xtick','','xcolor','w');
        if e > 1
            set(gca,'ytick','','ycolor','w');
        else
            ylabel('Full decoder');
        end
        
        dax(i) = subplot(4,length(epoches),e + 2*length(epoches)); i=i+1; %Plot decoding
        hx=plot(epoch_dt, squeeze(bpt.GOvNOGO(:,:,e)) - bpt_baseline.GOvNOGO );
        if length(hx) >= 8
            set(hx(8:end),'LineStyle','--');
        end
        set(gca,'xtick','','xcolor','w');
        if e > 1
            set(gca,'ytick','','ycolor','w');
        else
            ylabel('GO v NOGO decoder');
        end
        
        dax(i) = subplot(4,length(epoches),e + 3*length(epoches)); i=i+1; %Plot decoding
        hx=plot(epoch_dt, squeeze(bpt.LvR(:,:,e)) - bpt_baseline.LvR );
        if length(hx) >= 8
            set(hx(8:end),'LineStyle','--');
        end
        if e > 1
            set(gca,'ytick','','ycolor','w');
        else
            ylabel('L v R decoder');
            xlabel('Time relative to epoch');
        end
        
        if e > 1
            set(gca,'ytick','','ycolor','w');
        end
    end
    
    linkaxes(dax,'xy');
    ylim([-0.05 0.5]);
    legend(brainRegions,'Orientation','vertical'); legend boxoff;
    
    set(get(gcf,'children'),'box','off');
    
    set(gcf,'PaperOrientation','landscape');
    
    set(gcf,'Position',[0 0 1000 600]);
    print(gcf,fullfile('..','figures',decodingFiles(sess).name(1:end-4)),'-dpdf','-bestfit');
    
    savefig(gcf,fullfile('..','figures',decodingFiles(sess).name(1:end-4)));
    
    
    %Only get data for areas of interest
    [areas,idx,~] = intersect(brainRegions,ROIs);
    if ~isempty(areas)
        for epoch = 1:length(epoches)
            for a = 1:length(areas)
                
                relbpt.full(:,a,epoch,sess) = bpt.full(:,idx(a),epoch) - bpt_baseline.full;
                relbpt.gonogo(:,a,epoch,sess) = bpt.GOvNOGO(:,idx(a),epoch) - bpt_baseline.GOvNOGO;
                relbpt.lr(:,a,epoch,sess) = bpt.LvR(:,idx(a),epoch) - bpt_baseline.LvR;
            end
            
        end
    end
    
end

%% Compute average over sessions for these ROIs
avgbpt = structfun(@(f) nanmean(f,4), relbpt, 'uni', 0);

figure('color','w');

for epoch = 1:length(epoches)
    ax_full = subplot(3,length(epoches),epoch); hold on;
    ax_gonogo = subplot(3,length(epoches),epoch+length(epoches)); hold on;
    ax_lr = subplot(3,length(epoches),epoch+2*length(epoches)); hold on;
    for area = 1:length(ROIs)
        plt(1) = plot(ax_full,epoch_dt, avgbpt.full(:,area,epoch), '-');
        plt(2) = plot(ax_gonogo,epoch_dt, avgbpt.gonogo(:,area,epoch), '-');
        plt(3) = plot(ax_lr,epoch_dt, avgbpt.lr(:,area,epoch), '-');
        
        set(plt,'Color',cols(area,:),'Linewidth',1);
    end
    
    if epoch==1
        ylabel(ax_full,'full decoder');
        ylabel(ax_gonogo,'GO v NOGO decoder');
        ylabel(ax_lr,'L v R decoder');
    end
    
    
        xlabel(ax_lr,epoches{epoch});
    
end
linkaxes(get(gcf,'children'),'xy');
legend(ROIs);
legend boxoff; ylim([0 0.3]);

%% Plot GO v NOGO decoding for all regions 

decodingFiles = dir('../decoding/*.mat');

FULL=[];
GNG=[];
LR=[];
regions = {};
for sess = 1:length(decodingFiles)
    load( fullfile(decodingFiles(sess).folder,decodingFiles(sess).name) );
    
    b = bpt.full(:,:,strcmp(epoches,'stimulusOnTime')) - bpt_baseline.full;
    FULL = [FULL; b'];
    
    b = bpt.GOvNOGO(:,:,strcmp(epoches,'stimulusOnTime')) - bpt_baseline.GOvNOGO;
    GNG = [GNG; b'];
    
    b = bpt.LvR(:,:,strcmp(epoches,'stimulusOnTime')) - bpt_baseline.LvR;
    LR = [LR; b'];
    
    regions = [regions; brainRegions];
end

%Sort by region name
[~,idx]=sort(regions);

%Sort by size of decoding
[~,idx]=sort(max(FULL,[],2),'descend');

%Sort by onset time
[~,thres_idx]=max(FULL>0.1,[],2);
[~,idx]=sort(thres_idx,'descend');

figure('color','w');
ha = tight_subplot(1,3,0.02,0.05,[0.15 0.01]); 
for i = 1:3; hold(ha(i),'on'); end;

imagesc(ha(1),epoch_dt,[],FULL(idx,:));
line(ha(1),[0 0],get(ha(1),'ylim'),'LineWidth',1','Color',[0 0 0],'LineStyle','--');
set(ha(1),'YTickLabel',regions(idx),'ytick',1:length(idx));
xlabel(ha(1),'Time relative to first movement');

imagesc(ha(2),epoch_dt,[],GNG(idx,:));
line(ha(2),[0 0],get(ha(2),'ylim'),'LineWidth',1','Color',[0 0 0],'LineStyle','--'); 

imagesc(ha(3),epoch_dt,[],LR(idx,:));
line(ha(3),[0 0],get(ha(3),'ylim'),'LineWidth',1','Color',[0 0 0],'LineStyle','--');
set(gca,'YTickLabel','');
set(ha,'xlim',[-1 1]*0.5,'ylim',[0 length(regions)+1],'clim',[0 0.15], 'xticklabelmode','auto');
title(ha(1),'Full decoder');
title(ha(2),'Go v NoGo decoder');
title(ha(3),'Left v Right decoder');
