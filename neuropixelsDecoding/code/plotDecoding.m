function plotDecoding(eRef)
    load( fullfile('..','preproc',[eRef '.mat']) );
    load( fullfile('..','decoding',[eRef '.mat']) );
    
    disp(brainRegions);
%     return;
    f = figure('color','w','name',eRef);
    
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
        plot(epoch_dt, squeeze(bpt.full(:,e,:)) - bpt_baseline.full );
        set(gca,'xtick','','xcolor','w');
        if e > 1
            set(gca,'ytick','','ycolor','w');
        else
            ylabel('Full decoder');
        end
        
        dax(i) = subplot(4,length(epoches),e + 2*length(epoches)); i=i+1; %Plot decoding 
        plot(epoch_dt, squeeze(bpt.GOvNOGO(:,e,:)) - bpt_baseline.GOvNOGO );
        set(gca,'xtick','','xcolor','w');
        if e > 1
            set(gca,'ytick','','ycolor','w');
        else
            ylabel('GO v NOGO decoder');
        end
        
        dax(i) = subplot(4,length(epoches),e + 3*length(epoches)); i=i+1; %Plot decoding 
        plot(epoch_dt, squeeze(bpt.LvR(:,e,:)) - bpt_baseline.LvR );
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
    ylim([-0.05 0.3]);
    legend(brainRegions,'Orientation','vertical'); legend boxoff;
    
    set(get(gcf,'children'),'box','off');
    
    set(gcf,'PaperOrientation','landscape');
    
    set(gcf,'Position',[0 0 1000 600]);
    print(gcf,fullfile('..','figures',eRef),'-dpdf','-bestfit');
    
    savefig(gcf,fullfile('..','figures',eRef));
end