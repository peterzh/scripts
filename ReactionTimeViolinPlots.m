%% Violin plots showing the RT distribution for L and R choices, separated by stimulus. Combining sessions though.
% Pre-requisite: 
% http://uk.mathworks.com/matlabcentral/fileexchange/23661-violin-plots-for-plotting-multiple-distributions--distributionplot-m-
% subjects = {'Whipple','Morgan','Spemann','Hopkins','Eijkman'};
subjects = {'Whipple','Morgan','Spemann'};

% subjects = {'Chomsky'};
expRefs = cellfun(@(s)dat.listExps(s),subjects,'uni',0);

figure;
for s = 1:length(subjects)
    glms = cellfun(@(e)[GLM(e)],expRefs{s},'uni',0);
    
    params = dat.expParams(expRefs{s}); 
    
    responseWindow=[];
    for session = 1:length(params)
        try
            responseWindow(session) = params{session}.responseWindow;
        catch
            responseWindow(session) = NaN;
        end
    end
    glms(responseWindow~=1.5) = []; params(responseWindow~=1.5) = [];
    responseWindow = cellfun(@(st)(st.responseWindow),params);   
    performance = cellfun(@(g)mean(g.data.feedbackType==1),glms,'uni',0);
    taskDimensions = cellfun(@(g)(g.ContrastDimensions),glms);
    
%     many_contrasts = cellfun(@(g)size(unique(g.data.contrast_cond,'rows'),1)>8,glms,'uni',0);
    
    glms(cell2mat(performance)<0.6 | isnan(cell2mat(performance)) | taskDimensions==2)=[]; %remove discrim sessions

    structs = cellfun(@(g)[g.data],glms,'uni',0); %extract all data
    data = cell2mat(cellfun(@(st)[st.contrast_cond st.response st.repeatNum st.RT],structs,'uni',0)); %extract all data
    data(data(:,4)>1,:) = []; %exclude repeatNum > 1 data
    c = diff(data(:,[1,2]),[],2); 
%     tab = sortrows(tabulate(abs(c)),3); cVals = sort(tab(end-2:end,1)); %use the most common abs contrast values
%     cVals = sort([-cVals(2:end); cVals]);
%     tab = sortrows(tabulate(c),3); cVals = sort(tab(end-4:end,1)); %use the 6 most common contrast values
    
    cVals = [-1 -0.5 0 0.5 1];
    r = data(:,3);
    rt = data(:,5);
    
    rt_all = {};
    for response = 1:2
        for contrast = 1:length(cVals)
%             rt_i = rt(r == response & c == cVals(contrast));
            rt_all{response,contrast} = rt(r == response & c == cVals(contrast));
%             rt_all{response,contrast} = rt(r == response & abs(c-cVals(contrast))<0.05);
        end
    end
    
    subplot(length(subjects),3,3*s - 2);
    a=distributionPlot(rt_all(1,:),'histOpt',0,'widthDiv',[2 1],'histOri','left','color','g','showMM',0,'xValues',cVals,'globalNorm',2);
    b=distributionPlot(gca,rt_all(2,:),'histOpt',0,'widthDiv',[2 2],'histOri','right','color','b','showMM',0,'xValues',cVals,'globalNorm',2);
    
    if s == length(subjects)
        xlabel('Contrast'); ylabel('Reaction time [sec]');
        set(gca,'xtick',-1:0.25:1);
%         [LEGH,OBJH,OUTH,OUTM]=legend({'left choice','right choice'});
%         OBJH(4).FaceColor=[0 0 1];
%         legend boxoff;
    else
        set(gca,'xtick','','xcolor','w');
    end
    xlim([-1.2 1.2]);
    ylim([-0.1 2]); title(subjects{s});
    drawnow;
    set(gca,'box','off'); 
    
    
    cVals = [-1 -0.5 0 0.5 1];

    % Also plot median and quartiles
    
    hold on;
    
    Col = {'g','b'};
    for r = 1:2
        subplot(length(subjects),3,3*s -2 + r );
        hold on;
        for b=1:length(structs)
            D = structs{b};
            D = getrow(D,D.repeatNum==1);
            
            c = diff(D.contrast_cond,[],2);
            cVals = unique(c);
            rts = []; g = []; 
%             m = nan(length(cVals),1);
%             md = nan(size(m));
            
%             lq = nan(length(cVals),1); 
%             uq = nan(length(cVals),1);
            
            for contrast = 1:length(cVals)
                rt = D.RT(D.response == r & c==cVals(contrast));
                
                if length(rt)>5
%                 if ~isempty(rt)
                      m = median(rt);
                      md = mad(rt,1);
%                     m(contrast) = median(rt);
%                     md(contrast) = mad(rt);
%                     lq(contrast) = quantile(rt,.25);
%                     uq(contrast) = quantile(rt,.75);
                    
                    jtr_range = 0.03;
                    jtr= -0.5*jtr_range + jtr_range*rand;
                    line([cVals(contrast),cVals(contrast)]+jtr,[m-md,m+md]);
                    plot(cVals(contrast)+jtr,m,[Col{r} '.']);
%                     line([cVals(contrast),cVals(contrast)]+jtr,[lq(contrast),uq(contrast)]);
%                     plot(cVals(contrast)+jtr,m(contrast),[Col{r} '.']);
                else
%                     m(contrast) = NaN;
%                     lq(contrast) = NaN;
%                     uq(contrast) = NaN;
                    
                end
%                 rts = [rts;rt];
%                 g = [g; ones(length(rt),1)*contrast];
               
            end
            
            xlabel('Contrast');
            ylabel('Median RT +- MAD');
            
%             plot(cVals,m,[Col{r} '.'],'MarkerSize',10);
%             plot(cVals,lq,[Col{r} 'v'],'MarkerSize',5);
%             plot(cVals,uq,[Col{r} '^'],'MarkerSize',5);
%             boxplot(rts,g,'plotstyle','compact','positions',cVals(unique(g)),'labels',cVals(unique(g)),'symbol','','orientation','horizontal');
                     
        end
        xlim([-1 1]*1.1); ylim([0 1.5]);   
    end
    
end
set(gcf,'color','w');