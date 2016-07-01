%% Violin plots showing the RT distribution for L and R choices, separated by stimulus. Combining sessions though.
% Pre-requisite: 
% http://uk.mathworks.com/matlabcentral/fileexchange/23661-violin-plots-for-plotting-multiple-distributions--distributionplot-m-
subjects = {'Whipple','Morgan','Spemann','Hopkins','Eijkman'};
subjects = {'Chomsky'};
expRefs = cellfun(@(s)dat.listExps(s),subjects,'uni',0);

for s = 1:length(subjects)
    glms = cellfun(@(e)[GLM(e)],expRefs{s},'uni',0);
    
    params = dat.expParams(expRefs{s}); responseWindow = cellfun(@(st)(st.responseWindow),params);   
    performance = cellfun(@(g)mean(g.data.feedbackType==1),glms,'uni',0);
    
    glms(cell2mat(performance)<0.65 | isnan(cell2mat(performance)) | responseWindow~=1.5)=[]; %remove low performance sessions
    
    structs = cellfun(@(g)[g.data],glms,'uni',0); %extract all data
    data = cell2mat(cellfun(@(st)[st.contrast_cond st.response st.repeatNum st.RT],structs,'uni',0)); %extract all data
    data(min(data(:,[1,2]),[],2)>0,:) = [];    %exclude any discrimination trials
    data(data(:,4)>1,:) = []; %exclude repeatNum > 1 data
    c = diff(data(:,[1,2]),[],2); 
    tab = sortrows(tabulate(c),3); cVals = sort(tab(end-5:end,1)); %use the 6 most common contrast values
%     cVals = [-1 -0.5 0 0.5 1];
    r = data(:,3);
    rt = data(:,5);
    
    rt_all = {};
    for response = 1:2
        for contrast = 1:length(cVals)
            rt_i = rt(r == response & c == cVals(contrast));
            rt_all{response,contrast} = rt(r == response & c == cVals(contrast));
        end
    end
    
    subplot(length(subjects),1,s);
    distributionPlot(rt_all(1,:),'histOpt',0,'widthDiv',[2 1],'histOri','left','color','g','showMM',0,'xValues',cVals,'globalNorm',2);
    distributionPlot(gca,rt_all(2,:),'histOpt',0,'widthDiv',[2 2],'histOri','right','color','b','showMM',0,'xValues',cVals,'globalNorm',2);
    
    if s == length(subjects)
        xlabel('Contrast'); ylabel('Reaction time [sec]');
    else
%         set(gca,'xtick','','xcolor','w');
    end
    xlim([-1.2 1.2]);
    ylim([0 1.6]); title(subjects{s});
    drawnow;
    set(gca,'box','off'); 
end
set(gcf,'color','w');