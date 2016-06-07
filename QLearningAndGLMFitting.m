%% Many sessions concatenated 
models = {'nolearn','aB+aS','aB','aS'};

expRefs = {
    '2015-01-14_1_SS031_DA'
'2015-01-14_2_SS031_DA'
'2015-01-14_3_SS031_DA'
'2015-01-14_4_SS031_DA'
'2015-01-15_1_SS031_DA'
'2015-01-15_2_SS031_DA'
'2015-01-15_3_SS031_DA'
'2015-01-15_4_SS031_DA'
'2015-01-16_1_SS031_DA'
'2015-01-16_2_SS031_DA'
'2015-01-16_4_SS031_DA'
'2015-01-16_5_SS031_DA'
'2015-01-19_1_SS031_DA'
'2015-01-19_2_SS031_DA'
'2015-01-19_3_SS031_DA'
'2015-01-19_4_SS031_DA'
'2015-01-20_1_SS031_DA'
'2015-01-20_2_SS031_DA'
'2015-01-20_3_SS031_DA'
'2015-01-20_4_SS031_DA'
'2015-01-21_1_SS031_DA'
'2015-01-21_2_SS031_DA'
'2015-01-21_3_SS031_DA'
'2015-01-21_4_SS031_DA'
'2015-01-22_1_SS031_DA'
'2015-01-22_2_SS031_DA'
'2015-01-22_3_SS031_DA'
'2015-01-22_4_SS031_DA'
'2015-01-23_1_SS031_DA'
'2015-01-23_2_SS031_DA'
'2015-01-23_3_SS031_DA'
'2015-01-26_1_SS031_DA'
'2015-01-26_2_SS031_DA'
'2015-01-26_3_SS031_DA'
};

% q = Q(expRefs).fit;

fig_dir = 'B:\figures\GLM+Qlearning'; 
for b = 1:length(expRefs)
    q=Q(expRefs(b));
%     q=q.setModel('aB+aS').fit;
%     set(gcf, 'Position', get(0,'Screensize'));
% %     print(fullfile(fig_dir,[expRefs{b} '.pdf' ]),'-dpdf','-painters');
%     savefig(fullfile(fig_dir,[expRefs{b} '.fig' ]));
    close all;
    
    p=[];
    for m=1:length(models)
        p(:,m)=q.setModel(models{m}).crossvalidate;
    end
    
    f=figure;
    bar(mean(log2(p))); 
    set(gca,'XTickLabel',models);
    ylabel('Log_2 likelihood');
    title(expRefs{b},'interpreter','none');
    savefig(f,fullfile(fig_dir,[expRefs{b} 'modelComparison.fig' ]));
    
end

%% Fit QGLM to non-DA blocks and use their last W values as the initial values for DA blocks
%This is to counter the tendency for DA blocks to give weird initial W
%results

nonDA_expRefs = {'2015-01-16_1_SS031_DA';
    '2015-01-19_1_SS031_DA';
    '2015-01-20_1_SS031_DA';
    '2015-01-21_1_SS031_DA';
    '2015-01-22_1_SS031_DA';
    '2015-01-23_3_SS031_DA';
    '2015-01-26_2_SS031_DA'};

DA_expRefs = {
    '2015-01-14_1_SS031_DA';
    '2015-01-14_2_SS031_DA';
    '2015-01-14_3_SS031_DA';
    '2015-01-14_4_SS031_DA';
    '2015-01-15_1_SS031_DA';
    '2015-01-15_2_SS031_DA';
    '2015-01-15_3_SS031_DA';
    '2015-01-15_4_SS031_DA';
    '2015-01-16_2_SS031_DA';
    '2015-01-16_4_SS031_DA';
    '2015-01-16_5_SS031_DA';
    '2015-01-19_2_SS031_DA';
    '2015-01-19_3_SS031_DA';
    '2015-01-19_4_SS031_DA';
    '2015-01-20_2_SS031_DA';
    '2015-01-20_3_SS031_DA';
    '2015-01-20_4_SS031_DA';
    '2015-01-21_1_SS031_DA';
    '2015-01-21_2_SS031_DA';
    '2015-01-21_3_SS031_DA';
    '2015-01-21_4_SS031_DA';
    '2015-01-22_2_SS031_DA';
    '2015-01-22_3_SS031_DA';
    '2015-01-22_4_SS031_DA';
    '2015-01-23_1_SS031_DA';
    '2015-01-23_2_SS031_DA';
    '2015-01-26_1_SS031_DA';
    '2015-01-26_3_SS031_DA';
    };

%For each non-DA session, fit Qs and take the average the last 20 W values
%then take the average over each session and use that as the initial W for
%all DA sessions
w_end = nan(4,length(nonDA_expRefs));
for b = 1:length(nonDA_expRefs)
    q = Q(nonDA_expRefs(b)).setModel('aB').fit;

    %get last w values
    ax=get(gcf,'children');
    w = cell2mat(get(get(ax(end),'children'),'YData'));
    close gcf;
    w_end(:,b)=mean(w(:,end-20:end),2);
end
w_init_DA = flipud(mean(w_end,2));

for b = 1:length(DA_expRefs)
    q = Q(DA_expRefs(b)).setModel('aB');
    q.preset_winit = w_init_DA;
    q.fit;
ax=get(gcf,'children');
end


