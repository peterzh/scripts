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