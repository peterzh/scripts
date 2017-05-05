% %%%%%%%%%%%%% MODEL COMPARISON OVER MULTIPLE SESSIONS %%%%%%%%%%%%%%%%%%
saveDir = '\\basket.cortexlab.net\homes\peterzh\NeuroGLM-experimental\ModelFiles';

% expRefs = {'2015-05-28_1_Laveran',...
%     '2015-05-29_1_Laveran',...
%     '2015-05-30_1_Laveran',...
%     '2015-05-31_1_Laveran',...
%     '2015-06-01_1_Laveran',...
%     '2015-06-02_1_Laveran',...
%     '2015-06-03_1_Laveran',...
%     '2015-06-04_1_Laveran'};

% expRefs = {'2014-11-30_1_M140528NS1',...
%            '2014-12-01_1_M140528NS1',...
%            '2014-12-02_1_M140528NS1'};
% %
% expRefs = {'2015-01-24_1_M140617NS1',...
%            '2015-01-25_1_M140617NS1'};

% expRefs = {'2015-09-21_1_Hopkins','2015-09-21_2_Eijkman'};

% expRefs = {'2016-01-18_1_Morgan','2016-01-19_1_Morgan','2016-01-20_1_Morgan','2016-01-21_1_Morgan'};
% models = {'Offset','C','C-subset','C^N','C^N-subset','C50','C50-subset'};

% expRefs = {'2015-05-12_3_SS040','2015-05-18_5_SS040'};
expRefs = {'2015-05-21_3_SS040','2015-05-26_3_SS040','2015-05-12_3_SS040','2015-05-18_5_SS040'};
models = {'C^N-subset-2AFC','C^N-subset-Qlearning-noBias-2AFC','C^N-subset-Qlearning-2AFC'};


LLs = nan(length(expRefs),length(models));
for b = 1:length(expRefs)
    disp(['File: ' expRefs{b}]);
    %     g = GLM(expRefs{b});
    g = GLM(expRefs{b});
    for m = 1:length(models)
        g1 = g.setModel(models{m}).fitCV;
        LLs(b,m) = mean(-log2(g1.p_hat));
        save(fullfile(saveDir,[expRefs{b} '_' g1.modelString '_crossvalidated.mat']),'g1');
        disp('done');
    end
end

figure;
bar(LLs','stacked');
set(gca,'XTickLabel',models)
% for m = 1:length(models)
%     h(m)=subplot(1,length(models),m);
%     boxplot(cell2mat(LLs(:,m)),groupings,'plotstyle','compact','sym','r')
% end
% linkaxes(h,'y');


%% Comparing other kinds of models
expRefs = {'2016-02-22_1_Whipple','2016-02-25_1_Morgan','2016-02-19_1_Spemann'};
expRefs = dat.listExps('Spemann');
models = {'Offset','C50-subset','C50-subset-sharedS'};

LLs = nan(length(expRefs),length(models));
for b = 1:length(expRefs)
    disp(['File: ' expRefs{b}]);
    %     g = GLM(expRefs{b});
    try
        g = GLM(expRefs{b});
        if length(g.data.response)<150
            error('not enough data');
        end
        for m = 1:length(models)
            g1 = g.setModel(models{m}).fitCV(10);
            LLs(b,m) = mean(log2(g1.p_hat));
        end
    catch
    end
end
% LLs(isnan(LLs(:,1)),:)=[];

figure;
notBoxPlot(LLs);
% bar(LLs','stacked');
set(gca,'XTickLabel',models);

ylabel('log_2 likelihood');
xlabel('Model');

colormap(gray);