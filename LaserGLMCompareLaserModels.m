%% %%%%%%%%%%%%% C^N MODELS INCORPORATING LASER REGRESSORS %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expRef = '2015-11-19_1_Eijkman';
g = laserGLM(expRef);

models = {'C^N-subset-nolaser','C^N-subset-laser','C^N-subset-laser-offset','C^N-subset-laser-scale'};
model_labels = {'NoLaser','Full','Offset','Scale'};
for m = 1:length(models)
    g = g.setModel(models{m});
    for site = 1:size(g.inactivationSite,1)
        g1 = g.fitCV(site);
        save(fullfile(saveDir,[g1.expRef '_laser-' num2str(site) '_model-'  g1.modelString '.mat']),'g1');
    end
end

%% Plot
p = [];
NLL=[];
for m = 1:length(models)
    for site = 1:6
        try
            load(fullfile(saveDir,[expRef '_laser-' num2str(site) '_model-'  models{m} '.mat']));
            NLL(site,m) = -sum(log2(g1.p_hat));
        catch
            warning('site does not exist');
        end
    end
end
disp(NLL)
bar(NLL);
set(gca,'XTickLabel',model_labels);
ylabel('-SUM log_2 likelihood');