

%% %%%%%%%%%%%%% TESTING WITH SIMULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expRef = '2015-10-30_1_Hopkins'; %Used for getting input contrasts from only
% expRef = '2015-05-28_1_Laveran';
trueModel = 'C50-subset(N=1)';
trueModelParams = [-1 3 -2 6.1 0.3];

g = simulateGLM(expRef,trueModel,trueModelParams,1);
g.parameterBounds = [-inf -20 -inf -20 0;
                     +inf 20 +inf 20 20];

simexpRef = ['Simulation_' trueModel];
g.expRef = simexpRef;
g.parameterFits = g.trueParameters;
g.plotFit;

%% Run simulation sampling + fit many times to see whether parameter fits are settling into the true value.
% g.parameterStart = @()(2*randn(1,length(g.parameterLabels)));
reSampling = 1500;
params = nan(reSampling,length(g.parameterLabels));
g.parameterStart = @()(2*randn(1,length(g.parameterLabels)));

for i = 1:reSampling
    disp(i);
    g = g.resample.fit;
    params(i,:) = g.parameterFits;
end

figure;
for p = 1:length(g.parameterLabels)
    subplot(1,length(g.parameterLabels),p);
    hist(params(:,p));
    hold on;
    line([trueModelParams(p) trueModelParams(p)], get(gca,'ylim'));
    hold off;
    title(g.parameterLabels{p});
end

%% Run CV and save results
models = {'Offset','CL+CR-subset','C^N-subset','C^NL^NR-subset','C50-subset','C50-subset(N=1)'};
modelLabels = {'Offset','CL+CR sub','C^N sub','C^{NL,NR} sub','C50 sub','c50 sub n=1'};

for m = 1:length(models)
    g = g.setModel(models{m});
    g = g.fitCV;
    save(fullfile(saveDir,[g.expRef '_' g.modelString '_crossvalidated.mat']),'g');
    disp('done');
end

%% Plot
BitsperTrial = [];
BaselineEntropy = [];
ModelID = [];
for m = 1:length(models)
    p_hat=[];
    load(fullfile(saveDir,['Simulation_' trueModel '_' models{m} '_crossvalidated.mat']));
    p_hat = g.p_hat;
    nanIdx = find(isnan(p_hat));
    p_hat(nanIdx)=[]; %remove failed predictions
    
    bpt = -sum(log2(p_hat))/length(p_hat);
    BitsperTrial(m,1) = bpt;
    
    %Calculate baseline, if the model was guessing based only on the
    %fact that the mouse has a certain total ratio of L:R:NG
    tab = tabulate(g.data.response);
    tab = tab(:,3)/100;
    BaselineEntropy(m,1) = -sum(tab.*log2(tab));
    
end
figure;
bar(BitsperTrial,'grouped','EdgeColor','none');
ylabel('Bits per trial');
set(gca,'XTickLabel',modelLabels);
%ylim([0.9 1]);