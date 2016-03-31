%Load the file
% load('B:\stuff\Neurons_FR_DB_EJ008_Batch3.mat')
load('B:\stuff\Neurons_FR_DB_EJ008_Batch3_2.mat');
DATA = NeuronsFR_DB_EJ008DABatch3;
[unique,ia,ic]=unique(arrayfun(@(d)size(d.data,1),DATA));

D=struct;
numNeurons = length(DATA);

for n = 1:numNeurons
    numTrials = size(DATA(n).data,1);
    stim = DATA(n).data(:,2);
    stim = sort([stim zeros(numTrials,1)],2);
    stim = abs(stim);
    D.stimulus = stim;
    
    D.action = (DATA(n).data(:,3)==1)+1;
    
%     D.feedbackType = DATA(n).data(:,4);
    
    D.reward = DATA(n).data(:,4)==1;
    D.DA = DATA(n).data(:,5);
    
    q=Q(D).fit;
    w_init = q.fitWINIT(q.alpha,q.gamma);
    p.alpha=q.alpha;
    p.beta=q.beta;
    p.gamma=q.gamma;
    p.sL_init = w_init(1);
    p.bL_init = w_init(2);
    p.sR_init = w_init(3);
    p.bR_init = w_init(4);
    
    wt = q.calculateWT(p);
    xt = q.x;
    
    ph = q.calculatePHAT(wt,p.beta);
    bpt = q.calculateLOGLIK(ph)/length(q.data.action);
    loglik=bpt - q.guess_bpt;
    
    QL=[]; QR=[];
    for t = 1:numTrials
        QL(t) = wt{1}(:,t)'*xt{1}(:,t);
        QR(t) = wt{2}(:,t)'*xt{2}(:,t);
    end
    QL=QL'; QR=QR';
    
    %remove outliers in QL and QR
    QL(QL>(2.5*quantile(QL,0.75) - 1.5*quantile(QL,0.25)))=nan;
    QR(QR>(2.5*quantile(QR,0.75) - 1.5*quantile(QR,0.25)))=nan;
    
    neur.alignBeep = DATA(n).data(:,19);
    neur.alignStim = DATA(n).data(:,23);
    neur.alignResp = DATA(n).data(:,28);
   
    CAT.markerfill = {[.34 .26 .92],[.56 .51 .95],[.78 .75 .91],[1 1 1]*0.9,[.97 .75 .78],[.95 .51 .55],[.92 .26 .33]};
    CAT.markercolor = CAT.markerfill;
    CAT.markertype  = '.';
    CAT.markersize = 15;
    
    f=figure('name',['Neuron ' num2str(n) '. NumTrials ' num2str(numTrials)]); 
    subplot(3,3,1); scatterplot(QL,neur.alignBeep,'split',diff(D.stimulus,[],2),'CAT',CAT); ylabel('neur align Beep');
    subplot(3,3,2); scatterplot(QR,neur.alignBeep,'split',diff(D.stimulus,[],2),'CAT',CAT);
    subplot(3,3,3); scatterplot(QL-QR,neur.alignBeep,'split',diff(D.stimulus,[],2),'CAT',CAT);
    
    subplot(3,3,4); scatterplot(QL,neur.alignStim,'split',diff(D.stimulus,[],2),'CAT',CAT); ylabel('neur align Stim');
    subplot(3,3,5); scatterplot(QR,neur.alignStim,'split',diff(D.stimulus,[],2),'CAT',CAT);
    subplot(3,3,6); scatterplot(QL-QR,neur.alignStim,'split',diff(D.stimulus,[],2),'CAT',CAT);
    
    subplot(3,3,7); scatterplot(QL,neur.alignResp,'split',diff(D.stimulus,[],2),'CAT',CAT); ylabel('neur align Resp'); xlabel('QL');
    subplot(3,3,8); scatterplot(QR,neur.alignResp,'split',diff(D.stimulus,[],2),'CAT',CAT); xlabel('QR');
    subplot(3,3,9); scatterplot(QL-QR,neur.alignResp,'split',diff(D.stimulus,[],2),'CAT',CAT); xlabel('QL-QR');

    legend(cellfun(@(c)num2str(c),num2cell(sort(repmat(unique(diff(D.stimulus,[],2)),1,1))),'uni',0));
    
    savefig(f,['B:\figures\GLM+Qlearning\correlatingwithNeuralActivity\neuron' num2str(n) '_numTrials' num2str(numTrials) '_QlogLik' num2str(loglik) '.fig']);
    close all;
end