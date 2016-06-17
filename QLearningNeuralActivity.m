%% Get nonDA w_init values
%This is to counter the tendency for DA blocks to give weird initial W
%results

nonDA_expRefs = {'2015-01-16_1_SS031_DA';
    '2015-01-19_1_SS031_DA';
    '2015-01-20_1_SS031_DA';
    '2015-01-21_1_SS031_DA';
    '2015-01-22_1_SS031_DA';
    '2015-01-23_3_SS031_DA';
    '2015-01-26_2_SS031_DA'};

%For each non-DA session, fit Qs and take the average the last 20 W values
%then take the average over each session and use that as the initial W for
%all DA sessions
w_end = nan(4,length(nonDA_expRefs));
for b = 1:length(nonDA_expRefs)
    q = Q(nonDA_expRefs(b)).setModel('aB').fit;

    %get last w values
    w = [q.wt{1};q.wt{2}];
    w_end(:,b)=mean(w(:,end-20:end),2); 
end
w_init_DA = mean(w_end,2); %[sL,bL,sR,bR]

%% Fit on some random DA sessions
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

for b = 1:length(DA_expRefs)
    q = Q(DA_expRefs(b)).setModel('aB');
    q.preset_winit = w_init_DA;
    q = q.fit;
end

%% Now go to neural data and use those values

%Load the file
% load('B:\stuff\Neurons_FR_DB_EJ008_Batch3.mat')
load('\\basket.cortexlab.net\home\stuff\Neurons_FR_DB_EJ008_Batch3_2.mat');
DATA = NeuronsFR_DB_EJ008DABatch3;
% [uq,ia,ic]=unique(arrayfun(@(d)size(d.data,1),DATA));

D=struct;
numNeurons = length(DATA);
        
neur = struct;
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
    
    q=Q(D).setModel('aB');
    q.preset_winit = w_init_DA;
    q=q.fit;

    
%     %remove outliers in QL and QR
%     QL(QL>(2.5*quantile(QL,0.75) - 1.5*quantile(QL,0.25)))=nan;
%     QR(QR>(2.5*quantile(QR,0.75) - 1.5*quantile(QR,0.25)))=nan;

    neur(n).alignBeep = DATA(n).data(:,19);
    neur(n).alignStim = DATA(n).data(:,23);
    neur(n).alignResp = DATA(n).data(:,28);
    neur(n).dQ = [q.QL-q.QR]';
    neur(n).behav = D;
        
    
end
close all;

%% plot full scatter for each neuron
neurEdges = -50:10:50;
dqEdges = -5:1:5;

for n = 1:numNeurons
   
    CAT.markerfill = {[.34 .26 .92],[.56 .51 .95],[.78 .75 .91],[1 1 1]*0.9,[.97 .75 .78],[.95 .51 .55],[.92 .26 .33]};
    CAT.markercolor = CAT.markerfill;
    CAT.markertype  = '.';
    CAT.markersize = 15;
    
    f=figure('name',['Neuron ' num2str(n) '. NumTrials ' num2str(numTrials)]); 
    subplot(2,3,1); scatterplot(neur(n).dQ,neur(n).alignBeep,'split',diff(neur(n).behav.stimulus,[],2),'CAT',CAT); ylabel('neur align Beep');
    subplot(2,3,2); scatterplot(neur(n).dQ,neur(n).alignStim,'split',diff(neur(n).behav.stimulus,[],2),'CAT',CAT);ylabel('neur align Stim');
    subplot(2,3,3); scatterplot(neur(n).dQ,neur(n).alignResp,'split',diff(neur(n).behav.stimulus,[],2),'CAT',CAT); xlabel('QL-QR');ylabel('neur align Resp');
    legend(cellfun(@(c)num2str(c),num2cell(sort(repmat(unique(diff(neur(n).behav.stimulus,[],2)),1,1))),'uni',0));
    
    dQ_discr = discretize(neur(n).dQ,dqEdges);
    dQ_discr(neur(n).behav.DA==1)=dQ_discr(neur(n).behav.DA==1)+0.07;
    
    smoothing=0.01;
    DA = neur(n).behav.DA;
    subplot(2,3,4);
    plot(dQ_discr(DA==0),neur(n).alignBeep(DA==0),'ko'); hold on;
    plot(fit(dQ_discr(DA==0),neur(n).alignBeep(DA==0),'smoothingspline','smoothingparam',smoothing),'k-');
    plot(dQ_discr(DA==1),neur(n).alignBeep(DA==1),'ro');
    plot(fit(dQ_discr(DA==1),neur(n).alignBeep(DA==1),'smoothingspline','smoothingparam',smoothing),'r-');
    legend off; set(gca,'XTickLabel',dqEdges,'XTick',1:length(dqEdges));ylabel('neur align Beep'); xlim([0 length(dqEdges)]);
    
    subplot(2,3,5);
    
    plot(dQ_discr(DA==0),neur(n).alignStim(DA==0),'ko'); hold on;
    plot(fit(dQ_discr(DA==0),neur(n).alignStim(DA==0),'smoothingspline','smoothingparam',smoothing),'k-');
    plot(dQ_discr(DA==1),neur(n).alignStim(DA==1),'ro');
    plot(fit(dQ_discr(DA==1),neur(n).alignStim(DA==1),'smoothingspline','smoothingparam',smoothing),'r-');
    legend off; set(gca,'XTickLabel',dqEdges,'XTick',1:length(dqEdges));ylabel('neur align Stim');xlim([0 length(dqEdges)]);
    
    subplot(2,3,6);
    
    plot(dQ_discr(DA==0),neur(n).alignResp(DA==0),'ko'); hold on;
    plot(fit(dQ_discr(DA==0),neur(n).alignResp(DA==0),'smoothingspline','smoothingparam',smoothing),'k-');
    plot(dQ_discr(DA==1),neur(n).alignResp(DA==1),'ro');
    plot(fit(dQ_discr(DA==1),neur(n).alignResp(DA==1),'smoothingspline','smoothingparam',smoothing),'r-');
    legend off; set(gca,'XTickLabel',dqEdges,'XTick',1:length(dqEdges));ylabel('neur align Resp');xlim([0 length(dqEdges)]);
    
    savefig(f,['\\basket.cortexlab.net\home\figures\QGLM\neuron' num2str(n) '.fig']);
    close all;
end

%% plot summary of corr
figure;
coeffs = [cellfun(@(c)c(2),arrayfun(@(st)(corrcoef(st.dQ,st.alignBeep)),neur,'uni',0));
    cellfun(@(c)c(2),arrayfun(@(st)(corrcoef(st.dQ,st.alignStim)),neur,'uni',0));
    cellfun(@(c)c(2),arrayfun(@(st)(corrcoef(st.dQ,st.alignResp)),neur,'uni',0))];
imagesc(coeffs');
ylabel('Neuron'); caxis([-1 1]);
cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
    linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
colormap(cmap);
set(gca,'XTickLabel',{'alignBeep','alignStim','alignResp'},'XTick',1:3);
colorbar;


%% Check PSTHs also show contrast dependence as armin finds in this data
stim=[];
neur=[];
neurAll=[];stimAll=[];
for n = 1:numNeurons
    subplot(6,6,n);
    correct=DATA(n).data(:,4)==1;
    stim=DATA(n).data(correct,2);
    neur=DATA(n).data(correct,23);
    m=pivottable(stim,[],neur,'median');
    plot(unique(stim),m,'.','markersize',10);
%     xlabel('contrast');
%     ylabel('neural activity');
    neurAll=[neurAll;neur];
    stimAll=[stimAll;stim];
end

m=pivottable(stimAll,[],neurAll,'median');
figure; plot(unique(stim),m,'.','markersize',40);