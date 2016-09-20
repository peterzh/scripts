
%% 2AFC simulation + bootstrap estimation
%Generate data

N = 1000;

D = table;
D.x = linspace(-1,1,N)';
% D.x = -1 + 2*rand(N,1);

biases = -4:2:4;
figure;
psy = subplot(3,1,3); hold(psy,'on');
cols = get(gca,'ColorOrder');

for bi = 1:length(biases)
    
    %Generate data with known bias and sensitivity
    s = 10;
    zr = biases(bi) + s*D.x;
    pyr = exp(zr)./(1+exp(zr));
    D.y = binornd(1,pyr);
    D.y = ~D.y + 1; %Recode response because mnrfit expects Y to be 1 and 2 not 0 and 1;
    
    plot(psy,sort(D.x),sort(pyr),'.-');
    
    %Fit to complete dataset
    param_est = mnrfit(D.x,D.y);
    
    %Bootstrap resample
    numIter = 200;
    B = nan(numIter,2);
    for i = 1:numIter
        d = datasample(D,height(D));
        B(i,:) = mnrfit(d.x,d.y);
    end
        
    subplot(3,length(biases),bi);
    ax=plot(B(:,1),B(:,2),'o');
    ax.Color = cols(bi,:);
%     title(['bias: ' num2str(biases(bi))]);
    set(gca,'box','off');
    hold on; 
    plot(biases(bi),s,'k+','markersize',10,'linewidth',1); 
    plot(param_est(1),param_est(2),'ko','markersize',10,'linewidth',1); 
    hold off;
    r = corr(B(:,1),B(:,2));
    title(num2str(r));
    
    
    if bi == 1
        xlabel('Bias as offset'); ylabel('Sens');
    end
    
    subplot(3,length(biases),length(biases)+bi);
    ax=plot(B(:,1)./B(:,2),B(:,2),'o');
    ax.Color = cols(bi,:);
%     title(['bias: ' num2str(biases(bi))]);
    set(gca,'box','off');
    hold on; 
    plot(biases(bi)/s,s,'k+','markersize',10,'linewidth',1); 
    plot(param_est(1)/param_est(2),param_est(2),'ko','markersize',10,'linewidth',1); 
    hold off;
    r = corr(B(:,1)./B(:,2),B(:,2));
    title(num2str(r));
    
    if bi == 1
        xlabel('Bias as contrast'); ylabel('Sens');
    end
    
    if bi == length(biases)
        legend('\beta_{bootstrap}','\beta_{true}','\beta_{est}')
    end
    
    drawnow;
end

xlabel(psy,'x');
ylabel(psy,'p(y==1)');
title(psy,[num2str(N) ' trials']);
set(gcf,'color','w');

%% repeat 2AFC simulation + bootstrap estimation to test whether true param falls within 95CI 95% of the time
%Generate data

numRepeats = 100;
numTrials = 500;

D = table;
D.contrast_cond = [linspace(1,0,numTrials/2)' zeros(numTrials/2,1);
        zeros(numTrials/2,1) linspace(0,1,numTrials/2)'];

bias = -4;
sens = 10;

INCI = zeros(numRepeats,2);

for i = 1:numRepeats
    disp(i);
    zl = bias + sens*(D.contrast_cond(:,1) - D.contrast_cond(:,2));
    pL = exp(zl)./(1+exp(zl));
    D.response = binornd(1,1-pL) + 1;

    
    numIter = 200;
    B = nan(numIter,2);
    for j = 1:numIter
        d = datasample(D,height(D));
        
%         %Fit with mnrfit
        B(j,:) = mnrfit((d.contrast_cond(:,1) - d.contrast_cond(:,2)),d.response );
%         
        %Fit with my code
%         g = GLM(table2struct(d,'ToScalar',true));
%         B(j,:) = g.setModel('AFC_diff').fit.parameterFits;
    end
    
    ci = quantile(B(:,1),[0.025 0.975]);
    if ci(1) < bias && bias < ci(2)
        INCI(i,1) = 1;
    end
    
    ci = quantile(B(:,2),[0.025 0.975]);
    if ci(1) < sens && sens < ci(2)
        INCI(i,2) = 1;
    end
end

mean(INCI)

%% 2AFC simulation + bayesian inference

N = 1000;

D = struct;
D.contrast_cond = [linspace(1,0,N/2)' zeros(N/2,1);
        zeros(N/2,1) linspace(0,1,N/2)'];

biases = -4:2:4;
figure;
psy = subplot(3,1,3); hold(psy,'on');
cols = get(gca,'ColorOrder');

for bi = 1:length(biases)
    
    %Generate data with known bias and sensitivity
    s = 10;
    zl = biases(bi) + s*(D.contrast_cond(:,1) - D.contrast_cond(:,2));
    pL = exp(zl)./(1+exp(zl));
    D.response = ~binornd(1,pL) + 1;
    
    plot(psy,diff(D.contrast_cond,[],2),1-pL,'.-');
    
    g = bGLM(D);
    g = g.setModel('AFC_diff');
    g.metHastIter = 10000;
    
    P = g.posterior;
    param_est = mean(P);
        
    subplot(3,length(biases),bi);
    ax=plot(P(:,1),P(:,2),'.');
    ax.Color = cols(bi,:);
%     title(['bias: ' num2str(biases(bi))]);
    set(gca,'box','off');
    hold on; 
    plot(biases(bi),s,'k+','markersize',10,'linewidth',1); 
    plot(param_est(1),param_est(2),'ko','markersize',10,'linewidth',1); 
    hold off;
    r = corr(P(:,1),P(:,2));
    title(num2str(r));
    
    if bi == 1
        xlabel('Bias as offset'); ylabel('Sens');
    end
    
    subplot(3,length(biases),length(biases)+bi);
    ax=plot(P(:,1)./P(:,2),P(:,2),'.');
    ax.Color = cols(bi,:);
%     title(['bias: ' num2str(biases(bi))]);
    set(gca,'box','off');
    hold on; 
    plot(biases(bi)/s,s,'k+','markersize',10,'linewidth',1); 
    plot(param_est(1)/param_est(2),param_est(2),'ko','markersize',10,'linewidth',1); 
    hold off;
    r = corr(P(:,1)./P(:,2),P(:,2));
    title(num2str(r));
    
    if bi == 1
        xlabel('Bias as contrast'); ylabel('Sens');
    end
    
    if bi == length(biases)
        legend('Posterior','\beta_{true}','\beta_{est}')
    end
    
    drawnow;
end

xlabel(psy,'x');
ylabel(psy,'pR');
title(psy,[num2str(N) ' trials']);
set(gcf,'color','w');

%% Generate 2AUC data
N = 10;

BL = -3;
BR = -3;
SL = 10;
SR = 10;

D = struct;
D.contrast_cond = [linspace(1,0,N/2)' zeros(N/2,1);
    zeros(N/2,1) linspace(0,1,N/2)'];

ZL = BL + SL*D.contrast_cond(:,1);
ZR = BR + SR*D.contrast_cond(:,2);
pL = exp(ZL)./(1+exp(ZL)+exp(ZR));
pR = exp(ZR)./(1+exp(ZL)+exp(ZR));
pNG = 1./(1+exp(ZL)+exp(ZR));

R = mnrnd(1,[pL pR pNG]);
[~,D.response]=max(R,[],2);

g = bGLM(D); g.metHastIter = 10000;
g = g.setModel('C-subset');
post = g.fit;
disp(mean(post));

%% 2AUC simulations on DETECTION task

N = 500;

biases = -4:2:4;

figure;
psyL = subplot(5,3,13); hold(psyL,'on');
psyNG = subplot(5,3,14); hold(psyNG,'on');
psyR = subplot(5,3,15); hold(psyR,'on');
cols = get(gca,'ColorOrder');

for bi = 1:length(biases)
    
    %Generate data
    BL = -biases(bi);
    BR = biases(bi);
    SL = 10;
    SR = 10;
    
    D = table;
    D.contrast_cond = [linspace(1,0,N/2)' zeros(N/2,1);
        zeros(N/2,1) linspace(0,1,N/2)'];
       
    ZL = BL + SL*D.contrast_cond(:,1);
    ZR = BR + SR*D.contrast_cond(:,2);
    pL = exp(ZL)./(1+exp(ZL)+exp(ZR));
    pR = exp(ZR)./(1+exp(ZL)+exp(ZR));
    pNG = 1./(1+exp(ZL)+exp(ZR));
    
    R = mnrnd(1,[pL pR pNG]);
    [~,D.response]=max(R,[],2);
    
%     D.contrast_cond = bsxfun(@minus,D.contrast_cond,mean(D.contrast_cond));

    %Fit to entire dataset
    g = GLM(table2struct(D,'ToScalar',true));
    param_est = g.setModel('C-subset').fit.parameterFits;
    
    plot(psyL,diff(D.contrast_cond,[],2),pL,'.-');ylim([0 1]);
    plot(psyNG,diff(D.contrast_cond,[],2),pNG,'.-');ylim([0 1]);
    plot(psyR,diff(D.contrast_cond,[],2),pR,'.-');ylim([0 1]);

    %Bootstrap resample
    numIter = 100;
    B = nan(numIter,4);
    for i = 1:numIter
        d = datasample(D,height(D));
        g = GLM(table2struct(d,'ToScalar',true));
        g = g.setModel('C-subset').fit;
        
        B(i,:) = g.parameterFits;
    end
    
    % Plot parameter estimates against each other
    subplot(5,length(biases),bi);
    ax=plot(B(:,1),B(:,2),'o');
    ax.Color = cols(bi,:);
    set(gca,'box','off');
    hold on;
    plot(BL,SL,'k+','markersize',10,'linewidth',1);
    plot(param_est(1),param_est(2),'ko','markersize',10,'linewidth',1);
    hold off;
    r = corr(B(:,1),B(:,2));
    title(num2str(r));
    
    if bi == 1
        xlabel('Bias_L as offset'); ylabel('Sens_L');
    end
    
    
    subplot(5,length(biases),length(biases)+bi);
    ax=plot(B(:,3),B(:,4),'o');
    ax.Color = cols(bi,:);
    set(gca,'box','off');
    hold on;
    plot(BR,SR,'k+','markersize',10,'linewidth',1);
    plot(param_est(3),param_est(4),'ko','markersize',10,'linewidth',1);    
    hold off;
    r = corr(B(:,3),B(:,4));
    title(num2str(r));
    
    if bi == 1
        xlabel('Bias_R as offset'); ylabel('Sens_R');
    end
    
    
    subplot(5,length(biases),2*length(biases)+bi);
    ax=plot(B(:,1)./B(:,2),B(:,2),'o');
    ax.Color = cols(bi,:);
    set(gca,'box','off');
    hold on;
    plot(BL/SL,SL,'k+','markersize',10,'linewidth',1);
    plot(param_est(1)/param_est(2),param_est(2),'ko','markersize',10,'linewidth',1);
    hold off;
    r = corr(B(:,1)./B(:,2),B(:,2));
    title(num2str(r));
    
    if bi == 1
        xlabel('Bias_L as contrast'); ylabel('Sens_L');
    end
    
    
    subplot(5,length(biases),3*length(biases)+bi);
    ax=plot(B(:,3)./B(:,4),B(:,4),'o');
    ax.Color = cols(bi,:);
    set(gca,'box','off');
    hold on;
    plot(BR/SR,SR,'k+','markersize',10,'linewidth',1);
    plot(param_est(3)/param_est(4),param_est(4),'ko','markersize',10,'linewidth',1);
    hold off;
    r = corr(B(:,3)./B(:,4),B(:,4));
    title(num2str(r));
    
    if bi == 1
        xlabel('Bias_R as contrast'); ylabel('Sens_R');
    end
    
    if bi == length(biases)
        legend('\beta_{bootstrap}','\beta_{true}','\beta_{est}')
    end
    
   drawnow;
end

ylabel(psyL,'pL');
ylabel(psyNG,'pNG');
ylabel(psyR,'pR');
set(gcf,'color','w');

%% 2AUC simulations on DETECTION task BAYESIAN INFERENCE

N = 100;

biases = -4:2:4;

f=figure;
psyL = subplot(3,3,7); hold(psyL,'on');
psyNG = subplot(3,3,8); hold(psyNG,'on');
psyR = subplot(3,3,9); hold(psyR,'on');
cols = get(gca,'ColorOrder');

for bi = 1:length(biases)
    
    %Generate data
    BL = -biases(bi);
    BR = biases(bi);
    SL = 10;
    SR = 10;
    
    D = struct;
    D.contrast_cond = [linspace(1,0,N/2)' zeros(N/2,1);
        zeros(N/2,1) linspace(0,1,N/2)'];
       
    ZL = BL + SL*D.contrast_cond(:,1);
    ZR = BR + SR*D.contrast_cond(:,2);
    pL = exp(ZL)./(1+exp(ZL)+exp(ZR));
    pR = exp(ZR)./(1+exp(ZL)+exp(ZR));
    pNG = 1./(1+exp(ZL)+exp(ZR));
    
    R = mnrnd(1,[pL pR pNG]);
    [~,D.response]=max(R,[],2);
    
    %Fit to entire dataset
    g = bGLM(D);
    g = g.setModel('C-subset');
    g.metHastIter = 10000;
        
    figure(f);
    
    plot(psyL,diff(D.contrast_cond,[],2),pL,'.-'); ylim([0 1]);
    plot(psyNG,diff(D.contrast_cond,[],2),pNG,'.-');ylim([0 1]);
    plot(psyR,diff(D.contrast_cond,[],2),pR,'.-');ylim([0 1]);

    %Bootstrap resample
    B = g.posterior;
    
    % Plot parameter estimates against each other
    subplot(3,length(biases),bi);
    imagesc(corrcoef(B)); caxis([-1 1]);
    set(gca,'xtick',1:size(B,2),'ytick',1:size(B,2),'xticklabels',g.parameterLabels,'yticklabels',g.parameterLabels,'XTickLabelRotation',90);
   drawnow;
   
   figure; gplotmatrix(B);
   
end

ylabel(psyL,'pL');
ylabel(psyNG,'pNG');
ylabel(psyR,'pR');
set(gcf,'color','w');

cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
                    linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
            colormap(cmap);
