%% 
close all; clear all;

%% Set model configs/parameter
visContrast = [1 0; 0 1; 1 1; 0 0];
visLabels = {'CL','CR','CL=CR','C=[0 0]'};
respLabels = {'L','R','NG'};
r_labels = {'LV1' 'RV1', 'LS1' 'RS1' 'LM2' 'RM2'}';

p = struct;
p.weights = [-1 1 -1 -1 -1 1;
             1 -1 -1 -1 1 -1];
       
p.baselineActivity = 0;
p.InactivationScalingFactor = 0.01;

%% Simulate nonlaser and laser data
pedestals = [0 0.1 0.24 0.54];
% figure('color','w');
%psych curve: 2D representation
cols = [0 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.4940    0.1840    0.5560];
% for ped = 1:length(pedestals)
%     subplot(length(pedestals),1,ped); hold on;
%     set(gca,'colororder',cols);
% 
%     %Plot predictions
%     testCont = [linspace(1-pedestals(ped),0,100)' zeros(100,1); zeros(100,1) linspace(0,1-pedestals(ped),100)'] + pedestals(ped);
%     p_hat = nan(size(testCont,1),3);
%     for t=1:size(testCont,1)
%         p_hat(t,:) = simulateMechanisticFcn(testCont(t,:),p,0);
%     end
%     plot(diff(testCont,[],2),p_hat,'linewidth',0.5);
%     xlim([-1 1]); ylim([0 1]);
%     
%     
%     title(['pedestal= ' num2str(pedestals(ped))]);
%     axis square;
% end

figure('color','w','name','unilateral');
cmap = [ linspace(0,1,100)' linspace(0,1,100)' ones(100,1);
                ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];

cVals = [0 0.1 0.24 0.54];            
p_pred = nan(length(cVals),length(cVals),3);
for cl = 1:length(cVals)
    for cr = 1:length(cVals)
        p_pred(cl,cr,:) = simulateMechanisticFcn([cVals(cl) cVals(cr)],p,0);
    end
end

titles = {'p( Left | c)','p( Right | c)','p( NoGo | c)'};
for i=1:3
    h(i)=subplot(1,3,i);
    %                         p_plot = prop(:,:,i);
    %                         p_plot = [p_plot; nan(1,length(uniqueCR))];
    %                         p_plot = [p_plot, nan(length(uniqueCL)+1,1)];
    %                         pcolor([uniqueCR; 1],[uniqueCL; 1],p_plot); caxis([0 1]); shading('flat');
    %                         set(gca,'box','off');
    %                         set(gca,'xtick',0:0.2:1,'ytick',0:0.2:1);
    imagesc(p_pred(:,:,i)); caxis([0 1]);
    set(gca,'xtick',1:length(cVals),'ytick',1:length(cVals),'xticklabels',cVals,'yticklabels',cVals);
    set(gca,'YDir','normal','box','off');
    %                         xlabel('Contrast right');
    %                         ylabel('Contrast left');
    title(titles{i});
    axis square;
    %                         set(gca,'XTick','','YTick',0:0.1:0.5);
    if i > 1
        set(gca,'XTick',[],'ytick',[]);
    end
    
    if i == 1
        %                                 xlabel('Contrast right');
        ylabel('Contrast left');
    end
end


            
figure('color','w');           
for vc = 1:size(visContrast,1)
    
    p_nL = simulateMechanisticFcn(visContrast(vc,:),p,0);
    
    dP = [];
    for area = 1:6
        p_L = simulateMechanisticFcn(visContrast(vc,:),p,area,0);
        dP(area,:) = p_L - p_nL;
    end
    
    
    %Plot
    for r=1:3
        subplot(4,3,3*(vc-1)+ r);
        imagesc( flipud(reshape(dP(:,r),2,3)') );
        caxis([-1 1]*0.8);
        axis equal;
        xlim([0.5 2.5]);
        
        set(gca,'box','off','xtick','','xcolor','w','ytick','','ycolor','w');
        
        if r == 1
            set(gca,'ycolor','k'); ylabel(visLabels{vc});
        end
        
        if vc == size(visContrast,1)
            set(gca,'xcolor','k'); xlabel(respLabels{r});
        end
    end
    
end
colormap(cmap);

figure('color','w','name','bilateral');
for vc = 1:size(visContrast,1)

    p_nL = simulateMechanisticFcn(visContrast(vc,:),p,0);

    dP = [];
    i=1;
    for area = [2 4 6]
        p_L = simulateMechanisticFcn(visContrast(vc,:),p,area,1);
        dP(i,:) = p_L - p_nL;
        i = i+1;
    end
    
    
        %Plot
    for r=1:3
        subplot(4,3,3*(vc-1)+ r);
        imagesc( flipud(dP(:,r)) );
        caxis([-1 1]*0.8);
        axis equal;
        xlim([0.5 2.5]);

        set(gca,'box','off','xtick','','xcolor','w','ytick','','ycolor','w');
        
        if r == 1
            set(gca,'ycolor','k'); ylabel(visLabels{vc});
        end
        
        if vc == size(visContrast,1)
            set(gca,'xcolor','k'); xlabel(respLabels{r});
        end
    end
    
end
colormap(cmap);



%% Load datasets
o=omnibusLaserGLM('sparse_unilateral_2D',{'Spemann','Whipple','Morgan','Murphy','Chomsky'});
d1 = o.data{end};
d1.areaIdx(d1.laserIdx==0) = 0;
d1 = getrow(d1,d1.areaIdx<=6);
d1.bilateral = ones(size(d1.response))*0;

o2=omnibusLaserGLM('sparse_bilateral_2D',{'Spemann','Whipple','Morgan','Murphy'});
d2 = o2.data{end};
d2.areaIdx(d2.laserIdx==0) = 0;
d2 = getrow(d2,d2.areaIdx<=6);
d2.bilateral = ones(size(d2.response));

d3 = struct;
d3.stimulus = [d1.stimulus; d2.stimulus];
d3.response = [d1.response; d2.response];
d3.areaIdx = [d1.areaIdx; d2.areaIdx];
d3.bilateral = [d1.bilateral; d2.bilateral];

d4=getrow(d3,d3.areaIdx==0);

%% For a given dataset, go through each trial and get prediction of that trial

objective = @(vec) -simulateOptimise('loglik',d1,vec) + 0*sum(vec.^2);
LB = [-inf(1,13) 0];
UB = [inf(1,13) 0];

options = optiset('display','final','solver','NOMAD');
Opt = opti('fun',objective,'x0',zeros(1,14),'bounds',LB,UB,'options',options);
[vec,~,exitflag] = Opt.solve;

p = simulateOptimise('pvec2struct',vec);

if exitflag<1
    warning('Did not converge');
end