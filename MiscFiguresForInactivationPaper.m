%% Plot contingency table for 2D AUC task

 choice_cols = [0 0.4470 0.7410;
                0.8500 0.3250 0.0980;
                0.4940    0.1840    0.5560;
                0.4660    0.6740    0.1880];
     
cols = zeros(4,4,3);
for cl = 1:4
    for cr = 1:4
        
        if cl==1 && cr==1
            cols(cl,cr,:) = choice_cols(3,:);
        end
        
        if cl==cr && cl>1
            cols(cl,cr,:) = choice_cols(4,:);
        end
        
        if cl>cr
            cols(cl,cr,:) = choice_cols(1,:);
        end
        
        if cl<cr
            cols(cl,cr,:) = choice_cols(2,:);
        end
    end
end

figure;
image(cols);
set(gca,'ydir','normal','box','off');
set(gca,'XTickLabel',{'0','0.1','0.24','0.54'},'xtick',1:4,'yTickLabel',{'0','0.1','0.24','0.54'},'ytick',1:4)
xlabel('CR'); ylabel('CL'); axis square;
set(gcf,'color','w'); set(gca,'fontsize',20);

%%  black+white grid of rewarded choices
figure('color','w');
subplot(1,3,1);
Lcorr = (0.5*(triu(ones(4)) + triu(ones(4),1)))'; Lcorr(1)=0;
imagesc(Lcorr); caxis([0 1]);
set(gca,'XTickLabel',{'0','0.1','0.24','0.54'},'xtick',1:4,'yTickLabel',{'0','0.1','0.24','0.54'},'ytick',1:4)
set(gca,'ydir','normal','box','off');
xlabel('Contrast Right'); ylabel('Contrast Left'); axis square; title('Left choice')

subplot(1,3,2);
Rcorr = (0.5*(tril(ones(4)) + tril(ones(4),-1)))'; Rcorr(1)=0;
imagesc(Rcorr); caxis([0 1]);
set(gca,'xtick','','ytick','');
set(gca,'ydir','normal','box','off'); axis square; title('Right choice')

subplot(1,3,3);
imagesc([1 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0]); caxis([0 1]);
set(gca,'xtick','','ytick',''); 
set(gca,'ydir','normal','box','off'); axis square; title('NoGo')
colormap(flipud(gray));
set(get(gcf,'children'),'fontsize',20);

%% simple RT histogram
figure('color','w');
hist(o.data{end}.RT(o.data{end}.laserIdx==0 & o.data{end}.response<3),100);
set(gca,'box','off','ytick','','ycolor','w','fontsize',20); xlabel('Reaction Time (sec)');
xlim([0 1.5]);
%% Plot grid of inactivation points

idx = 1:size(o.inactivationCoords,1);
coordNum = 42;
idx(coordNum) = [];

kimg=imread('D:\kirkcaldie_brain_BW.png');
figure('color','w');
subplot(1,2,1); %unilateral

idx = 1:size(o.inactivationCoords,1);
coordNum = 42;
imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal'); hold on;
set(imX,'AlphaData',0.5);
s1=scatter(o.inactivationCoords(idx,2),o.inactivationCoords(idx,1),200);
s2=scatter(o.inactivationCoords(coordNum,2),o.inactivationCoords(coordNum,1),200,'filled');
set([s1, s2],'CData',[0 0 1]);
s1.LineWidth=1; s1.CData=[0.4 0.4 1];
title([num2str(length(o.names)-1) ' mice, ' num2str(length(o.data{end}.response)) ' trials']);
xlabel('Medial-Lateral axis'); ylabel('Anterior-Posterior axis');

coordNum = [42 45];
subplot(1,2,2);
imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal'); hold on;
set(imX,'AlphaData',0.5);
s1=scatter(o.inactivationCoords(idx,2),o.inactivationCoords(idx,1),200);
s2=scatter(o.inactivationCoords(coordNum,2),o.inactivationCoords(coordNum,1),200,'filled');
set([s1, s2],'CData',[0 0 1]);
s1.LineWidth=1; s1.CData=[0.4 0.4 1];
title([num2str(length(o2.names)-1) ' mice, ' num2str(length(o2.data{end}.response)) ' trials']);
xlabel('Medial-Lateral axis'); ylabel('Anterior-Posterior axis');

%% Plot dots onto kirkcaldie brain to determine the right colours
img = imread('D:\kirkcaldie_brain.png');
figure; imagesc(linspace(-5,5,100),linspace(4,-6,100),img);
set(gca,'ydir','normal'); hold on;

O = o2;
numSites = size(O.inactivationCoords,1);
s=scatter(O.inactivationCoords(:,2),O.inactivationCoords(:,1),400,'filled');
text(O.inactivationCoords(:,2),O.inactivationCoords(:,1),cellfun(@(c)num2str(c),num2cell(1:numSites),'uni',0));


cols = [69, 198, 234;
        65, 140, 202;
        234, 155, 196;
        176, 115, 175;
        182, 216, 150;
        252, 200, 102;
        245, 242, 159;
        243, 237, 73;
        0.5 0.5 0.5]/255;
    
    %unilateral
areaID = [1 1 2 3 3 2 1 1,...
          1 1 2 3 3 2 1 1,...
          5 4 4 3 3 4 4 5,...
          5 5 5 3 3 5 5 5,...
          5 5 6 7 7 6 5 5,...
          6 6 7 7 6 6,...
          7 7 7 7,...
          8 8,...
          9];

%     %bilateral
areaID = [3 2 1 1,...
          3 2 1 1,...
          3 4 4 5,...
          3 5 5 5,...
          7 6 5 5,...
          7 6 6,...
          7 7,...
          8];
s=scatter(O.inactivationCoords(:,2),O.inactivationCoords(:,1),400,cols(areaID,:),'filled');
s.MarkerEdgeColor=[0 0 0];
text(O.inactivationCoords(:,2),O.inactivationCoords(:,1),cellfun(@(c)num2str(c),num2cell(1:numSites),'uni',0));

%% See if bilateral effects can be explained as a sum of unilateral effects

unilat_p=[];
for site = 1:size(o2.inactivationCoords,1)
    coord = o2.inactivationCoords(site,:);
    
    idx = (o.inactivationCoords(:,1) == coord(1) & o.inactivationCoords(:,2) == coord(2)) ...
        | (o.inactivationCoords(:,1) == coord(1) & o.inactivationCoords(:,2) == -coord(2));
    
%     bilat_p = o2.fitData.params{end}(site,:);
    unilat_p(site,:) = sum(o.fitData.params{end}(idx,:),1);
end

figure('color','w');
labels = {'bL','bR','sL','sR'};
for p=1:4
subplot(2,2,p);
plot(unilat_p(:,p),o2.fitData.params{end}(:,p),'o'); axis equal;
xl = xlim; yl = ylim;

hold on; ezplot('y=x');
xlim(xl); ylim(yl);
xlabel('Unilat L + Unilat R parameter');
ylabel('Bilateral parameter');
set(gca,'box','off');
title(labels{p});
end

%try comparing on the level of predicted probabilities

p_nL = mean(o2.fitData.nonLaserParams{end});
[ZL,ZR] = o2.getModel('bias+sub_sens','biasAsContrast',5);

c = [0.5 0.5];
figure('color','w');
a=subplot(1,2,1); hold(a,'on'); ezplot('y=x'); xlim([0 1]); ylim([0 1]); axis square; title('L'); xlabel('Unilat L+R'); ylabel('Bilat');
b=subplot(1,2,2); hold(b,'on'); ezplot('y=x'); xlim([0 1]); ylim([0 1]); axis square; title('R'); xlabel('Unilat L+R'); ylabel('Bilat');
for site = 1:size(o2.inactivationCoords,1)
    zl = ZL(p_nL,o2.fitData.params{end}(site,:),c);
    zr = ZR(p_nL,o2.fitData.params{end}(site,:),c);
    pLb = exp(zl)/(1+exp(zl)+exp(zr));
    pRb = exp(zr)/(1+exp(zl)+exp(zr));
    
    zl = ZL(p_nL,unilat_p(site,:),c);
    zr = ZR(p_nL,unilat_p(site,:),c);
    pLu = exp(zl)/(1+exp(zl)+exp(zr));
    pRu = exp(zr)/(1+exp(zl)+exp(zr));
    
    plot(a,pLu,pLb,'o');
    plot(b,pRu,pRb,'o');
end