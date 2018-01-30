%Script for checking the effects of a light baffle

%1.5mW 1.5seconds  DC laser sessions
sessions = {'2017-10-26_2_Nyx'
    '2017-10-26_1_Vin';
    '2017-10-26_1_Keynes';
'2017-10-27_1_Keynes';
'2017-10-27_4_Nyx';
'2017-10-27_1_Vin';
'2017-10-30_2_Vin';
'2017-10-30_1_Keynes';
'2017-11-01_3_Nyx';
'2017-11-01_1_Vin';
'2017-11-01_1_Keynes';
'2017-11-03_1_Nyx';
'2017-11-03_1_Vin';
'2017-11-03_1_Keynes';
};

%1.5mW 40Hz 1.5seconds sine wave sessions
laserRegion = {'LeftVIS','RightVIS','LeftM2','RightM2','LeftFrontOutside','RightFrontOutside','LeftBackOutside','RightBackOutside'};
figure('color','w');

Dall = struct;
for sess = 1:length(sessions)
    D = loadData(sessions{sess});
    
    %Only look at zero contrast conditions
    D = getrow(D, D.stimulus(:,1)==0 & D.stimulus(:,2)==0);
%     D = getrow(D, D.stimulus(:,1)==D.stimulus(:,2) & D.stimulus(:,1)>0);
%     
    nL = getrow(D,D.laserType==0);
    [pNL, pci] = binofit(sum(nL.response==[1 2 3],1), length(nL.response));
%     
%     subplot(length(sessions),3,3*sess - 2); hold on;
%     lx = line([0 length(laserRegion)+1],[1 1]*p(1));
%     lx.LineStyle = '--';
% %     hx = errorbar(0,p(1),p(1)-pci(1,1),pci(1,2)-p(1)); hx.Marker='.'; hx.MarkerSize=20;
%     ylabel(sessions{sess}); 
    
%     subplot(length(sessions),3,3*sess - 1); hold on;
%     lx = line([0 length(laserRegion)+1],[1 1]*p(2));
%     lx.LineStyle = '--';
%     
%     subplot(length(sessions),3,3*sess - 0); hold on;
%     lx = line([0 length(laserRegion)+1],[1 1]*p(3));
%     lx.LineStyle = '--';
    
    for region = 1:length(laserRegion)
        
        subplot(length(laserRegion),1,region); hold on;
        
        L = getrow(D,D.laserRegion==laserRegion{region});
        p = binofit(sum(L.response==[1 2 3],1), length(L.response));
        
        plot(sess, p(2)-pNL(2),'ro');
        plot(sess, p(1)-pNL(1),'bo');
        
%         subplot(length(sessions),3,3*sess - 2);
%         hx = errorbar(sess,p(1),p(1)-pci(1,1),pci(1,2)-p(1)); hx.Marker='.'; hx.MarkerSize=20;
%         hx.Color=[0 0 0];
%         subplot(length(sessions),3,3*sess - 1);
%         hx = errorbar(region,p(2),p(2)-pci(2,1),pci(2,2)-p(2)); hx.Marker='.'; hx.MarkerSize=20;
%         
%         subplot(length(sessions),3,3*sess - 0);
%         hx = errorbar(region,p(3),p(3)-pci(3,1),pci(3,2)-p(3)); hx.Marker='.'; hx.MarkerSize=20;
        ylabel(laserRegion{region});
    end
    
    drawnow;
    
    
    Dall = addstruct(Dall,D);
end
axs = get(gcf,'children');

for ax = 1:length(axs)
    line(axs(ax),[axs(ax).XLim],[0 0]);
end
set(axs,'ylim',[-1 1]*0.4);
% set(axs(1:3),'xtick',[1:length(laserRegion)],'xticklabel',laserRegion,'xticklabelrotation',45);
% 
% 
% 
% figure('color','w');
% nL = getrow(Dall,Dall.laserType==0);
% [p, pci] = binofit(sum(nL.response==[1 2 3],1), length(nL.response));
% 
% subplot(1,3,1); hold on;
% lx = line([0 length(laserRegion)+1],[1 1]*p(1));
% lx.LineStyle = '--';
% 
% subplot(1,3,2); hold on;
% lx = line([0 length(laserRegion)+1],[1 1]*p(2));
% lx.LineStyle = '--';
% 
% subplot(1,3,3); hold on;
% lx = line([0 length(laserRegion)+1],[1 1]*p(3));
% lx.LineStyle = '--';
% 
% for region = 1:length(laserRegion)
%     L = getrow(Dall,Dall.laserRegion==laserRegion{region});
%     [p, pci] = binofit(sum(L.response==[1 2 3],1), length(L.response));
%     
%     subplot(1,3,1);
%     hx = errorbar(region,p(1),p(1)-pci(1,1),pci(1,2)-p(1)); hx.Marker='.'; hx.MarkerSize=20;
%     
%     subplot(1,3,2);
%     hx = errorbar(region,p(2),p(2)-pci(2,1),pci(2,2)-p(2)); hx.Marker='.'; hx.MarkerSize=20;
%     
%     subplot(1,3,3);
%     hx = errorbar(region,p(3),p(3)-pci(3,1),pci(3,2)-p(3)); hx.Marker='.'; hx.MarkerSize=20;
%     
% end
% axs = get(gcf,'children');
% set(axs,'ylim',[0 1]);
% set(axs(1:3),'xtick',[1:length(laserRegion)],'xticklabel',laserRegion,'xticklabelrotation',45);

