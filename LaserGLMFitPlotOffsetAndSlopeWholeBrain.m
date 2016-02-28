
%% %%%%%%%%%%%%% CREATE MAP OF GLM FIT SLOPE&OFFSET FROM GLM FITS  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% go through many sessions of a mouse, fit GLM with and without laser and plot offset and slope at c=0
subject = 'Eijkman';
expRefs = dat.listExps(subject);
model = 'C^N-subset-laser';

laserCoord = [];
pSlopes = [];
pOffsets = [];

badFiles = {};
for f = 1:length(expRefs)
    try
        g = laserGLM(expRefs{f}).setModel(model);
        
        for site = 1:size(g.inactivationSite,1)
            g1 = g.fit(site);
        
            %get non-laser offset
            offset_nolaser = g1.calculatePhat(g1.parameterFits,[0 0 0]);
            offset_laser = g1.calculatePhat(g1.parameterFits,[0 0 1]);
            
            %get non-laser slope
            slopes_nolaser = calculateSlopeC0(g1,0);
            slopes_laser = calculateSlopeC0(g1,1);
            
            laserCoord = [laserCoord; g1.inactivationSite];

            pOffsets = [pOffsets; 100*offset_laser(1:2)./offset_nolaser(1:2)];

            pSlopes = [pSlopes; 100*slopes_laser./slopes_nolaser];
        end
        
    catch err
        %         keyboard;
        disp([expRefs{f} ': ' err.message])
        
        if strcmp(err.message(1:20),'Reference to non-exi')
            badFiles = [badFiles; expRefs{f}];
        end
    end
    
end

disp('Done!');

%save to file
save(fullfile(saveDir,[subject '_Slope&OffsetMap.mat']),'laserCoord','pOffsets','pSlopes');

%% Plot 
subject = 'Eijkman';
load(fullfile(saveDir,[subject '_Slope&OffsetMap.mat']));

newlaserCoord = laserCoord;
if length(unique(laserCoord,'rows')) < length(laserCoord)
    [uniqueSites,~,siteID]=unique(laserCoord,'rows');
     
    for i = 1:length(uniqueSites)
        idx=find(laserCoord(:,1) == uniqueSites(i,1) & laserCoord(:,2) == uniqueSites(i,2));
        if length(idx) > 1
            for c = 1:length(idx)
                [x,y]=pol2cart(c*2*pi./length(idx),0.1);
                newlaserCoord(idx(c),1) = laserCoord(idx(c),1) + y;
                newlaserCoord(idx(c),2) = laserCoord(idx(c),2) + x;
            end
        end
    end
end

cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
        linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];

dotsize = 200;

figure; 

subplot(2,2,1);
scatter(newlaserCoord(:,2),newlaserCoord(:,1),dotsize,pOffsets(:,1),'filled')
title('Left');
ylabel('% of non-laser Offset');
colormap(cmap);
hold on; scatter(newlaserCoord(:,2),newlaserCoord(:,1),dotsize,'MarkerEdgeColor',[0 0 0]); hold off
caxis([0 200])
axis equal

subplot(2,2,2);
scatter(newlaserCoord(:,2),newlaserCoord(:,1),dotsize,pOffsets(:,2),'filled')
colormap(cmap);
hold on; scatter(newlaserCoord(:,2),newlaserCoord(:,1),dotsize,'MarkerEdgeColor',[0 0 0]); hold off
title('Right');
caxis([0 200])
axis equal


subplot(2,2,3);
scatter(newlaserCoord(:,2),newlaserCoord(:,1),dotsize,pSlopes(:,1),'filled')
colormap(cmap);
hold on; scatter(newlaserCoord(:,2),newlaserCoord(:,1),dotsize,'MarkerEdgeColor',[0 0 0]); hold off
ylabel('% of non-laser Slope');
caxis([0 200])
axis equal

subplot(2,2,4);
scatter(newlaserCoord(:,2),newlaserCoord(:,1),dotsize,pSlopes(:,2),'filled')
colormap(cmap);
hold on; scatter(newlaserCoord(:,2),newlaserCoord(:,1),dotsize,'MarkerEdgeColor',[0 0 0]); hold off
caxis([0 200]);
axis equal