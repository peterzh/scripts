
%% %%%%%% CREATE MAP OF LASER GLM MODEL LOG LIKELIHOOD DIFFERENCES %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subject = 'Eijkman';
expRefs = dat.listExps(subject);
saveDir = '\\basket.cortexlab.net\homes\peterzh\NeuroGLM-experimental\ModelFiles';
models = {'C^N-subset-nolaser','C^N-subset-laser','C^N-subset-laser-offset','C^N-subset-laser-scale'};

pSlopes = [];
pOffsets = [];

badFiles = {};
for f = 1:length(expRefs)
    try
        g = laserGLM(expRefs{f});
        
        for m = 1:length(models)
            g = g.setModel(models{m});

            for site = 1:size(g.inactivationSite,1)
                g1 = g.fitCV(site);
                save(fullfile(saveDir,[g1.expRef '_laser-' num2str(site) '_model-'  g1.modelString '.mat']),'g1');
            end
        end

    catch err
        %         keyboard;
        disp([expRefs{f} ': ' err.message])
        
%         if strcmp(err.message(1:20),'Reference to non-exi')
%             badFiles = [badFiles; expRefs{f}];
%         end
    end
end

disp('Done!');

%% Get and plot values
subject = 'Eijkman';
expRefs = dat.listExps(subject);
models = {'C^N-subset-nolaser','C^N-subset-laser','C^N-subset-laser-offset','C^N-subset-laser-scale'};
modelLabels = {'NoLaser','Full','Offset','Scale'};

laserCoord = [];
measure = [];
for f = 1:length(expRefs)
    for site = 1:6

        try
            NLL = [];
            for m = 1:length(models)
                load(fullfile(saveDir,[expRefs{f} '_laser-' num2str(site) '_model-'  models{m} '.mat']));
                NLL(1,m) = -sum(log2(g1.p_hat))/length(g1.data.response);
            end
            
             if g1.inactivationSite(1) <= -2 && g1.inactivationSite(2) >= +1.5
%                  keyboard;
                 disp(g1.expRef);
                 disp(site);
             end
% 
%             if strcmp(g1.expRef,'2015-10-07_1_Hopkins')
%                 keyboard;
%             end
            OUTPUT = NLL;
            
            if ~isinf(OUTPUT) 
                measure = [measure; OUTPUT];
                laserCoord = [laserCoord; g1.inactivationSite];
            end
            
        catch
%             warning('site does not exist');
        end

    end
end

% Plot values

% Check for duplicate sites 
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

allPairs = nchoosek(1:length(models),2);

figure;
for i = 1:length(allPairs)
    a = allPairs(i,1);
    b = allPairs(i,2);
    diff = measure(:,b) - measure(:,a);
    subplot(2,3,i);
    scatter(newlaserCoord(:,2),newlaserCoord(:,1),dotsize,diff,'filled');
    colormap(cmap);
    hold on; scatter(newlaserCoord(:,2),newlaserCoord(:,1),dotsize,'MarkerEdgeColor',[0 0 0]); hold off
    caxis([-max(abs(diff)) max(abs(diff))])
    axis equal
    title(['{\color{red}' modelLabels{b} '}' ' vs ' '{\color{blue}' modelLabels{a} '}'], 'FontSize', 20);
    colorbar;
end

%% %%%%%% CREATE MAP OF MISCELLANEOUS SESSION INFORMATION %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subject = 'Eijkman';
expRefs = dat.listExps(subject);

laserCoord = [];
measure = [];
for f = 1:length(expRefs)
    for site = 1:6

        try
            load(fullfile(saveDir,[expRefs{f} '_laser-' num2str(site) '_model-'  'C^N-subset-nolaser' '.mat']));
            D = getrow(g1.data,g1.data.laserIdx==0);
            nogo = D.contrast_cond(:,1) == 0 & D.contrast_cond(:,2) == 0;
            tab = tabulate(D.response(nogo==1));
            
            OUTPUT = tab(1:2,3)';
            
            if ~isinf(OUTPUT)
                measure = [measure; OUTPUT];
                laserCoord = [laserCoord; g1.inactivationSite];
            end
            
        catch
        end

    end
end

% Check for duplicate sites and average the measure at those sites
%Check for repeat sites and take average if so
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

labels = {'Left FA [%]','Right FA [%]'};
dotsize = 200;
figure;
for i = [1,2]
    subplot(1,2,i);
    scatter(newlaserCoord(:,2),newlaserCoord(:,1),dotsize,measure(:,i),'filled');
    hold on; scatter(newlaserCoord(:,2),newlaserCoord(:,1),dotsize,'MarkerEdgeColor',[0 0 0]); hold off
    caxis([0 max(abs(measure(:,i)))])
    axis equal
    colorbar;
    title(labels{i});
end
