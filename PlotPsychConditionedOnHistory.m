%% Grab data
subject = {'Hopkins','Eijkman','Murphy','Spemann','Morgan','Whipple'};
expRefs = dat.listExps(subject)';
expRefs = vertcat(expRefs{:});

alldat = struct;
for f = 1:length(expRefs)
    try
        try
            g = laserGLM(expRefs{f});
        catch
            g = GLM(expRefs{f});
            g.data.laserIdx = zeros(length(g.data.response),1);
            g.data.laser = nan(length(g.data.response),2);
        end
        
        if length(g.data.laserIdx) ~= length(g.data.response)
            error('Mismatch laser + behav data length');
        end
        
        g.data = structfun(@(f)f(1:end-20,:),g.data,'uni',0); %remove last 20 trials
        
        alldat = addstruct(alldat,g.data);
        disp(expRefs{f});
        
    catch
        %         warning('error');
    end
end
D = getrow(alldat,alldat.laserIdx==0);
disp('Done!');

remove = D.contrast_cond(:,1)>0 & D.contrast_cond(:,2)>0;
D=structfun(@(f)f(~remove,:),D,'uni',0);


%% Grab subsets and plot
HISTORY = 3;

contrast_groups = [-1.1 -0.5 -0.001 0.001 0.5 1.1];
contrast_groups_labs = {'High L','Low L','Zero','Low R','High R'};
response_labs = {'Chose L','Chose R','NG'};
ic = discretize(diff(D.contrast_cond,[],2),contrast_groups);

hist = struct;
hist.contrast_group = circshift(ic,HISTORY);
hist.response = circshift(D.response,HISTORY);
hist.feedbackType = circshift(D.feedbackType,HISTORY);

figure;
a=1;
for prevR = [2 3 1]
    for prevC = 1:(length(contrast_groups)-1)
        try
            E = getrow(D,hist.contrast_group==prevC & hist.response==prevR);
            
            subplot(3,(length(contrast_groups)-1),a);
            GLM(E).setModel('C^N-subset').fit.plotFit;
            ax=get(gca,'children');
            set(ax(1:3),'linewidth',3);
%             ax(1).delete;
            ax(4).delete;
            ax(5).LData=zeros(1,length(ax(5).XData));
            ax(5).UData=zeros(1,length(ax(5).XData));

            ax(6).LData=zeros(1,length(ax(6).XData));
            ax(6).UData=zeros(1,length(ax(6).XData));

%             ax(5).delete;
%             ax(6).delete;
            xlim([-1 1]);
            
            h=getrow(hist,hist.contrast_group==prevC & hist.response==prevR);
            if h.feedbackType(1)==1
                set(gca,'Color',[0 1 0 0.05]);
            else
                set(gca,'Color',[1 0 0 0.05]);
            end
            
            set(gca,'box','off');
            xlabel(''); ylabel('');
            title(['n=' num2str(length(E.response))]);
            %                 keyboard;
            
            if prevC == 1
                ylabel(response_labs{prevR})
            end
            
            if prevR == 1
                xlabel(contrast_groups_labs{prevC});
            end
        catch
        end
        a=a+1;
    end
end

set(gcf,'color','white');

