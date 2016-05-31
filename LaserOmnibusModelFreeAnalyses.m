% Prerequisites: expRefs and laserGLM objects from the fitting script
g=g.setModel('C50');

for s = 1:numSubjects
    c0=sum(l(s).data.contrast_cond,2)==0;
    cL=diff(l(s).data.contrast_cond,[],2)<0;
    cR=diff(l(s).data.contrast_cond,[],2)>0;
        
    figure('name',expRefs{s,1});
    numSessions = max(l(s).data.sessionID);
    for session = 1:numSessions
        cF = (1*cL + 2*cR + 3*c0).*(l(s).data.sessionID==session).*(l(s).data.laserIdx==0);
        
        params = [sess(s).Biases(session,1),...
                  sess(s).CLeft(session,1),...
                  sess(s).CRight(session,1),...
                  sess(s).Biases(session,2),...
                  sess(s).CLeft(session,2),...
                  sess(s).CRight(session,2),...
                  c50_fits{s}(session,:)];
        %Actual choices        
        for stim = 1:3
            choices = l(s).data.response(cF==stim);
            tab=[sum(choices==1) sum(choices==2) sum(choices==3)]';
            [p,pci]=binofit(tab,length(choices));
            
           
            subplot(numSessions,3,3*session - 3 + stim);
            errorbar(1:3,p,p-pci(:,1),pci(:,2)-p,'rs');
%             plot(p,'s'); 
            ylim([0 1]);
            
            avStim = mean(l(s).data.contrast_cond(cF==stim,:));
            hold on;
            plot(g.calculatePhat(params,avStim));
            hold off;
        end
    end
end