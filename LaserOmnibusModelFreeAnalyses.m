% Prerequisites: expRefs and laserGLM objects from the fitting script
g=g.setModel('C50');

for s = 1:numSubjects
    
        
    figure('name',expRefs{s,1});
    numSessions = max(l(s).data.sessionID);
    for session = 1:numSessions
        D = getrow(l(s).data,l(s).data.laserIdx==0 & l(s).data.sessionID==session);
        params = [sess(s).Biases(session,1),...
                  sess(s).CLeft(session,1),...
                  sess(s).CRight(session,1),...
                  sess(s).Biases(session,2),...
                  sess(s).CLeft(session,2),...
                  sess(s).CRight(session,2),...
                  c50_fits{s}(session,:)];
              
%         D.p_hats = g.calculatePhat(params,D.contrast_cond);
        C = unique(D.contrast_cond(:,1));
        
        actual = nan(length(C),length(C),3);
        actual_ci = cell(length(C),length(C));
        pred = nan(length(C),length(C),3);
        for cl=1:length(C)
            for cr=1:length(C)
                idx = D.contrast_cond(:,1)==C(cl) & D.contrast_cond(:,2)==C(cr);
                r = D.response(idx);
                tab=[sum(r==1) sum(r==2) sum(r==3)]';
                
                if ~isempty(r)
                    [p,pci] = binofit(tab,length(r));
                    actual(cl,cr,1:3) = p;
                    actual_ci{cl,cr} = pci;
                end
                pred(cl,cr,1:3) = g.calculatePhat(params,[C(cl) C(cr)]);%D.p_hats(find(idx,1,'first'),:);
            end
        end
        
        c0=[]; cE=[]; cL=[]; cR=[];
        c0_pred = squeeze(pred(1,1,:));
        c0_act =  squeeze(actual(1,1,:));
        c0_act_ci = actual_ci{1,1};
        for r = 1:3
            pre = pred(:,:,r);
            act = actual(:,:,r);
            cE(r,1) = mean(diag(pre(2:end,2:end)));
            cE(r,2) = nanmean(diag(act(2:end,2:end)));
            cL(r,1) = mean(pre([2,3,4,7,8,12]));
            cL(r,2) = nanmean(act([2,3,4,7,8,12]));
            cR(r,1) = mean(pre([5,9,10,13,14,15]));
            cR(r,2) = nanmean(act([5,9,10,13,14,15]));
        end

        %Actual choices        

           
            subplot(numSessions,4,4*session - 3);
            plot(c0_pred); title('C=0');
            hold on; errorbar(1:3,c0_act,c0_act-c0_act_ci(:,1),c0_act_ci(:,2)-c0_act,'rs'); hold off;
            subplot(numSessions,4,4*session - 2);
            plot(cE(:,1)); title('CL=CR');
            hold on; plot(cE(:,2),'s');hold off;
            subplot(numSessions,4,4*session - 1);
            plot(cL(:,1)); title('C=L');
            hold on; plot(cL(:,2),'s');hold off;
            subplot(numSessions,4,4*session);
            plot(cR(:,1)); title('C=R');
            hold on; plot(cR(:,2),'s');hold off;
            
%             errorbar(1:3,p,p-pci(:,1),pci(:,2)-p,'rs');
            ylim([0 1]);
            
%             avStim = mean(l(s).data.contrast_cond(cF==stim,:));
%             avStim = g.calculatePhat(params,l(s).data.contrast_cond(cF==stim,:));
%             avP = mean(g.calculatePhat(params,l(s).data.contrast_cond(cF==stim,:)));
%             hold on;
%             plot(avP);
%             hold off;
        
    end
end