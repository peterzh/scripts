% Prerequisites: expRefs and laserGLM objects from the fitting script
%% Averaging over individual samples
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
        plot(c0_pred); if session==1; title('C=0'); end;
        hold on; errorbar(1:3,c0_act,c0_act-c0_act_ci(:,1),c0_act_ci(:,2)-c0_act,'rs'); hold off;
        subplot(numSessions,4,4*session - 2);
        plot(cE(:,1)); if session==1; title('CL=CR');end;
        hold on; plot(cE(:,2),'s');hold off;
        subplot(numSessions,4,4*session - 1);
        plot(cL(:,1)); if session==1; title('C=L');end;
        hold on; plot(cL(:,2),'s');hold off;
        subplot(numSessions,4,4*session);
        plot(cR(:,1)); if session==1; title('C=R');end;
        hold on; plot(cR(:,2),'s');hold off;
        
        %             errorbar(1:3,p,p-pci(:,1),pci(:,2)-p,'rs');
        ylim([0 1]);
        ax = get(gcf,'children');
        set(ax,'Xticklabel','');
        set(ax(1),'Xticklabel',{'L','R','NG'},'Xtick',1:3);
        %             avStim = mean(l(s).data.contrast_cond(cF==stim,:));
        %             avStim = g.calculatePhat(params,l(s).data.contrast_cond(cF==stim,:));
        %             avP = mean(g.calculatePhat(params,l(s).data.contrast_cond(cF==stim,:)));
        %             hold on;
        %             plot(avP);
        %             hold off;
        
    end
end

%% Combining over stimulus quadrants
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
        
        D.p_hats = g.calculatePhat(params,D.contrast_cond);
        
        stim_labels = {'C=0','CL=CR','C=L','C=R'};
        stim = cell(1,4);
        stim{1} = sum(D.contrast_cond,2)==0;
        stim{2} = ~stim{1} & (D.contrast_cond(:,1)==D.contrast_cond(:,2));
        stim{3} = diff(D.contrast_cond,[],2)<0;
        stim{4} = diff(D.contrast_cond,[],2)>0;
        
        for stim_type = 1:4
            r = D.response(stim{stim_type});
            tab = [sum(r==1) sum(r==2) sum(r==3)]';
            [p,pci] = binofit(tab,sum(tab));
            
            subplot(numSessions,4,4*session - 4 + stim_type);
            errorbar(1:3,p,p-pci(:,1),pci(:,2)-p,'rs'); %plot actual
            hold on;
            plot(mean(D.p_hats(stim{stim_type},:)));
            
            hold off;
            
            if session==1
                title(stim_labels{stim_type});
            end
        end
    end
    
    ax=get(gcf,'children');
    set(ax,'Xticklabel','');
    set(ax(1),'Xticklabel',{'L','R','NG'},'xtick',1:3);
end