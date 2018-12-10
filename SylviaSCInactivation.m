clear all; close all;

load('D:\Downloads\SC_inactivation_data.mat');
outcome = rmfield(outcome, 'date');
outcome.stimulus = outcome.contrast;
outcome = structfun(@(f) f(1:length(outcome.response),:), outcome, 'uni', 0);

%% Pooled session plots

laserTypes = [0 1; 1 0; 1 1];
laserTypeLabels = {'Right','Left','Bilat'};

for i = 1:3
    nonLaser = getrow(outcome, outcome.laserType(:,1)==0 & outcome.laserType(:,2)==0);
    fig1 = GLM(nonLaser).setModel('C50-subset').fit.plotFit;

    Laser = getrow(outcome, outcome.laserType(:,1)==laserTypes(i,1) & outcome.laserType(:,2)==laserTypes(i,2));
    fig2 = GLM(Laser).setModel('C50-subset').fit.plotFit;
    
    axs1 = get(fig1,'children');
    axs2 = get(fig2,'children');
    
    for ax = 1:length(axs1)
        nLstuff = get(axs1(ax),'children');
        set(nLstuff,'color',[0 0 0]);
        
        
        Lstuff = get(axs2(ax),'children');
        set(Lstuff,'color',[1 0 0]);
        copyobj(Lstuff, axs1(ax))
    end
    
    close(fig2)
    set(fig1,'name',laserTypeLabels{i});
end

%% Per session

laserTypes = [0 1; 1 0; 1 1];
laserTypeLabels = {'Right','Left','Bilat'};

for sess = 1:max(outcome.session)
    outcome_sess = getrow(outcome, outcome.session==sess);
    
    for i = 1:3
        nonLaser = getrow(outcome_sess, outcome_sess.laserType(:,1)==0 & outcome_sess.laserType(:,2)==0);
        fig1 = GLM(nonLaser).setModel('C50-subset').fit.plotFit;
        
        Laser = getrow(outcome_sess, outcome_sess.laserType(:,1)==laserTypes(i,1) & outcome_sess.laserType(:,2)==laserTypes(i,2));
        
        if ~isempty(Laser.response)
            fig2 = GLM(Laser).setModel('C50-subset').fit.plotFit;
            
            axs1 = get(fig1,'children');
            axs2 = get(fig2,'children');
            
            for ax = 1:length(axs1)
                nLstuff = get(axs1(ax),'children');
                set(nLstuff,'color',[0 0 0]);
                
                
                Lstuff = get(axs2(ax),'children');
                set(Lstuff,'color',[1 0 0]);
                copyobj(Lstuff, axs1(ax))
            end
            
            close(fig2)
            set(fig1,'name',['Session ' num2str(sess) '  Inactivation: ' laserTypeLabels{i}]);

        else
            close(fig1);
        end
    end
end