function decoding(eRef)
    %Performs choice decoding, and saves decoding data for plotting later
    numT = 100;
    epoch_dt = linspace(-1,+1,numT);
    
    load( fullfile('..','preproc',[eRef '.mat']) );
    
    %Crossvalidate behavioural model
    D = struct('stimulus',behav.stimulus,'response',double(behav.response));
    g = GLM(D).setModel('C50-subset').fitCV(10);
    pL = g.p_hat(:,1);    pR = g.p_hat(:,2);    pNG = g.p_hat(:,3);
    
    bpt_baseline = struct;
    likelihood_Full = pL.*(g.data.response==1) + pR.*(g.data.response==2) + pNG.*(g.data.response==3);
    bpt_baseline.full = mean( log2( likelihood_Full ) );
    likelihood_GOvNOGO = (1 - pNG).*(g.data.response<3) + pNG.*(g.data.response==3);
    bpt_baseline.GOvNOGO = mean( log2( likelihood_GOvNOGO ) );
    
    likelihood_LvR = (pL./(1-pNG)).*(g.data.response==1) + (pR./(1-pNG)).*(g.data.response==2) + 10*(g.data.response==3);
    likelihood_LvR(likelihood_LvR==10) = [];
    bpt_baseline.LvR = mean( log2( likelihood_LvR ) );
    
    %Centre and whiten each neuron's firing rates
%     spikes.firingRatesNorm = zscore(spikes.firingRates')';
    
    brainRegions = unique(population.region);
    
    tab = table;
    tab.Y = D.response;
    tab.offset = behav.GLM_offset;
    
    
    bpt = struct;
    bpt.full = nan(numT, length(epoches), length(brainRegions));
    bpt.GOvNOGO = nan(numT, length(epoches), length(brainRegions));
    bpt.LvR = nan(numT, length(epoches), length(brainRegions));
    
    numberOfRuns = length(brainRegions)*length(epoches)*numT;
    
    tic;
    i = 1;
    for region = 1:length(brainRegions)
        fprintf('Region: %d/%d\n',region,length(brainRegions))
        
        pop = population( strcmp(population.region,brainRegions{region}) ,:);
        %Compute epoched firing rate for each neuron
        for ep = 1:length(epoches)
            fprintf('\tEpoch: %d/%d\n',ep,length(epoches));
            t_zero = behav.(epoches{ep});
            for t = 1:numT %for each timestep in that epoch
                fprintf('\t\tTime: %d/%d',t,numT);
                queryT = t_zero + epoch_dt(t);
                
                %Interpolate to get the firingRate of all neurons at this time
                tab.X = interp1(pop.Properties.UserData.fr_times ,pop.firingRates',queryT);
                tab.X = [tab.X zeros(size(tab.X))];
                
                %Remove any entries with NaN
                tabCleaned = tab(~isnan(tab.X(:,1)),:);
                
%                 fit = glmnet(tab.X,tab.Y,'multinomial',glmnetSet(struct('offset',tab.offset)));
                
                %Using this to decode choice!
                writetable(tabCleaned,'C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\data.csv');
                exitflag = system('"C:\Program Files\R\R-3.4.0\bin\R" CMD BATCH C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\glmnet_crossval.R');
                assert(exitflag==0,'Error running GLMNET via R interface');
                
                phat = dlmread('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\phat.dat',',');
                lik = phat(:,1).*(tab.Y==1) + phat(:,2).*(tab.Y==2) + phat(:,3).*(tab.Y==3);
                bpt.full(t, ep, region) = mean(log2(lik));
                
                lik = (1 - phat(:,3)).*(tab.Y<3) + phat(:,3).*(tab.Y==3);
                bpt.GOvNOGO(t, ep, region) = mean(log2(lik));
                
                lik = (phat(:,1)./(1-phat(:,3))).*(tab.Y==1) + (phat(:,2)./(1-phat(:,3))).*(tab.Y==2) + 10*(tab.Y==3);
                lik(lik==10) = [];
                bpt.LvR(t, ep, region) = mean(log2(lik));
                
%                 bpt.full(t, ep, region) = dlmread('C:\Users\Peter\Documents\MATLAB\GLMNET_DATA\bpt.dat','\t');
                fprintf('\t%6.3f',bpt.full(t, ep, region)-bpt_baseline.full);
                
                %Print very approximate time left of computation
                timePerRun = toc/i;
                timeLeft = timePerRun*( numberOfRuns-i );
                timeLeft = datestr(timeLeft/(60*60*24),'HH:MM:SS');
                fprintf('\t Time Left: %s\n',timeLeft);
                i = i+1;
            end
            pause(0.5);
        end
    end
    save(['../decoding/' eRef '.mat'],'bpt_baseline','brainRegions','epoches','epoch_dt','bpt')
end