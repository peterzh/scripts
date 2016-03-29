%% Grab data
subject = {'Eijkman'};
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
        
        alldat = addstruct(alldat,g.data);
        disp(expRefs{f});
        
    catch
        %         warning('error');
    end
end
D = getrow(alldat,alldat.laserIdx==0);
disp('Done!');

%% Grab subsets and plot
HISTORY = 1;

contrast_groups = [-0.6 -0.3 -0.001 0.001 0.3 0.6];
ic = discretize(diff(D.contrast_cond,[],2),contrast_groups);

hist = struct;
hist.contrast_group = circshift(ic,HISTORY);
hist.response = circshift(D.response,HISTORY);
hist.feedbackType = circshift(D.feedbackType,HISTORY);

for feedback = [-1 1]
    figure;
    a=1;
    for prevR = 1:3
        for prevC = 1:(length(contrast_groups)-1)
            try
                E = getrow(D,hist.feedbackType==feedback & hist.contrast_group==prevC & hist.response==prevR);
                
                subplot(3,(length(contrast_groups)-1),a);
                GLM(E).setModel('C^N-subset').fit.plotFit;
            catch
                warning('no trials of this type');
            end
            a=a+1;
        end
    end
end
