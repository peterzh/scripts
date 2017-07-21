function [expRefs] = LaserIdentifySessions(subject,region)
% subject = {'Hopkins','Eijkman','Morgan','Whipple','Murphy','Spemann'};
expRefs = dat.listExps(subject)';
expRefs = vertcat(expRefs{:});

% region = 'RightVisual';

switch(region)
    case 'LeftVisual'
        AP = [-5.5 -2];
        ML = [-4 -1];
    case 'RightVisual'
        AP = [-5.5 -2];
        ML = [1 4];
    case 'LeftBarrel'
        AP = [-1.5 0.5];
        ML = [-4 -1.5];
    case 'RightBarrel'
        AP = [-1.5 0.5];
        ML = [1.5 4];
    otherwise
        error('Choose a region');
end

a=1;
nexpRefs={};

coords = [];

for b = 1:length(expRefs)
    try
        l=laserGLM(expRefs{b});
        
        for s=1:size(l.inactivationSite)
            site = l.inactivationSite(s,:);
            
            if AP(1) < site(2) && site(2) < AP(2) && ML(1) < site(1) && site(1) < ML(2)
                disp(expRefs{b});
                nexpRefs{a,1} = expRefs{b};
                nexpRefs{a,2} = s;
                
%                 coords = [coords; site];
                a=a+1;
            end
        end
    catch
    end
end

disp('Done!');
expRefs = nexpRefs;
figLabel = region;

% unique(coords(:,1:2),'rows')
end