%Reassign choices based on whether there was an acceleration component
%after stimulus onset

%variable=dat
dat.new_response = dat.response;
for i = 1:length(dat.response)
    t = dat.wheel_stimulusOn_timesteps(i,:);
    p = dat.wheel_stimulusOn(i,:);

    p1 = smooth(p,50,'lowess');
    accel = gradient(gradient(p1));
    
    
    if dat.response(i)<3
%         plotyy(t,p1,t,accel);
        
        if dat.response(i)==1
            tMax = t( find(accel==max(accel)) );
        elseif dat.response(i)==2
            tMax = t( find(accel==min(accel)) );
        end
        
%         hold on;
%         lx=line([tMax tMax], [-1 1]*1000); lx.Color = [0 0 0];
%         hold off;
        
        if tMax<0.05
            dat.new_response(i) = 3;
%              title([num2str(dat.response(i)) ' -> 3'])
        else
%              title(dat.response(i))
        end
        
%         keyboard;
    end
end