function varargout = galvoPulseInactivation(what,varargin)

switch(what)
    case 'getData'
        otherSessionCriteria = varargin{1};
        S = d.Session & otherSessionCriteria & 'project_id="galvoPulse_unilateral"' & 'num_trials>150' & 'session_date>"2018-03-18"';
        expRefs = fetchn(S,'concat(session_date,"_",session_num,"_",mouse_name)->expRef');
        
        D = struct;
        for e = 1:length(expRefs)
            dd = loadData(expRefs{e});
            dd = structfun(@(x)(x(6:(end-14),:)),dd,'uni',0); %trim first 5 trials and last 15
            dd.laserRegion(dd.laserCoord(:,1)==0.5 & dd.laserCoord(:,2)==-2) = 'RightRSP';
            dd.laserRegion(dd.laserCoord(:,1)==-0.5 & dd.laserCoord(:,2)==-2) = 'LeftRSP';
            dd.laserRegion(dd.laserCoord(:,1)==-3.5 & dd.laserCoord(:,2)==-0.5) = 'LeftS1';
            dd.laserRegion(dd.laserCoord(:,1)==3.5 & dd.laserCoord(:,2)==-0.5) = 'RightS1';
            dd.laserRegion(dd.laserCoord(:,1)==1 & dd.laserCoord(:,2)==-0.5) = 'RightM1';
            dd.laserRegion(dd.laserCoord(:,1)==-1 & dd.laserCoord(:,2)==-0.5) = 'LeftM1';
            %             dd.sessionID = ones(length(dd.response),1)*e;
            
            if any(dd.laserType>0)
                D = addstruct(D,dd);
            end
            
            D = getrow(D,D.repeatNum==1);
        end
        
        varargout = {D};
        
    case 'InactivationTiming'
        %Plot effect of the laser at different times
        D = varargin{1};
        %         D = galvoPulseInactivation('getData');
        
        D.CL_gt_CR = D.stimulus(:,1) > D.stimulus(:,2);
        D.CR_gt_CL = D.stimulus(:,2) > D.stimulus(:,1);
        
        nonLaser = getrow(D,D.laserType==0 & (D.CR_gt_CL | D.CL_gt_CR) );
        nonLaser_perf = mean(nonLaser.feedbackType);
        
        inactivationRegions = {'VIS','M2','S1'};
        
        f = figure('color','w');        
        f.UserData.inactivationRegions = inactivationRegions;
        f.UserData.nonLaser_perf = nonLaser_perf;
        f.UserData.window_msec = 0.1;
        
        for laserPos = 1:length(inactivationRegions)
            e = getrow(D, (D.laserRegion == ['Left' inactivationRegions{laserPos}] & D.CR_gt_CL) | (D.laserRegion == ['Right' inactivationRegions{laserPos}] & D.CL_gt_CR) );
            
            times = e.laserOnset;
            %Add a bit of random noise to the times to allow smoothdata
            %function to work
            times = times + randn(length(times),1)/100000;
            
            fb = e.feedbackType;
            [times,sortIdx]=sort(times);
            fb = fb(sortIdx);
            
            f.UserData.times{laserPos} = times;
            f.UserData.feedback{laserPos} = fb;
        end
        
        f.WindowScrollWheelFcn = @(a,b)galvoPulseInactivation('InactivationTimingPlot_callback_scrollCount',a,b);
        galvoPulseInactivation('InactivationTimingPlot_callback_draw',f);
        
    case 'InactivationTimingPlot_callback_draw'
        figObj=varargin{1};
        delete(get(figObj,'children'));
        
                hold on;

        lx = line([-1 1],[1 1]*figObj.UserData.nonLaser_perf);
        set(lx,'Linestyle','--','Color',[0.5 0.5 0.5]);
                
        for laserPos = 1:length(figObj.UserData.inactivationRegions)
            times = figObj.UserData.times{laserPos};
            fb = figObj.UserData.feedback{laserPos};
% 
%             fb_smooth = smoothdata(fb,'movmean',figObj.UserData.window_msec,'SamplePoints',times);
%             plot(times,fb_smooth);
%             
%             
%             
            
            %Create sliding window
            tSteps = linspace(times(1),times(end),200)';
            fb_smoothed = nan(size(tSteps));
            fb_sm_ci = [];
            for t = 1:length(tSteps)
                idx = logical( WithinRanges(times, tSteps(t) + [-0.5 0.5]*figObj.UserData.window_msec) );
                fbI = fb(idx);
                [fb_smoothed(t), fb_sm_ci(t,:)] = binofit( sum(fbI), length(fbI));
            end
            
%             axL=plot(tSteps,fb_smoothed, '.-','markersize',20);
            axL=plot(tSteps,fb_smoothed, 'o-');
            %only highlight points which are significant
            sigIdx = ~( fb_sm_ci(:,1) < figObj.UserData.nonLaser_perf & figObj.UserData.nonLaser_perf < fb_sm_ci(:,2) );
            axL.MarkerIndices = find(sigIdx);

        end
        hold off;
        
        xlim([-0.3 0.3]);
        ylim([0.4 1]);
        
         title({'Inactivation contralateral to stimulus side',sprintf('Smoothing=%.0f msec',1000*figObj.UserData.window_msec)});
        xlabel('Time of laser relative to stimulus onset');
        ylabel('pCorrect');
        
               
    case 'InactivationTimingPlot_callback_scrollCount'
        figObj=varargin{1};
        scrollEvent=varargin{2};
        
        figObj.UserData.window_msec = figObj.UserData.window_msec + scrollEvent.VerticalScrollCount*scrollEvent.VerticalScrollAmount/1000;
        figObj.UserData.window_msec = max(figObj.UserData.window_msec, 0.0001);
        
        galvoPulseInactivation('InactivationTimingPlot_callback_draw',figObj)
end

end