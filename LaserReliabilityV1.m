subject = {'Eijkman'};
expRefs = dat.listExps(subject)';
expRefs = vertcat(expRefs{:});

expRefs = {
    '2016-01-27_1_Hopkins';
    '2016-01-29_1_Hopkins';
    '2016-02-01_1_Hopkins';
    '2016-02-01_2_Hopkins';
    '2016-02-01_3_Hopkins';
    '2016-02-01_4_Hopkins';
    '2016-02-02_1_Hopkins';
    '2016-02-03_1_Hopkins';
    '2016-02-05_1_Hopkins';
    '2016-02-05_2_Hopkins';
    '2016-02-05_3_Hopkins';
    '2016-02-05_4_Hopkins';
    '2016-02-05_5_Hopkins';
    '2016-02-05_6_Hopkins';
    '2016-02-09_1_Hopkins';
    '2016-02-09_2_Hopkins';
    '2016-02-11_1_Hopkins';
    '2016-02-12_1_Hopkins';
    '2016-02-15_1_Hopkins';
    '2016-02-16_1_Hopkins';
    '2016-02-17_1_Hopkins';
    '2016-02-18_1_Hopkins';
    '2016-02-19_1_Hopkins';
    '2016-02-19_2_Hopkins';
    '2016-02-22_1_Hopkins';
    '2016-02-23_1_Hopkins';
    '2016-02-24_1_Hopkins';
    '2016-02-26_1_Hopkins';
    '2016-02-29_1_Hopkins';
    '2016-02-29_2_Hopkins';
    '2016-02-29_3_Hopkins';
    '2016-03-01_1_Hopkins';
    '2016-03-01_2_Hopkins';
    };

% expRefs = {
%         '2016-01-27_1_Eijkman';
%     '2016-01-29_1_Eijkman';
%     '2016-02-01_1_Eijkman';
%     '2016-02-02_1_Eijkman';
%     '2016-02-03_1_Eijkman';
%     '2016-02-09_1_Eijkman';
%     '2016-02-11_1_Eijkman';
%     '2016-02-12_1_Eijkman';
%     '2016-02-15_1_Eijkman';
%     '2016-02-16_1_Eijkman';
%     '2016-02-16_2_Eijkman';
%     '2016-02-16_3_Eijkman';
%     '2016-02-17_1_Eijkman';
%     '2016-02-18_1_Eijkman';
%     '2016-02-19_1_Eijkman';
%     '2016-02-22_1_Eijkman';
%     '2016-02-24_1_Eijkman';
%     '2016-02-24_2_Eijkman';
%     };

model = 'C50-subset';
D=struct;
a=1;
figure;
for b=1:length(expRefs)
    try
        tr=dat.loadBlock(expRefs{b}).trial;
        c=[tr.condition];
        laser=[c.rewardOnStimulus]'; laser=(laser(:,2)>0);
        
        if sum(laser)==0 
            error('no laser trials in this session');
        end
        
        empty=arrayfun(@(i)isempty(tr(i).responseMadeID),1:length(tr))';
        
        E=struct;
        E.contrast_cond = [c(~empty).visCueContrast]';
        E.response = [tr(~empty).responseMadeID]';
        E.repeatNum = [c(~empty).repeatNum]';
        E.laserIdx = laser(~empty);
        
        subplot(3,4,a)
        gNolaser=GLM(getrow(E,E.laserIdx==0)).setModel(model).fit;
        gNolaser.plotFit;
        l=get(gca,'children');
        delete(l(4:6))
        set(l(1:3),'linestyle','--');
        
        
        hold on;
        glaser=GLM(getrow(E,E.laserIdx==1)).setModel(model).fit;
        glaser.plotFit;
        l=get(gca,'children');
        delete(l(4:6))
        hold off;
        
        drawnow;
        a=a+1;
    catch
    end
end