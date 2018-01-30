%% Light baffle sessions

% % Cone sat ontop  %1.5mW 1.5sec
names = {'Nyx','Vin','Keynes'};
expRefs = {{'2017-10-26_2_Nyx'; '2017-10-27_4_Nyx'; '2017-11-01_3_Nyx'; '2017-11-03_1_Nyx'},...
           {'2017-10-26_1_Vin'; '2017-10-27_1_Vin'; '2017-10-30_2_Vin'; '2017-11-01_1_Vin'; '2017-11-03_1_Vin';},...
           {'2017-10-26_1_Keynes'; '2017-10-27_1_Keynes'; '2017-10-30_1_Keynes'; '2017-11-01_1_Keynes'; '2017-11-03_1_Keynes'}};
%        
%     

%Session with inactivation outside the brain but no light baffle?




% % Implant %1.5mW 1.5sec
names = {'Nyx','Vin','Keynes','Heinz'};
expRefs = {{'2017-11-23_1_Nyx';'2017-11-27_1_Nyx';'2017-11-28_1_Nyx';'2017-11-29_1_Nyx';'2017-11-30_1_Nyx'},...
           {'2017-11-23_1_Vin';'2017-11-27_1_Vin';'2017-11-28_4_Vin';'2017-11-29_1_Vin';'2017-11-30_1_Vin'},...
           {'2017-11-23_1_Keynes';'2017-11-27_1_Keynes';'2017-11-28_2_Keynes';'2017-11-29_1_Keynes';'2017-11-30_1_Keynes'},...
           {'2017-11-27_2_Heinz';'2017-11-28_1_Heinz';'2017-11-29_1_Heinz';'2017-11-30_1_Heinz'}};
%      
%Full light isolation %1.5mW 1.5sec
names = {'Nyx','Vin','Keynes','Heinz'};
expRefs = {{'2017-12-02_1_Nyx';'2017-12-03_1_Nyx';'2017-12-04_1_Nyx'},...
           {'2017-12-02_1_Vin';'2017-12-03_1_Vin';'2017-12-04_1_Vin'},...
           {'2017-12-02_1_Keynes';'2017-12-03_1_Keynes';'2017-12-04_1_Keynes'},...
           {'2017-12-02_1_Heinz';'2017-12-03_1_Heinz';'2017-12-04_1_Heinz'}};

%Full light isolation with brief pulses
% names = {'Nyx','Vin','Keynes','Heinz'};
tab = readtable('C:\Users\Peter\Desktop\LocalSessionData\SESSION_TABLE.csv');
eRefs = tab.expRef(contains(tab.notes,'Light isolation with brief pulse. Exclude from all other analyses!'));
names = {'all'};
expRefs = {eRefs};


o=omnibusLaserGLM(expRefs,names);
o=o.fit('BO_NC50','bias+sens+^c50');
o.plotFit_Laser(1);


%% Old
o=omnibusLaserGLM('galvo_unilateral_2D_tdelayPulse',{'Nyx','Vin'});