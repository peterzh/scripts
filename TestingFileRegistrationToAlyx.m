alyxInstance = alyx.loginWindow;

%% Create basic prereqs
d = struct('name','fileserver');
alyx.postData(alyxInstance, 'data-repository-type', d);

d = struct('name','zserver','path','\\zserver.cortexlab.net\');
alyx.postData(alyxInstance, 'data-repository', d);

d = struct('name','block');
alyx.postData(alyxInstance, 'dataset-types', d);


%% Create session and subsession
d = struct;
d.subject = 'Vin';
d.procedures = {'Behavior training/tasks'};
d.narrative = 'auto-generated session';
d.start_time = datestr(now,31);
d.type = 'Base';
session = alyx.postData(alyxInstance, 'sessions', d);

d = struct;
d.subject = 'Vin';
d.procedures = {'Behavior training/tasks'};
d.narrative = 'auto-generated session';
d.start_time = datestr(now,31);
d.type = 'Experiment';
d.parent_session = session.url;
d.number = 3;
subsession = alyx.postData(alyxInstance, 'sessions', d);

%% Register a file for that latest subclass
filePath = '\\zserver.cortexlab.net\Data2\Subjects\Vin\2017-09-01\1\2017-09-01_1_Vin_Block.mat';
alyx.registerFile('Vin',[],'Block',filePath,'zserver',alyxInstance)


%% Find today's base session
subject = 'Vin';
thisDate = alyx.datestr(now);

sessions = alyx.getData(alyxInstance, ['sessions?type=Base&subject=' subject]);  
latest_base = sessions{end};
alyx.datestr(alyx.datenum(latest_base.start_time))