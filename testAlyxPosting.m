%% Alyx instance
alyxInstance = alyx.loginWindow();

%% create datarepo type 'fileserver'
d = struct;
d.name = 'fileserver';

wa = alyx.postData(alyxInstance, 'data-repository-type', d);


%% Create datarepo zserver
d = struct;
d.name = 'zserver';
d.path = '\\zserver.cortexlab.net\Data';
d.repository_type = 'fileserver';

wa = alyx.postData(alyxInstance, 'data-repository', d);

%% Create dataset types
d = struct('name','Block');
wa = alyx.postData(alyxInstance, 'dataset-types', d);

d = struct('name','Timeline');
wa = alyx.postData(alyxInstance, 'dataset-types', d);

d = struct('name','Parameters');
wa = alyx.postData(alyxInstance, 'dataset-types', d);

%% Create a session to link data to
d = struct;
d.subject = 'Nyx';
d.users = {'Peter','Nick'};
d.location = 'Ephys room 2.1.03';
d.procedures = "Behavior training/tasks";
d.narrative = 'nothing';
wa = alyx.postData(alyxInstance, 'sessions', d);


%% try registerFile function on test file
file = '\\zserver.cortexlab.net\Data\expInfo\Nyx\2017-06-09\1\2017-06-09_1_Nyx_Block.mat';
sessionUUID = 'd712f492-b8e2-44ce-a7c4-c976a0567597';
alyx.registerFile(sessionUUID,'Block',file,'zserver',alyxInstance);
