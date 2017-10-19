%% Alyx instance
alyxInstance = alyx.loginWindow();

%% create datarepo type 'fileserver'
d = struct;
d.name = 'fileserver';

wa = alyx.postData(alyxInstance, 'data-repository-type', d);


%% Create datarepo zserver
d = struct;
d.name = 'zserver';
d.path = '\\zserver.cortexlab.net\Data\';
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
d.subject = subject;
%           d.users = {alyxInstance.username};
d.procedures = {'Behavior training/tasks'};
d.narrative = 'auto-generated session';
d.start_time = datestr(floor(now),31);
w = alyx.postData(alyxInstance, 'sessions', d);
wa = alyx.postData(alyxInstance, 'sessions', d);


%% try registerFile function on test file
file = '\\zserver.cortexlab.net\Data\expInfo\Nyx\2017-05-11\1\2017-05-11_1_Nyx_Block.mat';
alyx.registerFile('Nyx',[],'Block',file,'zserver',alyxInstance);
