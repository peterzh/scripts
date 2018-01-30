% Demo ALF file registration to alyx
%% Get instance
warning('MAKE SURE YOURE ON ALYX-DEV FOR THIS DEMO');
alyxInstance = alyx.loginWindow();

%% Create data formats
d = struct('name','npy',...
            'description','npy-formatted square numerical array',...
            'alf_filename', '*.*.npy',...
            'matlab_loader_function', 'readNPY',...
            'python_loader_function', 'numpy.load');
alyx.postData(alyxInstance, 'data-formats', d);

d = struct('name','mj2',...
            'description','mj2 movie',...
            'alf_filename', '*.*.mj2',...
            'matlab_loader_function', 'VideoReader',...
            'python_loader_function', 'ffmpeg');
alyx.postData(alyxInstance, 'data-formats', d);

d = struct('name','notData',...
            'description','file format for datasets with no associated files (eg parent datasets)',...
            'alf_filename', 'NA',...
            'matlab_loader_function', 'NA',...
            'python_loader_function', 'NA');
alyx.postData(alyxInstance, 'data-formats', d);

d = struct('name','mat',...
            'description','matlab native file format',...
            'alf_filename', '*.*.mat',...
            'matlab_loader_function', 'load',...
            'python_loader_function', 'scipy.io.loadmat');
alyx.postData(alyxInstance, 'data-formats', d);

%% Create parent datasetTypes
parents = {'cwFeedback';
    'cwGoCue';
    'cwResponse';
    'cwStimOn';
    'cwTrials';
    'eye';
    'face';
    'lickSignal';
    'licks';
    'passiveBeep';
    'passiveStimOn';
    'passiveTrials';
    'passiveValveClick';
    'passiveWhiteNoise';
    'sparseNoise';
    'spikes';
    'spontaneous';
    'wheel';
    'wheelMoves'};

for p = 1:length(parents)
    d = struct('name',parents{p},...
               'description',sprintf('Parent dataset %s',parents{p}) );
    alyx.postData(alyxInstance, 'dataset-types', d);

end

%% Create child dataset types
d = struct('name','timestamps','description','timestamps','alf_filename', '*.timestamps.*');
alyx.postData(alyxInstance, 'dataset-types', d);

d = struct('name','times','description','times of an event','alf_filename', '*.times.*');
alyx.postData(alyxInstance, 'dataset-types', d);

d = struct('name','intervals','description','start/stop intervals','alf_filename', '*.intervals.*');
alyx.postData(alyxInstance, 'dataset-types', d);

d = struct('name','eye.area','description','area of pupil','alf_filename', 'eye.area.npy');
alyx.postData(alyxInstance, 'dataset-types', d);

d = struct('name','eye.blink','description','area of pupil','alf_filename', 'eye.blink.npy');
alyx.postData(alyxInstance, 'dataset-types', d);

d = struct('name','eye.xPosition','description','x position of pupil','alf_filename', 'eye.xPosition.npy');
alyx.postData(alyxInstance, 'dataset-types', d);

d = struct('name','eye.yPosition','description','y position of pupil','alf_filename', 'eye.yPosition.npy');
alyx.postData(alyxInstance, 'dataset-types', d);

d = struct('name','eye.movie','description','movie of the eye','alf_filename', 'eye.movie.mj2');
alyx.postData(alyxInstance, 'dataset-types', d);

d = struct('name','lickSignal.trace','description','raw lick trace','alf_filename', 'lickSignal.trace.npy');
alyx.postData(alyxInstance, 'dataset-types', d);

d = struct('name','wheel.position','description','wheel position','alf_filename', 'wheel.position.npy');
alyx.postData(alyxInstance, 'dataset-types', d);

d = struct('name','wheel.velocity','description','wheel velocity','alf_filename', 'wheel.velocity.npy');
alyx.postData(alyxInstance, 'dataset-types', d);

d = struct('name','wheelMoves.type','description','classified type of movmeent (L,R,flinch,other)','alf_filename', 'wheelMoves.type.npy');
alyx.postData(alyxInstance, 'dataset-types', d);

d = struct('name','sparseNoise.positions','description','x and y positions at which stimuli were shown','alf_filename', 'sparseNoise.positions.npy');
alyx.postData(alyxInstance, 'dataset-types', d);

d = struct('name','cwFeedback.type','description','whether it is positive or negative','alf_filename', 'cwFeedback.type.npy');
alyx.postData(alyxInstance, 'dataset-types', d);

d = struct('name','cwResponse.choice','description','which choice was made, 1, 2, or 3','alf_filename', 'cwResponse.choice.npy');
alyx.postData(alyxInstance, 'dataset-types', d);

d = struct('name','cwStimOn.contrastLeft','description','contrast of left stimulus','alf_filename', 'cwStimOn.contrastLeft.npy');
alyx.postData(alyxInstance, 'dataset-types', d);

d = struct('name','cwStimOn.contrastRight','description','contrast of right stimulus','alf_filename', 'cwStimOn.contrastRight.npy');
alyx.postData(alyxInstance, 'dataset-types', d);

d = struct('name','cwTrials.inclTrials','description','boolean suggesting which trials to include/exclude in analysis','alf_filename', 'cwTrials.inclTrials.npy');
alyx.postData(alyxInstance, 'dataset-types', d);

d = struct('name','cwTrials.repNum','description','the repetition number of the trial, i.e. if it is the same condition as the previous trial(s)','alf_filename', 'cwTrials.repNum.npy');
alyx.postData(alyxInstance, 'dataset-types', d);


%% Register ALF directory
%First create session and subsession
alfDir = '\\zserver.cortexlab.net\Data\Subjects\Cori\2016-12-18\alf';
sessionURL = 'https://alyx-dev.cortexlab.net/sessions/5e15d67b-106e-4314-bd07-a7870a96cf0b';

alyx.registerALFtoAlyx(alfDir,sessionURL,alyxInstance);
