% Defines some configurations for the pulseInactivation analysis
rawDataDir = '\\zserver.cortexlab.net\Data2\Subjects';
rootDir = fullfile( fileparts(mfilename('fullpath')), '..' );
preprocDir = fullfile(rootDir,'preproc');
figuresDir = fullfile(rootDir,'figures');

subjects = {'Vin','Nyx','Keynes'};

makePlots = false;