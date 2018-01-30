
%% Create geom.csv file, marking the electrode coords
ephysDir = 'Z:\Subjects\Cori\2016-12-18\ephys_V1';
cd(ephysDir);
mkdir('mountainsort');

coords = readNPY(fullfile(ephysDir,'sorting','channel_positions.npy'));
map = readNPY(fullfile(ephysDir,'sorting','channel_map.npy'));

writemda16i(map+1,fullfile(ephysDir,'mountainsort','channelSubset.mda'));
dlmwrite(fullfile(ephysDir,'mountainsort','geom.csv'),coords,',');

binName = dir(fullfile(ephysDir,'*.ap_CAR.bin'));

mdaCmd = ['mp-run-process pyms.extract_timeseries ',...
          '--channels_array=channelSubset.mda ',...
          '--timeseries=', binName.name,' ',...
          '--timeseries_out=raw.mda ',...
          '--timeseries_dtype=int16 ',...
          '--timeseries_num_channels=385\n'];
sortCmd = 'mlp-run mountainsort3.mlp sort --raw=raw.mda --geom=geom.csv --firings_out=firings.mda --_params=params.json --curate=true\n';

%Write bash file
f = fopen(fullfile(ephysDir,'mountainsort','runSort'),'w');
fprintf(f, mdaCmd);
fprintf(f, sortCmd);
fclose(f);

%Write params.json file
f = fopen(fullfile(ephysDir,'mountainsort','params.json'),'w');
fprintf(f, '{"samplerate":30000, "detect_sign":-1, "adjacency_radius":100}');
fclose(f);