
binFile = 'E:\raw_data\Cori_start.bin';
info = dir(binFile);
sample_rate = 30000;
numChannels = 385;
encoding = 'int16';

totalNumSamples = info.bytes/2;
samplesPerChannel = totalNumSamples/numChannels;
outputFile = 'E:\raw_data\header.mda';


%% Create header file
fout = fopen(outputFile,'w');
fwrite(fout,-4,'int32'); %-4 for int16
fwrite(fout,2,'int32'); %2 for number of bytes in each entry (?)
fwrite(fout,2,'int32'); %2 for number of dimensions in raw data (num channels by num timepoints)
fwrite(fout,numChannels,'int32'); 
fwrite(fout,samplesPerChannel,'int32'); 
fclose(fout);

% Now use 'cat header.mda Radnitz..... > full.mda' on unix