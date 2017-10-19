%% example script with everything: part 1, finding movements

mouseName = 'Radnitz'; thisDate = '2017-01-08'; expNum = 3;

load(dat.expFilePath(mouseName, thisDate, expNum, 'block', 'master'));
rawpos = block.inputSensorPositions;
rawposT = block.inputSensorPositionTimes;

% resample to an even 1kHz sampling rate
Fs = 1000; 
t = rawposT(1):1/Fs:rawposT(end);
pos = interp1(rawposT, rawpos, t);

params.makePlots = true;
[moveOnsets, moveOffsets] = wheel.findWheelMoves3(pos, t, Fs, params);

%% example script with everything: part 2, classifying movements

tr = block.trial;
tr = tr(1:block.numCompletedTrials);
respTime = [tr.responseMadeTime];
intStartTime = [tr.interactiveStartedTime];
resp = [tr.responseMadeID];

hasTurn = resp==1|resp==2; 
resp = resp(hasTurn);
intStartTime = intStartTime(hasTurn); 
respTime = respTime(hasTurn);

moveType = wheel.classifyWheelMoves(t, pos, moveOnsets, moveOffsets, intStartTime, respTime, resp);

clear dm; dm.moveOnsets = moveOnsets; dm.moveOffsets = moveOffsets; dm.moveType = moveType;
plotWheel(t, pos, dm);