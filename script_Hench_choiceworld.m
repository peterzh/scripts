

%% load imaging data


mouseName = 'Hench';
thisDate = '2017-06-12';
expNum = 2; 
nSV = 500;
Fs = 35;

expPath = dat.expPath(mouseName, thisDate, expNum, 'main', 'master');
expRoot = fileparts(expPath);


%%


[U, V, t, mimg] = quickHemoCorrect(expPath, nSV);

%% if already corrected:

[U, V, t, mimg] = quickLoadUVt(expPath, nSV);


%% correlation map
pixelCorrelationViewerSVD(U, V)

%% a movie of the session
quickMovieWithVids(mouseName, thisDate, expNum)

%% load timeline and block

load(dat.expFilePath(mouseName, thisDate, expNum, 'Timeline', 'master'));
load(dat.expFilePath(mouseName, thisDate, expNum, 'block', 'master'));

%% other way to do behavioral variables
% addpath('F:\Dropbox\ucl\code\behavioralTask\burgessTask');
% computeAndSaveBehavioralVars(mouseName, thisDate, expNum, expNum)

%% Some things from the block

tr = block.trial; tr = tr(1:block.numCompletedTrials);
cond = [tr.condition];
stimOn = [tr.stimulusCueStartedTime];
vcc = [cond.visCueContrast];
contrastLeft = vcc(1,:);
contrastRight = vcc(2,:);
choice = [tr.responseMadeID];
goCue = [tr.interactiveStartedTime];
responseTime = [tr.responseMadeTime];
reactionTime = [tr.responseMadeTime]-goCue;
feedback = [tr.feedbackType];
repNum = [cond.repeatNum];
trialStarts = [tr.trialStartedTime];
trialEnds = [tr.trialEndedTime];

[contrastCondDefs, ia, contrastCondsRaw] = unique(vcc', 'rows');

sw = block.stimWindowUpdateTimes;

%% get alignment between block and timeline
pd = Timeline.rawDAQData(:,2);
tt = Timeline.rawDAQTimestamps;

pdFlips = schmittTimes(tt,pd, [3 5]);

figure; 
plot(sw, ones(size(sw)), '.');
hold on; 
plot(pdFlips, ones(size(pdFlips))+1, '.');
ylim([-3 6])

blockToTL = makeCorrection(pdFlips(2:end-1), sw, false);


%% quick plot of wheel traces in this session

st = blockToTL(stimOn);

tt = Timeline.rawDAQTimestamps;
wh = Timeline.rawDAQData(:, strcmp({Timeline.hw.inputs.name}, 'rotaryEncoder'));

wh = wh(tt<max(st)+5); tt = tt(tt<max(st)+5);

wh = wheel.correctCounterDiscont(wh);
Fs = Timeline.hw.daqSampleRate;
whv = wheel.computeVelocity2(wh, 0.01, Fs);


win = -0.2:0.002:0.5;

winSamps = bsxfun(@plus, st, win);

whTrig = interp1(tt, wh, winSamps);
whTrig = bsxfun(@minus, whTrig, whTrig(:,find(win<0,1,'last')));

figure; 
% colors = colorsLeftRight('all');
colors = [1 0 0; 0 0 1; 0 0 0];
for r = 1:3
    plot(win, whTrig(feedback==1 & repNum==1 & choice==r,:), 'Color', colors(r,:))
    hold on;
end

%% detect movement onsets

% clear params
% params.makePlots = true;
% newFs = 500;
% newt = 0:1/newFs:tt(end);
% newwh = interp1(tt,wh,newt);
% [mon, moff] = wheel.findWheelMoves3(newwh(:), newt(:), newFs, params);

Fs = 1000;
[moveOnsets, moveOffsets] = wheel.findWheelMoves3(wh, tt, Fs, []);

resp = choice; hasTurn = resp==1|resp==2; resp = resp(hasTurn);
intStartTime = blockToTL(goCue(hasTurn)); respTime = blockToTL(responseTime(hasTurn));
moveType = wheel.classifyWheelMoves(tt, wh, moveOnsets, moveOffsets, intStartTime, respTime, resp);
clear dm; dm.moveOnsets = moveOnsets; dm.moveOffsets = moveOffsets; dm.moveType = moveType;
plotWheel(tt, wh, dm);

%% pick move onsets
son = blockToTL(stimOn);
mon = arrayfun(@(x)moveOnsets(find(moveOnsets>x,1)), son);
rTime = blockToTL(responseTime);
countMove = arrayfun(@(x)mon(x)<rTime(x), 1:numel(mon));



%% look at movement onsets

win = -0.2:0.002:0.5;

winSamps = bsxfun(@plus, mon, win);

whTrig = interp1(tt, wh, winSamps);
whTrig = bsxfun(@minus, whTrig, whTrig(:,find(win<0,1,'last')));

figure; 
% colors = colorsLeftRight('all');
colors = [1 0 0; 0 0 1; 0 0 0];
for r = 1:3
    plot(win, whTrig(feedback==1 & repNum==1 & choice==r & countMove,:), 'Color', colors(r,:))
    hold on;
end

%% rt distribution
rt = [mon-son]';

figure; 
for r = 1:3
    subplot(3,1,r);
    hist(rt(feedback==1 & repNum==1 & choice==r & countMove), 0:0.01:0.5);
    xlim([0 0.5]);
end

% answer is: between about 150 and 225ms in this session! Could be shifted
% back by 50ms with a less careful detection of onset time. 

%% look at an event-triggered activity
% stim onset, just for single stimulus trials
inclByRT = rt>0.125 & (rt<0.5 | choice==3);
inclTrials = contrastRight==0 | contrastLeft==0 & inclByRT & repNum==1 & feedback==1;
eventTimes = blockToTL(stimOn(inclTrials));
eventLabels = contrastRight(inclTrials)-contrastLeft(inclTrials);

%%
pixelTuningCurveViewerSVD(U, V, t, eventTimes, eventLabels, [-0.5 2.5])
set(gcf,'name','Stimulus aligned (U and V)');

%% dff
[newU, newV] = dffFromSVD(U, V, mimg);

%% derivative?
dVold = diff(newV, [], 2); 
dtold = t(1:end-1)+mean(diff(t))/2;

%% another way to do derivative: with differentiation filter
%https://uk.mathworks.com/help/signal/ug/take-derivatives-of-a-signal.html
Nf = 50; 
Fpass = 8.5; 
Fstop = 10;
Fs = 1/mean(diff(t));

d = designfilt('differentiatorfir','FilterOrder',Nf, ...
    'PassbandFrequency',Fpass,'StopbandFrequency',Fstop, ...
    'SampleRate',Fs);

% fvtool(d,'MagnitudeDisplay','zero-phase','Fs',Fs)

dV = filter(d,newV')'; %*Fs;

delay = mean(grpdelay(d));
dt = t(1:end-delay);

dV(:,1:delay) = [];

dt(1:delay) = [];
dV(:,1:delay) = [];


%%
newUth = newU;
mask = double(mimg>6000);
for sv = 1:size(newUth,3)
    newUth(:,:,sv) = newUth(:,:,sv).*mask;
end

%%

% pixelTuningCurveViewerSVD(U, dV, dt, eventTimes, sign(eventLabels), [-0.5 1.5])
% pixelTuningCurveViewerSVD(newUth, dV, dt, eventTimes, sign(eventLabels), [-0.3 0.8])
pixelTuningCurveViewerSVD(newUth, dV, dt, eventTimes, eventLabels, [-0.3 0.8])
set(gcf,'name','Stimulus aligned (dF/F and derivative signal)');

%% for movement onset
inclTr = feedback==1 & repNum==1 & countMove;

pixelTuningCurveViewerSVD(newU, dV, dt, mon(inclTr), choice(inclTr), [-1.5 1])
set(gcf,'name','Movement aligned (dF/F and derivative signal)');

%% pull out a frame

epoch = 'ms270';

ch = get(gca, 'Children');
cax = caxis();
im = ch(end);
f = figure; imagesc(im.CData); set(f, 'Position', [61 271 685 603]);
colormap(colormap_blueblackred)
axis image; axis off; caxis(cax);
% addAllenCtxOutlines(bregma, lambda, 'w')
% print(f, ['widefield' epoch], '-dpng', '-r200')

% ch = get(gca, 'Children');
% cax = caxis();
% im = ch(end);
f = figure; imagesc(im.CData); set(f, 'Position', [33.0000  491.2857  677.7143  602.8571]);
colormap(colormap_blueblackred)
axis image; axis off; caxis(cax);
addAllenCtxOutlines(bregma, lambda, 'w')
% print(f, ['widefield' epoch '_withOutlines'], '-dpng', '-r200')


%% contour of the activity
brainIm = im.CData;

thresh = 120;

bimt = brainIm>thresh;
se = strel('disk',25);
bimt = imerode(imdilate(bimt, se), se);


%%
figure; hold on;


[c,h] = contour(double(bimt), [0.5 0.5]); delete(h);
ii = 1;
while ii<size(c,2)
    n = c(2,ii);
    x = c(1,ii+1:ii+n);
    y = c(2,ii+1:ii+n);
    x = smooth(x,50,'loess');
    y = smooth(y,50,'loess');
    x = [x;x(1)];
    y = [y;y(1)];    
    
    plot(x,y, 'LineWidth', 2.0, 'Color', 'r');
    fill(x,y, 'r', 'FaceAlpha', 0.5);
    ii = ii+n+1;

end
addAllenCtxOutlines(bregma, lambda, 'k')
set(gca, 'YDir', 'reverse');
axis image; axis off;

%% add inactivation coords

apDir = lambda-bregma; apDir = apDir./norm(apDir);

pixSize = 0.0217; % mm/pix. This is for PCO edge 5.5 with 0.6x mag (as kilotrode)

ccfbregma = allenCCFbregma()/100/pixSize;

coords = standardBrainCoords;

hold on; 
lineColor = 'b';
    
    % these are in 10um voxel coordinates, so first convert to mm, then to
    % pixels 
    cx = coords(1,:)'/pixSize;
    cy = -coords(2,:)'/pixSize;    
    
    T = affineMat.rot(-atan(apDir(2)/apDir(1)));
    
    newc = T*[cx cy ones(size(cx))]';
    cx = newc(1,:)'; cy = newc(2,:)';
    
    cx = cx+bregma(2); cy = cy+bregma(1);
    
    plot(cx,cy, 'o', 'Color', lineColor, 'MarkerSize', 10, 'MarkerFaceColor', lineColor);
    hold on;


%% next to do:
% - try hemocorrect!
% - pull out traces (below)
% - check RTs of this session
% - check which trials are included?
% - extract contour of the S1/M1 area and overlay with cortex outlines and
% inactivation locations


%% more detailed stimulus onset analysis

% 1. try upsampling to look at time course better, first for a single pixel
% and single contrast condition. 
% 2. Then for a selection of pixels. 

Fs = 500;
calcWin = [-0.3 0.8];
winSamps = calcWin(1):1/Fs:calcWin(2);

[avgPeriEventV, winSamps, periEventV, sortedLabels] = ...
    eventLockedAvgSVD(newU, dV, dt, eventTimes, sign(eventLabels), calcWin);

%%
px = [349 249];

% periEventV is nTrials x nSv x nTimepoints
inclTr = find(sortedLabels==-1); nTr = numel(inclTr);
thisU = squeeze(U(px(1), px(2), :));
thisTraces = zeros(nTr, numel(winSamps));
for c = 1:nTr
    thisTraces(c,:) = thisU'*squeeze(periEventV(inclTr(c),:,:));
end

figure; 
subplot(2,1,1); 
imagesc(winSamps, [], thisTraces)
box off; ylabel('trial number');
subplot(2,1,2);
plot(winSamps, median(thisTraces))
box off; xlabel('time from visual stimulus onset (s)'); ylabel('dF, arb units')

%% a few pixels

clear px; p = 1;
% px(p).name = 'VISp'; px(p).xy = [414 184]; p = p+1;
% px(p).name = 'VISpm'; px(p).xy = [399 205]; p = p+1;
% px(p).name = 'VISam'; px(p).xy = [331 185]; p = p+1;
% px(p).name = 'MOs'; px(p).xy = [139 223];p = p+1;
% px(p).name = 'MOp'; px(p).xy = [228 207];p = p+1;
px(p).name = 'VISp'; px(p).xy = [410 364]; p = p+1;
px(p).name = 'VISal'; px(p).xy = [342 450]; p = p+1;
% px(p).name = 'VISpm'; px(p).xy = [353 359]; p = p+1;
% px(p).name = 'VISam'; px(p).xy = [320 358]; p = p+1;
px(p).name = 'MOs'; px(p).xy = [131 292];p = p+1;
px(p).name = 'MOp'; px(p).xy = [209 333];p = p+1;


inclTr = find(sortedLabels==-1); nTr = numel(inclTr);


figure;
for p =1:numel(px)
    thisU = squeeze(newUth(px(p).xy(1), px(p).xy(2), :));
    thisTraces = zeros(nTr, numel(winSamps));
    for c = 1:nTr
        thisTraces(c,:) = thisU'*squeeze(periEventV(inclTr(c),:,:));
    end
    
    % rectify % not much effect, but doesn't hurt
    thisTraces(thisTraces<0) = 0;
    
    % try subtracting baselines? % not much effect, but doesn't hurt
    thisTraces = bsxfun(@minus, thisTraces, mean(thisTraces(:,winSamps<0),2));
        
    hold on; 
    plot(winSamps, mean(thisTraces)./max(mean(thisTraces)))
%     plot(winSamps, mean(thisTraces))
end
box off; xlabel('time from visual stimulus onset (s)'); ylabel('dF/dt, arb units')
xlim([-0.1 0.4])
hold on; 
xl = xlim(); plot(xl, [0 0], 'k--');
yl = ylim(); plot([0 0],yl,  'k--');
legend({px.name})
% makepretty;

%% plot where the points are
figure; 
imagesc(mimg); hold on; colormap gray; axis image; axis off;
for p =1:numel(px)
    h = plot(px(p).xy(2), px(p).xy(1), 'o');
    set(h, 'MarkerFaceColor', get(h, 'Color'));
end
addAllenCtxOutlines(bregma, lambda, 'w')

%%
% Hm, not so compelling in Hench - more movement artifacts. 
% 1. Try hemo-correction (copy dataset to zubjects)
% 2. Check out reichstein, maybe there's a good session. 

%% try showing mean image with atlas outlines overlaid

figure; 
imagesc(mimg)
axis image
axis off
colormap gray
hold on; 

bregma = [218 261];
lambda = [400 271];
apDir = lambda-bregma; apDir = apDir./norm(apDir);

pixSize = 0.0217; % mm/pix

ccfbregma = allenCCFbregma()/100/pixSize;

load('J:\allen\ctxOutlines.mat');

% allPos = standardBrainCoords(); % from optogenetic experiments so far
% allPos = expandedBrainCoords(); % filling this implant better

% plot(allPos(1,:)*100+bregma(3), -allPos(2,:)*100+bregma(1), 'bo', 'MarkerFaceColor','b');
% plot(bregma(2), bregma(1), 'ro', 'MarkerFaceColor', 'r'); % bregma
% plot(bregma(2)+4/pixSize*apDir(2), bregma(1)+4/pixSize*apDir(1), 'ro', 'MarkerFaceColor', 'r'); % lambda



hold on; 
for q = 1:numel(coords) % coords is from ctxOutlines.mat
    
    % these are in 10um voxel coordinates, so first convert to mm, then to
    % pixels 
    cx = coords(q).x/100/pixSize;
    cy = coords(q).y/100/pixSize;
    
    % to do this transformation, first subtract bregma to zero, then
    % rotate, then add back the other bregma
    
    cx = cx-ccfbregma(3); cy = cy-ccfbregma(1);
    
%     T1 = affineMat.trans(bregma(2)-ccfbregma(3),bregma(1)-ccfbregma(1));
    T2 = affineMat.rot(-atan(apDir(2)/apDir(1)));
%     T2 = affineMat.rot(0);
    
    newc = T2*[cx cy ones(size(cx))]';
    cx = newc(1,:)'; cy = newc(2,:)';
    
    cx = cx+bregma(2); cy = cy+bregma(1);
    
    plot(cx,cy, 'LineWidth', 1.0, 'Color', 'w');
    hold on;
end


load('J:\allen\wellEdges.mat');
offsetMM = [0.2 10.5]; % AP, LR
wellEdges = wellEdges/pixSize;
T2 = affineMat.rot(pi/2);
newc = T2*[wellEdges ones(size(wellEdges,1),1)]';
wx = newc(1,:)'; wy = newc(2,:)';
wx = wx+offsetMM(2)/pixSize; wy = wy+offsetMM(1)/pixSize;
hold on; 
plot(wx, wy, 'g');

%% try to find latency to peak for each pixel

% reconstruct each frame in order and determine whether the value in each
% pixel is bigger than a running max

win = [0 0.2];

maxVal = ones(size(mimg))*-Inf;
maxInd = ones(size(mimg));

% inclTr = find(sortedLabels==-1); nTr = numel(inclTr);
inclCond = 1;

Fs = 500;
calcWin = [-0.3 0.8];
winSamps = calcWin(1):1/Fs:calcWin(2);

[avgPeriEventV, winSamps, periEventV, sortedLabels] = ...
    eventLockedAvgSVD(newU, dV, dt, eventTimes, sign(eventLabels), winSamps);

inclSamps = find(winSamps>win(1) & winSamps<win(2));
tVals = winSamps(inclSamps);
gw = gausswin(10)./sum(gausswin(10));
useV = squeeze(avgPeriEventV(inclCond,:,:));
useVsm = conv2(1, gw, useV, 'same');
figure; 
for ss = 1:numel(inclSamps)
    fr = svdFrameReconstruct(newUth, useVsm(:,inclSamps(ss)));
    
    fr = conv2(gw, gw, fr, 'same');
    
    maxInd(fr>maxVal) = tVals(ss);
    maxVal(fr>maxVal) = fr(fr>maxVal);
    
    subplot(1,3,1);
    if ss==1
        im1 = imagesc(fr); caxis([-1 1]*3e-3); axis image; axis off;
        colormap(hot);
    else
        set(im1, 'CData', fr);
    end
    subplot(1,3,2);
    if ss==1
        im3 = imagesc(maxVal); caxis([-1 1]*3e-3); axis image; axis off;
%         colormap(colormap_blueblackred);
    else
        set(im3, 'CData', maxVal);
    end
    subplot(1,3,3);
    
    if ss==1
        im2 = imagesc(maxInd); caxis([0.1 0.2]); axis image; axis off;
    else
        set(im2, 'CData', maxInd);
    end
    drawnow;
    
end


% hm, can't really get this to look nice. :/



