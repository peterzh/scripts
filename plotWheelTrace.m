%Plot wheel trace of behaviour
 'Eijkman'     '2015-09-15'    1
[~,meta] = loadData('2018-05-01_1_Ochoa');
bl = load(meta.blockFile);

times = bl.block.inputs.wheelTimes';
wheel = bl.block.inputs.wheelValues';
eventTimes = bl.block.events.stimulusOnTimes';

%Wheel position at each eventTime
wheelBaseline = interp1(times, wheel, eventTimes);


%Wheel position at relative windows
timeSteps = linspace(-1,1,100);
wheelAligned = interp1(times, wheel, eventTimes + timeSteps);

%subtract wheel position at each event t = 0
wheelAligned = wheelAligned - wheelBaseline;


Lchoices = bl.block.events.responseValues == -1;
Rchoices = bl.block.events.responseValues == 1;
NGchoices = bl.block.events.responseValues == 0;


figure;
axis; hold on;
plot(timeSteps, wheelAligned(Lchoices,:), 'b-');
plot(timeSteps, wheelAligned(Rchoices,:), 'r-');
plot(timeSteps, wheelAligned(NGchoices,:), 'k-');

hx=errorbar( timeSteps, mean(wheelAligned), std(wheelAligned));
hx.Color = [1 1 1];

