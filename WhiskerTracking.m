function [varargout] = WhiskerTracking(filename,fig)
%Procedure which takes in a whisking data file and calculates the average
%movement of the whiskers along the measurement axis.
%INPUT: 
    %-filename (e.g. '//zserver2/Data/EyeCamera/Fibiger/20150903/whisk1')
    %-figure flag (0/1)
%OUTPUT:
    %A struct containing fields:
    %-Timestep used for the analysis
    %-Vector containing the shifts in position
    
%the number of elements over which a cross-correlation comparison is made
dt=3;

% Load data 
disp('Loading data');
datFull = readWhisker(filename);

% Select subset of data, detrend it, and plot if fig == 1
disp('Detrending');
dat = -bsxfun(@minus,datFull,median(datFull,2));

% CROSSCORR METHOD
%Try cross-correlation method for tracing the change in position through time
%At each timestep, calculate the crosscorrelation with the following
%timestep and determine the lag which gives the highest correlation. That
%lag is taken to be the change in position over that time.
disp(['Calculating cross-correlations, using dt=' num2str(dt)]);
t_steps = 1:dt:(length(dat)-dt);
LAGS=nan(length(t_steps),1);

a=1;
for t = t_steps
    %Cross correlate adjacent timepoints
    [xcf,lags,~]=crosscorr(dat(:,t),dat(:,t+dt));
    
    %Find lag with maximum cross-correlation
    [~,idx]=max(xcf);
    
    %quality checking: if the lag is unrealistically high, discard
    LAGS(a)=lags(idx);
    a=a+1;
end

if fig==1
    figure;
    a=subplot(2,1,1);
    imagesc(dat);
    b=subplot(2,1,2);
    plot(t_steps,LAGS);
    set(gca,'YDir','reverse');
    ylabel(['Shifts along the measurement axis over dt=' num2str(dt)])
    linkaxes([a,b],'x');
end

W       = struct;
W.dt    = dt;
W.dP    = LAGS;

varargout = {W};

end