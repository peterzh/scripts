%% Get expRefs
expt = 'sparse_bilateral_2D';

expRefs = {
    'Spemann',...
    LaserOmnibusCheckSessions('Spemann',expt);
    
    'Murphy',...
    LaserOmnibusCheckSessions('Murphy',expt);
    
    'Morgan',...
    LaserOmnibusCheckSessions('Morgan',expt);
    
    'Whipple',...
    LaserOmnibusCheckSessions('Whipple',expt)
    };

numSubjects = size(expRefs,1);
disp('done');

%% Load timing data
clear l;
cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
        linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
    
for s = 1:numSubjects
    for b = 1:size(expRefs{s,2},1)
        block=dat.loadBlock(expRefs{s,2}{b});
        trial=block.trial;
        timeStart = [trial.stimulusCueStartedTime]';
        l = laserGLM(expRefs{s,2}{b});
        
        if length(timeStart) > size(l.data.response,1)
            timeStart(end)=[];
        end
        
        idx = l.data.laserIdx==3;
        timeStart = timeStart(idx)-1;
        
        vidFilename = ['\\zserver\Data\EyeCamera\' expRefs{s,1} '\' expRefs{s,2}{b}(1:10) '\' expRefs{s,2}{b}(12) '\face.mj2'];
        v = VideoReader(vidFilename);
        
        scaleFactor = 0.2;
        frames = 120;
        imgC = cell(frames,length(timeStart));
        for offset=1:frames
            for t = 1:length(timeStart)
                v.CurrentTime = timeStart(t)+offset/v.FrameRate;
                imgC(offset,t)={imresize(v.readFrame,scaleFactor)}; 
            end
            disp(offset);
        end
        
        numEle = length(timeStart);
        div = floor(sqrt(numEle));
        while mod(numEle,div)~=0
            numEle = numEle - 1;
        end
        
        while 1==1
            disp('Starting');
            for offset=1:frames
                image(cell2mat(reshape(imgC(offset,1:numEle),div,[])))
                pause(1/v.FrameRate);
            end
        end
        
       keyboard; 
    end
end

%% Load eye data and track pupil diameter
clear l;
    
for s = 1:numSubjects
    for b = 1:size(expRefs{s,2},1)
        block=dat.loadBlock(expRefs{s,2}{b});
        trial=block.trial;
        timeStart = [trial.stimulusCueStartedTime]';
        l = laserGLM(expRefs{s,2}{b});
        
        if length(timeStart) > size(l.data.response,1)
            timeStart(end)=[];
        end
        
        idx = l.data.laserIdx==3;
        timeStart = timeStart(idx)-1;
        
        vidFilename = ['\\zserver\Data\EyeCamera\' expRefs{s,1} '\' expRefs{s,2}{b}(1:10) '\' expRefs{s,2}{b}(12) '\eye.mj2'];
        v = VideoReader(vidFilename);
        
        scaleFactor = 0.2;
        frames = 120;
        imgC = cell(frames,length(timeStart));
        for offset=1:frames
            for t = 1:length(timeStart)
                v.CurrentTime = timeStart(t)+offset/v.FrameRate;
                imgC(offset,t)={imresize(v.readFrame,scaleFactor)}; 
            end
            disp(offset);
        end
        
        numEle = length(timeStart);
        div = floor(sqrt(numEle));
        while mod(numEle,div)~=0
            numEle = numEle - 1;
        end
        
        while 1==1
            disp('Starting');
            for offset=1:frames
                image(cell2mat(reshape(imgC(offset,1:numEle),div,[])))
                pause(1/v.FrameRate);
            end
        end
        
       keyboard; 
    end
end