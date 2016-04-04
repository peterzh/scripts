%Script for analysing laser manipulator drift
l = laserGLM('2016-04-01_2_test');
v = VideoReader('\\zserver2\Data\EyeCamera\TEST\2016-04-01\5\face.mj2');
img = imrotate(v.readFrame,180);

refImg = ones(100,100);
refImg(50,50)=0; 
refImg(20,[45 55])=0;
refImg(50,[15 85])=0;
refImg(90,[15 85])=0;

[mP,fP]=cpselect(refImg,img,'Wait',true);
tform = fitgeotrans(fP,mP,'projective');
imgW = imwarp(img,tform);
imshowpair(refImg,imgW,'montage');

%% Warp image and find peaks
v = VideoReader('\\zserver2\Data\EyeCamera\TEST\2016-04-01\5\face.mj2');
f=1; i=1; laserPos=[];
skip=10;
stack = nan(size(imgW,1),size(imgW,2),round(v.Duration*v.FrameRate/skip));
while hasFrame(v)
    imgW=imwarp(imrotate(readFrame(v),180),tform);
    
    if f>skip
        p=FastPeakFind(imgW,250);
        if ~isempty(p)
            laserPos(:,i)=p(1:2);
        else
            laserPos(:,i)=[nan;nan];
        end
%         imshow(imgW);
        scatter(laserPos(1,:),-laserPos(2,:),'.');
        drawnow;
        stack(:,:,i)=imgW;
        f=0; i=i+1;
    end
    f=f+1;
end
imagesc(max(stack,[],3));
disp('done!');
% scatter(laserPos(1,:),-laserPos(2,:),'o');
