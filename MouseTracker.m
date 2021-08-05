% Used for track the center of mass of the mouse or one specific color.

clear
close all

%% Read the video file
[fname pname]=uigetfile('*.*','Get the video file');
cd(pname)
currentVideo=VideoReader(fname);

%% Defination of parameters
tracking_color = 'black'; % Switch between 'black', 'red', 'green', 'blue'
% 'black' refers to the center of mass of the mouse
trimVideoLeft=340; % Number of pixels to be trimmed from the left
trimVideoRight=300; % Number of pixels to be trimmed from the right
trimVideoWidth=currentVideo.Width-(trimVideoLeft+trimVideoRight);
trimVideoTop=10; % Number of pixels to be trimmed from the top
trimVideoBottom=10; % Number of pixels to be trimmed from the bottom
trimVideoHeight=currentVideo.Height-(trimVideoTop+trimVideoBottom);
bining = 6; % Frame bining factor for processing the video

currentVideo.CurrentTime=5; % The starting time of video processing in s
i=0;
j=0;
vidObj = VideoWriter('test.avi','Motion JPEG AVI');
vidObj.FrameRate=25; % Hz
open(vidObj);

%% Tracking, plotting and writing video
while hasFrame(currentVideo)&&i<=10000
    currentFrame=readFrame(currentVideo);
    j=j+1;
    if mod(j,bining)==1 
        i=i+1 % Index indicator
        
        imageRed=rgbImage(:,:,1);                           % Red channel
        imageGreen=rgbImage(:,:,2);                         % Green channel
        imageBlue=rgbImage(:,:,3);                          % Blue channel
        grayImage=im2double(rgb2gray(currentFrame));        % Gray scale
        imageRedRatio=imageRed.*imageRed./imageGreen./imageBlue;
        imageGreenRatio=imageGreen.*imageGreen./imageBlue./imageRed-imageRed-imageBlue;
        imageBlueRatio=imageBlue.*imageBlue./imageGreen./imageRed-imageRed-imageRed;

        switch tracking_color
            case black
                untrimImage = grayImage;
            case red
                untrimImage = imageRedRatio;
            case green
                untrimImage = imageGreenRatio;
            case blue
                untrimImage = imageBlueRatio;
        end
        
        trimImage=imcrop(untrimImage, [trimVideoLeft trimVideoTop trimVideoWidth trimVideoHeight]);
        inverseImage=1-trimImage(:,:)-.65;
        rectifiedImage=inverseImage.*(inverseImage>0);
        spherifiedImage=imresize(imresize(rectifiedImage,[trimVideoHeight,trimVideoWidth]/20 ),[trimVideoHeight,trimVideoWidth]);
        [r, c]=size(spherifiedImage);
        biggerImage=zeros(r+400,c+400);
        biggerImage(201:200+r,201:200+c)=spherifiedImage;
        filteredImage=bpass(biggerImage,1,100);
        croppedImage=filteredImage(201:200+r,201:200+c);
        pk=pkfnd(filteredImage,max(max(filteredImage))*0.9,150); % May need to tune the threshold for different videos
        
        objectLocations(i,1)=pk(1,1)-200;
        objectLocations(i,1)=objectLocations(i,1)+trimVideoLeft;
        objectLocations(i,2)=pk(1,2)-200;
        objectLocations(i,2)=objectLocations(i,2)+trimVideoTop;
        
        if i~=1
            xdis(i)=objectLocations(i,1)-objectLocations(i-1,1);
            ydis(i)=objectLocations(i,2)-objectLocations(i-1,2);
            dis(i)=(xdis(i).^2+ydis(i).^2).^(1/2);
        else
            dis=0;
        end
         
        if dis(i)>400
            objectLocations(i,1)=objectLocations(i-1,1);
            objectLocations(i,2)=objectLocations(i-1,2);
        end 
                
        time=i/1;
        time=num2str(time);
        time=strcat(time,' s');

        fig=figure;
        imshow(currentFrame)
        caxis auto 
        hold on
        plot(objectLocations(:,1),objectLocations(:,2),'r.-','MarkerSize',20)
        pause(0.01)

        NIRMov(i)=getframe(fig);
        close(fig)      
        writeVideo(vidObj,NIRMov(i));
    
        save trajectory.dat objectLocations -ascii
    end
end

% Close the file.
close(vidObj);
