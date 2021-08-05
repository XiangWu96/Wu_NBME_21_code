% Written by Guosong Hong on April 8, 2015
% Modified and simplified by Guosong Hong on July 7, 2015
% New clustering methods added by Guosong Hong on July 8, 2015

% Initialization
clear;
%close all;
set(0,'DefaultFigureWindowStyle','docked')

%%%%%%%%%%%%%%%%%%%%%%%%%PARAMETER INITIALIZATIONS%%%%%%%%%%%%%%%%%%%%%%%%%
% Change PC indices here
% Number 1 represents the first principal component, number 2 the second
% principal component, so on and so forth.;o,
PCIndex1=1;
PCIndex2=2;

% Specify PCA cluster method
% Three clustering methods are available, being 'sharp', 'blur' and
% 'overlap'.
PCA_method='sharp';

% Save option
saveOrNot='yes';

% Specify threshold for clustering method of 'blur'
testorThreshold_blur = 0.7;

% Change number of clusters here
% Maximum number is 7 based on the following 'colorOptions'; however, one
% can always generate almost unlimited number of color options using RGB
% array.
numClus=2;

% Change color codes here
% k - black
% r - red
% b - blue
% g - green
% y - yellow
% c - cyan
% m - magenta
colorOptions=['r' 'g' 'b' 'k' 'y' 'c' 'm'];
toSaveNameOptions={'red' 'green' 'blue' 'black' 'yellow' 'cyan' 'magenta'};

% Change timestamp here
timeStamp='050216_Ch14';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load sorted spikes file
[fname pname]=uigetfile('*.*','Get sorted spike traces data');
cd(pname)
tracesData=load(fname);
% allData=load(fname);
% tracesDataWOTime=allData.data.waveform{1,32}'; % for monkey recording
% tracesData=[[0:0.1:3.1]' tracesDataWOTime];

% PCA analysis
tracesDataWOTime=tracesData;
tracesDataWOTime(:,1)=[];
[coeff score]=pca(tracesDataWOTime');
forScatterPlot(:,1)=score(:,PCIndex1);
forScatterPlot(:,2)=score(:,PCIndex2);
% forScatterPlot(:,2)=forScatterPlot(:,2)+forScatterPlot(:,1)*(+0.209); % positive: CCW; negative: CW
% forScatterPlot(:,1)=(forScatterPlot(:,1)*1-9.672)/1; % Center of Mass Realignment
% forScatterPlot(:,2)=(forScatterPlot(:,2)*1.0-11.561)/1; % Center of Mass Realignment

figure
scatter(forScatterPlot(:,1),forScatterPlot(:,2),'k.');
axis([-300 300 -200 200])

% Initialize storage variables
selected_indices=cell(1,numClus);
clusterPlot=cell(1,numClus);
overlayPlot=cell(1,numClus);

% Clustering
switch PCA_method
    case 'sharp'
        for i=1:numClus
            h=impoly();
            nodes = getPosition(h);
            selected_indices{1,i}=inpoly(forScatterPlot,nodes);
            clusterPlot{1,i}=forScatterPlot(selected_indices{1,i},:);
            overlayPlot{1,i}=tracesData';
            overlayPlot{1,i}=overlayPlot{1,i}([true selected_indices{1,i}']',:);
        end
    case 'blur'
        for i=1:numClus
            h=impoly();
            nodes = getPosition(h);
            selected_indices_in=inpoly(forScatterPlot,nodes);
            h=impoly();
            nodes = getPosition(h);
            selected_indices_out=inpoly(forScatterPlot,nodes);
            numIndices=size(selected_indices_in,1);
            selected_indices_final=or(selected_indices_in,(rand(numIndices,1)>testorThreshold_blur).*xor(selected_indices_in,selected_indices_out));
            clusterPlot{1,i}=forScatterPlot(selected_indices_final,:);
            overlayPlot{1,i}=tracesData';
            overlayPlot{1,i}=overlayPlot{1,i}([true selected_indices_final']',:);
        end
    case 'overlap'
        for i=1:numClus
            h=impoly();
            nodes = getPosition(h);
            selected_indices{1,i}=inpoly(forScatterPlot,nodes);
            clusterPlot{1,i}=forScatterPlot;
            overlayPlot{1,i}=tracesData';
        end
        for j=size(forScatterPlot,1):-1:1
            guideArray=[];
            for k=1:numClus
                guideArray=[guideArray selected_indices{1,k}(j)];
            end
            while size(find(guideArray),2)>1
                randIndex=ceil(rand*numClus);
                guideArray(randIndex)=0;
            end
            for ii=1:numClus
                if guideArray(ii)==0
                    clusterPlot{1,ii}(j,:)=[];
                    overlayPlot{1,ii}(j+1,:)=[];
                end
            end
        end
end

% Plot clusters in new colors
figure
for i=1:numClus
    scatter(clusterPlot{1,i}(:,1),clusterPlot{1,i}(:,2),32,strcat(colorOptions(i),'.'));
    hold on
    switch saveOrNot
        case 'yes'
            toSaveFileName=strcat(toSaveNameOptions{i},'_cluster_scatter_',timeStamp,'.dat');
            toSaveMatrix=clusterPlot{1,i};
            command=sprintf('save %s toSaveMatrix -ascii', toSaveFileName);
            eval(command);
    end
end
set(gca,'FontSize',24,'Linewidth',2,'XLim',[-300 300],'YLim',[-200 200],'XTick',[-300:100:300],'YTick',[-200:100:200],'box','off')
xlabel(strcat('Principal Componen',['t' ' ' num2str(PCIndex1)]),'fontsize',36,'FontName','Arial','FontWeight','bold')
ylabel(strcat('Principal Componen',['t' ' ' num2str(PCIndex2)]),'fontsize',36,'FontName','Arial','FontWeight','bold')
%title('Spike Clusters by Principal Component Analysis','fontsize',36,'FontName','Arial','FontWeight','bold')

% Plot spike overlay in new colors
figure
for i=1:numClus
    for j=2:1:size(overlayPlot{1,i},1)
        plot(overlayPlot{1,i}(1,:),overlayPlot{1,i}(j,:),strcat(colorOptions(i),'-'),'Linewidth',2);
        hold on
    end
    switch saveOrNot
        case 'yes'
            toSaveFileName=strcat(toSaveNameOptions{i},'_cluster_spikes_',timeStamp,'.dat');
            toSaveMatrix=overlayPlot{1,i};
            command=sprintf('save %s toSaveMatrix -ascii', toSaveFileName);
            eval(command);
    end
end
xlabel('Time [ms]','fontsize',36,'FontName','Arial','FontWeight','bold')
ylabel('Voltage [uV]','fontsize',36,'FontName','Arial','FontWeight','bold')
set(gca,'FontSize',24,'Linewidth',2,'XLim',[0 3],'YLim',[-300 300],'XTick',[0:0.5:3],'YTick',[-300:100:300],'box','off')