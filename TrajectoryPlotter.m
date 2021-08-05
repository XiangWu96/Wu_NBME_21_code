% Plots the trajectory of the mouse . The trajectory should have been
% generated by MouseTracker.m by analyzing a digital camera video.

% Tip: remember to change the pixelsize, which can be derived by plotting
% one single picture from the video and matches the its dimensions with
% real word objects

clear all
clc

%% Defination of parameters
pixelsize=0.067;        % Unit: cm

%% Read the file

[fname, pname]=uigetfile('*.*','please select the tracjectory file to open');
cd(pname)
a=load(fname);
rescaled_a=a.*pixelsize;

figure;
lineplot = plot(rescaled_a(:,1),rescaled_a(:,2),'Linewidth', 1);
lineplot.Color=[0,0,0,0.1];
hold on
scatterplot=scatter(rescaled_a(:,1),rescaled_a(:,2),1, 'filled','k');
scatterplot.MarkerFaceAlpha = 0.1;
axis ij
axis image
xlim([20,110])
ylim([0,70])

xlabel('x (cm)','fontsize',15,'FontName','Arial','FontWeight','bold')
ylabel('y (cm)','fontsize',15,'FontName','Arial','FontWeight','bold')
set(gca,'FontSize',15,'Linewidth',2)
title(fname)