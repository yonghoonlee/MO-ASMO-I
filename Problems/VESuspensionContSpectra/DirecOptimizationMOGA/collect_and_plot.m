close all; clear; clc;
set(0,'DefaultTextInterpreter','latex'); % change the text interpreter
set(0,'DefaultLegendInterpreter','latex'); % change the legend interpreter
set(0,'DefaultAxesTickLabelInterpreter','latex'); % change the tick interpreter
FntSz = 16;

% Retrieve results of 2-mode Maxwell model
load('VESusMaxwell2ModeDirectOptim/RunDO_Maxwell2Mode.mat');
Maxwell2.pop = [];
Maxwell2.scr = [];
for i = 1:size(pophistory,3)
    Maxwell2.pop = [Maxwell2.pop; pophistory(:,:,i)];
    Maxwell2.scr = [Maxwell2.scr; scrhistory(:,:,i)];
end

% Retrieve results of Log-normal model
load('VESusContSpectraDirectOptim/RunDO_LogNormal.mat');
LogNormal.pop = [];
LogNormal.scr = [];
for i = 1:size(pophistory,3)
    LogNormal.pop = [LogNormal.pop; pophistory(:,:,i)];
    LogNormal.scr = [LogNormal.scr; scrhistory(:,:,i)];
end

% Retrieve results of BSW model
load('VESusBSWSpectraDirectOptim/RunDO_BSW.mat');
BSW.pop = [];
BSW.scr = [];
for i = 1:size(pophistory,3)
    BSW.pop = [BSW.pop; pophistory(:,:,i)];
    BSW.scr = [BSW.scr; scrhistory(:,:,i)];
end

hf1 = figure('Color',[1 1 1]);
hp1 = plot(Maxwell2.scr(:,1),Maxwell2.scr(:,2),'r.','MarkerSize',10); hold on;
hp2 = plot(LogNormal.scr(:,1),LogNormal.scr(:,2),'b.','MarkerSize',10); hold on;
hp3 = plot(BSW.scr(:,1),BSW.scr(:,2),'m.','MarkerSize',10); hold on;
hax = gca;
hax.FontSize = FntSz;
xlabel('$f_1=max(abs(\ddot{x}_1(t)))$ (comfort)','FontSize',FntSz);
ylabel('$f_2=max(x_2(t)-z(t))-min(x_2(t)-z(t))$ (handling)','FontSize',FntSz);
title('Objective function space','FontSize',FntSz);
legend([hp1,hp2,hp3],{'2-mode Maxwell model','Log-normal model','BSW model'},'FontSize',FntSz);
eval(['export_fig ','objective_space', ' -pdf']);

