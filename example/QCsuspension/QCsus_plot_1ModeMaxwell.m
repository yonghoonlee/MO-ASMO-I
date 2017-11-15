close all
clear;
clc;

% 1 MODE Case 351
[mpath,mname] = fileparts(mfilename('fullpath'));
load(fullfile(mpath,'solution','QCsus_1ModeMaxwell351_final.mat'),'result');
x_D1_351 = cell2mat(result.data.c07_PoolXFea(end));
f_D1_351 = cell2mat(result.data.c08_PoolHffFFea(end));
[x_D1sort351,f_D1sort351,n_D1didx351] = ndSort(x_D1_351,f_D1_351);

% 1 MODE Case 451
[mpath,mname] = fileparts(mfilename('fullpath'));
load(fullfile(mpath,'solution','QCsus_1ModeMaxwell451_final.mat'),'result');
x_D1_451 = cell2mat(result.data.c07_PoolXFea(end));
f_D1_451 = cell2mat(result.data.c08_PoolHffFFea(end));
[x_D1sort451,f_D1sort451,n_D1didx451] = ndSort(x_D1_451,f_D1_451);

fg1 = figure('Color',[1 1 1]);
plot(f_D1_351(:,1),f_D1_351(:,2),'k.'); hold on;
plot(f_D1_451(:,1),f_D1_451(:,2),'r.'); hold on;
