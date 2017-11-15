close all
clear;
clc;

% 1 MODE
[mpath,~] = fileparts(mfilename('fullpath'));
load(fullfile(mpath,'solution','QCsus_1ModeMaxwell451_final.mat'),'result');
x_D1 = cell2mat(result.data.c07_PoolXFea(end));
f_D1 = cell2mat(result.data.c08_PoolHffFFea(end));
[x_D1sort,f_D1sort,n_D1didx] = ndSort(x_D1,f_D1);

% 2 MODE
[mpath,~] = fileparts(mfilename('fullpath'));
load(fullfile(mpath,'solution','QCsus_2ModeMaxwell451_final.mat'),'result');
x_D2 = cell2mat(result.data.c07_PoolXFea(end));
f_D2 = cell2mat(result.data.c08_PoolHffFFea(end));
[x_D2sort,f_D2sort,n_D2didx] = ndSort(x_D2,f_D2);

% 3 MODE
[mpath,~] = fileparts(mfilename('fullpath'));
load(fullfile(mpath,'solution','QCsus_3ModeMaxwell451_final.mat'),'result');
x_D3 = cell2mat(result.data.c07_PoolXFea(end));
f_D3 = cell2mat(result.data.c08_PoolHffFFea(end));
[x_D3sort,f_D3sort,n_D3didx] = ndSort(x_D3,f_D3);

% 4 MODE
[mpath,~] = fileparts(mfilename('fullpath'));
load(fullfile(mpath,'solution','QCsus_4ModeMaxwell451_final.mat'),'result');
x_D4 = cell2mat(result.data.c07_PoolXFea(end));
f_D4 = cell2mat(result.data.c08_PoolHffFFea(end));
[x_D4sort,f_D4sort,n_D4didx] = ndSort(x_D4,f_D4);

% 10 MODE
[mpath,~] = fileparts(mfilename('fullpath'));
load(fullfile(mpath,'solution','QCsus_10ModeMaxwell451_final.mat'),'result');
x_D10 = cell2mat(result.data.c07_PoolXFea(end));
f_D10 = cell2mat(result.data.c08_PoolHffFFea(end));
[x_D10sort,f_D10sort,n_D10didx] = ndSort(x_D10,f_D10);

% LN RESULTS
[mpath,~] = fileparts(mfilename('fullpath'));
load(fullfile(mpath,'solution','QCsus_LogNormal451_final.mat'),'result');
x_LN = cell2mat(result.data.c07_PoolXFea(end));
f_LN = cell2mat(result.data.c08_PoolHffFFea(end));
[x_LNsort,f_LNsort,n_LNdidx] = ndSort(x_LN,f_LN);

% Plot 1: Designs in linear-scaled objective function space
fg1 = figure('Color',[1 1 1]);
fg1.Position = [60 120 560 420];
ph1 = plot(f_D1(:,1),f_D1(:,2),'k.'); hold on;
ph2 = plot(f_D2(:,1),f_D2(:,2),'r.'); hold on;
ph3 = plot(f_D3(:,1),f_D3(:,2),'b.'); hold on;
ph4 = plot(f_D4(:,1),f_D4(:,2),'g.'); hold on;
ph10 = plot(f_D10(:,1),f_D10(:,2),'c.'); hold on;
phLN = plot(f_LN(:,1),f_LN(:,2),'m.'); hold on;
legend([ph1,ph2,ph3,ph4,ph10,phLN],{'1-mode Maxwell','2-mode Maxwell',...
    '3-mode Maxwell','4-mode Maxwell','10-mode Maxwell','Log-Normal'},...
    'Location','northeast');

% Plot 2: Designs in log-scaled objective function space
fg2 = figure('Color',[1 1 1]);
fg2.Position = [620 120 560 420];
ph1 = loglog(f_D1(:,1),f_D1(:,2),'k.'); hold on;
ph2 = loglog(f_D2(:,1),f_D2(:,2),'r.'); hold on;
ph3 = loglog(f_D3(:,1),f_D3(:,2),'b.'); hold on;
ph4 = loglog(f_D4(:,1),f_D4(:,2),'g.'); hold on;
ph10 = loglog(f_D10(:,1),f_D10(:,2),'c.'); hold on;
phLN = loglog(f_LN(:,1),f_LN(:,2),'m.'); hold on;
legend([ph1,ph2,ph3,ph4,ph10,phLN],{'1-mode Maxwell','2-mode Maxwell',...
    '3-mode Maxwell','4-mode Maxwell','10-mode Maxwell','Log-Normal'},...
    'Location','southwest');

% Plot 3: Pareto set comparison
fg3 = figure('Color',[1 1 1]);
fg3.Position = [1180 120 560 420];
ph1 = plot(f_D1sort(n_D1didx==1,1),f_D1sort(n_D1didx==1,2),'k.'); hold on;
ph2 = plot(f_D2sort(n_D2didx==1,1),f_D2sort(n_D2didx==1,2),'r.'); hold on;
ph3 = plot(f_D3sort(n_D3didx==1,1),f_D3sort(n_D3didx==1,2),'b.'); hold on;
ph4 = plot(f_D4sort(n_D4didx==1,1),f_D4sort(n_D4didx==1,2),'g.'); hold on;
ph10 = plot(f_D10sort(n_D10didx==1,1),f_D10sort(n_D10didx==1,2),'c.'); hold on;
phLN = plot(f_LNsort(n_LNdidx==1,1),f_LNsort(n_LNdidx==1,2),'m.'); hold on;
legend([ph1,ph2,ph3,ph4,ph10,phLN],{'1-mode Maxwell','2-mode Maxwell',...
    '3-mode Maxwell','4-mode Maxwell','10-mode Maxwell','Log-Normal'},...
    'Location','southwest');

