%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Open figure windows
%===============================================================================
set(0,'DefaultAxesTickLabelInterpreter','latex');
set(0,'DefaultColorbarTickLabelInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultPolaraxesTickLabelInterpreter','latex');
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultTextarrowshapeInterpreter','latex');
set(0,'DefaultTextboxshapeInterpreter','latex');
set(0,'DefaultFigurePosition',[20 150 560 420]);
%-------------------------------------------------------------------------------
close all;
%-------------------------------------------------------------------------------
fg1 = figure('Color',[1 1 1]);
fg1.Position = [10 550 560 420];
%-------------------------------------------------------------------------------
fg2 = figure('Color',[1 1 1]);
fg2.Position = [450 500 560 420];
%-------------------------------------------------------------------------------
fg3 = figure('Color',[1 1 1]);
fg3.Position = [890 450 560 420];
%-------------------------------------------------------------------------------
fg4 = figure('Color',[1 1 1]);
fg4.Position = [1330 400 560 420];
%===============================================================================
