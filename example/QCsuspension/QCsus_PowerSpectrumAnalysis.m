%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Quarter car suspension road input analysis
%===============================================================================
clear;
close all;
clc;

plotPreparation;

[mpath,mname] = fileparts(mfilename('fullpath'));
load(fullfile(mpath,'IRI_737b.mat'), 'road*');
vlist = [5,10,15,30,45,60,80];
samprate = 4;

v = vlist(1);
t = [road_x/v];
z = [road_z];
if (length(vlist) > 1)
    for k = 2:length(vlist)
        v = vlist(k);
        t = [t, 2*t(end) - t(end-1) + road_x/v];
        z = [z, road_z];
    end
end

tu = [0:(10^(-samprate)):round(t(end),samprate)]';
zu = [interp1(t,z,tu,'pchip')]';
fs = 10^(samprate);
segLen = 10^(samprate);
nfft = length(zu);
[P,F] = pwelch(zu,nfft,0,nfft,fs,'power');

% Plot Figure 01: Power spectrum density of road input
try % Open figure window
    figure(fg1);
catch
    fg1 = figure('Color',[1 1 1]);
end
semilogx(F,10*log10(P));
ax = gca;
ax.FontSize = 16;
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

% export plot
eval(['export_fig ''', fullfile(mpath,'plot',[mname,'_road_psd']), ''' -pdf']);
