%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Plot Figure 01-02: Variables
%===============================================================================
try figure(fg1); % Open figure window
catch, fg1 = figure('Color',[1 1 1]); end; fg1.Position = [10 800 560 220];
hold off;
try figure(fg2); % Open figure window
catch, fg1 = figure('Color',[1 1 1]); end; fg2.Position = [440 790 560 220];
hold off;
%-------------------------------------------------------------------------------
% Prepare and normalize data from all stored high fidelity results
xplot = [c07_PoolXFea; c27_valSurXFea];
fplot = [c08_PoolHffFFea; c29_valHffFFea];
[nx,mx] = size(xplot);
[nf,mf] = size(fplot);
xfplot = sortrows([xplot, fplot], mx+1);
xplot = xfplot(:,1:mx);
xplot = (xplot - repmat(reshape(prob.bound.xlb,1,mx),nx,1)) ...
    ./(repmat(reshape(prob.bound.xub,1,mx),nx,1) ...
        - repmat(reshape(prob.bound.xlb,1,mx),nx,1));
fplot = xfplot(:,(mx+1):end);
fplot = (fplot - repmat(reshape(prob.bound.flb,1,mf),nf,1)) ...
    ./(repmat(reshape(prob.bound.fub,1,mf),nf,1) ...
        - repmat(reshape(prob.bound.flb,1,mf),nf,1));
cm = 0.6 + 0.4*plasma(nx);
%-------------------------------------------------------------------------------
% Plot design varaibles for all stored high fidelity results
figure(fg1);
for idx = 1:nx
    plot(1:mx, xplot(idx,:), '-', 'Color', cm(idx,:)); hold on;
end
axis([1, mx, 0, 1]);
ax = gca; ax.XTick = linspace(1,mx,mx); ax.YTick = [];
xlabel('design variable', 'FontSize', prob.plotpareto.fontsize);
ax.FontSize = prob.plotpareto.fontsize;
title(['[iteration ', num2str(k), ']'], ...
    'FontSize', (prob.plotpareto.fontsize - 2));
%-------------------------------------------------------------------------------
% Plot objective function variables for all stored high fidelity results
figure(fg2);
for idx = 1:nf
    plot(1:mf, fplot(idx,:), '-', 'Color', cm(idx,:)); hold on;
end
axis([1, mf, 0, 1]);
ax = gca; ax.XTick = linspace(1,mf,mf); ax.YTick = [];
xlabel('objective function variable', 'FontSize', prob.plotpareto.fontsize);
ax.FontSize = prob.plotpareto.fontsize;
%-------------------------------------------------------------------------------
% Prepare and normalize data from current predicted Pareto set
xplot = c18_parSurXFea;
fplot = c19_parSurFFea;
[nx,mx] = size(xplot);
[nf,mf] = size(fplot);
xfplot = sortrows([xplot, fplot], mx+1);
xplot = xfplot(:,1:mx);
xplot = (xplot - repmat(reshape(prob.bound.xlb,1,mx),nx,1)) ...
    ./(repmat(reshape(prob.bound.xub,1,mx),nx,1) ...
        - repmat(reshape(prob.bound.xlb,1,mx),nx,1));
fplot = xfplot(:,(mx+1):end);
fplot = (fplot - repmat(reshape(prob.bound.flb,1,mf),nf,1)) ...
    ./(repmat(reshape(prob.bound.fub,1,mf),nf,1) ...
        - repmat(reshape(prob.bound.flb,1,mf),nf,1));
cm = plasma(nx);
%-------------------------------------------------------------------------------
% Plot design variables for current predicted Pareto set
figure(fg1);
for idx = 1:nx
    plot(1:mx, xplot(idx,:), '-', 'Color', cm(idx,:), ...
        'LineWidth', 2); hold on;
end
axis([1, mx, 0, 1]); ax = gca; ax.XTick = linspace(1,mx,mx);
%-------------------------------------------------------------------------------
% Plot objective function variables for current predicted Pareto set
figure(fg2);
for idx = 1:nf
    plot(1:mf, fplot(idx,:), '-', 'Color', cm(idx,:), ...
        'LineWidth', 2); hold on;
end
axis([1, mf, 0, 1]); ax = gca; ax.XTick = linspace(1,mf,mf);
%-------------------------------------------------------------------------------
if (prob.control.plotexport)
    figure(fg1);
    eval(['export_fig ''', ...
        fullfile( ...
            prob.control.plotpath, [ ...
                prob.control.case, '_fig01_iter', num2str(k,'%04d')] ...
        ), ''' -pdf']);
    figure(fg2);
    eval(['export_fig ''', ...
        fullfile( ...
            prob.control.plotpath, [ ...
                prob.control.case, '_fig02_iter', num2str(k,'%04d')] ...
        ), ''' -pdf']);
end
%===============================================================================
