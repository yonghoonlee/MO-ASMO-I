%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Plot Figure 04: Pareto frontier evolution
%===============================================================================
try figure(fg4); % Open figure window
catch, fg4 = figure('Color',[1 1 1]); end; fg4.Position = [1330 770 560 220];
%-------------------------------------------------------------------------------
cm = flipud(plasma(k));
th = linspace(1,3,k);
figure(fg4);
hold off;
for idx = 1:k
    % Prepare data
    fplot = cell2mat(R.data.c19_parSurFFea(idx));
    [nf,mf] = size(fplot);
    fplot = sortrows(fplot, 1);
    %---------------------------------------------------------------------------
    ph{idx} = plot(fplot(:,1), fplot(:,2), '-', 'LineWidth', th(idx), ...
        'Color', cm(idx,:)); hold on;
end
if (numel(prob.plotpareto.range) ~= 0)
    axis(prob.plotpareto.range);
end
ax = gca; ax.FontSize = prob.plotpareto.fontsize;
xlabel('$f_1$: objective-1', 'FontSize', prob.plotpareto.fontsize);
ylabel('$f_2$: objective-2', 'FontSize', prob.plotpareto.fontsize);
title(['[iteration ', num2str(k), ']'], ...
	'FontSize', (prob.plotpareto.fontsize - 2));
colormap(cm);
caxis([0.5 k+0.5]);
cb = colorbar;
cb.Ticks = 1:k;
cb.Label.String = 'iteration number';
cb.Label.Interpreter = 'latex';
%-------------------------------------------------------------------------------
if (prob.control.plotexport)
    figure(fg4);
    eval(['export_fig ''', ...
        fullfile( ...
            prob.control.plotpath, [ ...
                prob.control.case, '_fig04_iter', num2str(k,'%04d')] ...
        ), ''' -pdf']);
end
%===============================================================================
