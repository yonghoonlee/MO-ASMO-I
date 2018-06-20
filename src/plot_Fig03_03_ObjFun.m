%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Plot Figure 03: Pareto frontier in the objective space and validation point
%===============================================================================
try figure(fg3); % Open figure window
catch, fg3 = figure('Color',[1 1 1]); end; fg3.Position = [870 780 560 220];
hold off;
%-------------------------------------------------------------------------------
% Prepare data from current predicted Pareto set
xplot1 = c18_parSurXFea;
fplot1 = c19_parSurFFea;
[nx1,mx1] = size(xplot1);
[nf1,mf1] = size(fplot1);
xfplot1 = sortrows([xplot1, fplot1], mx1+1);
xplot1 = xfplot1(:,1:mx1);
fplot1 = xfplot1(:,(mx1+1):end);
%-------------------------------------------------------------------------------
% Prepare data from validation sample data set
xplot2 = c27_valSurXFea;
fplot2 = c28_valSurFFea;
hplot2 = c29_valHffFFea;
[nx2,mx2] = size(xplot2);
[nf2,mf2] = size(fplot2);
[nh2,mh2] = size(hplot2);
xfhplot2 = sortrows([xplot2, fplot2, hplot2], mx2+mf2+1);
xplot2 = xfhplot2(:,1:mx2);
fplot2 = xfhplot2(:,(mx2+1):(mx2+mf2));
hplot2 = xfhplot2(:,(mx2+mf2+1):end);
%-------------------------------------------------------------------------------
% Plot objective function space
if (nx2>0)
    cm = plasma(nx2);
    switch lower(prob.plotpareto.type)
        case 'pareto2d'
            ph1 = plot(fplot1(:,1), fplot1(:,2), 'k-', 'MarkerSize', 10, ...
                'LineWidth', 2); hold on;
            for idx = 1:nx2
                plot(fplot2(idx,1), fplot2(idx,2), 'o', 'MarkerSize', 11, ...
                    'LineWidth', 2.5, 'Color', (cm(idx,:) + [0 0 0])/2);
                hold on;
                plot(hplot2(idx,1), hplot2(idx,2), 's', 'MarkerSize', 10, ...
                    'Color', cm(idx,:), 'MarkerFaceColor', cm(idx,:));
                hold on;
                plot([fplot2(idx,1) hplot2(idx,1)], ...
                    [fplot2(idx,2) hplot2(idx,2)], ':', ...
                    'LineWidth', 1.5, 'Color', (cm(idx,:) + [0 0 0])/2);
                hold on;
            end
            ph2 = plot(fplot2(1,1), fplot2(1,2), 'o', 'MarkerSize', 11, ...
                'LineWidth', 2.5, 'Color', cm(1,:)); hold on;
            ph3 = plot(hplot2(1,1), hplot2(1,2), 's', 'MarkerSize', 10, ...
                'Color', cm(1,:), 'MarkerFaceColor', cm(1,:)); hold on;
            if (numel(prob.plotpareto.range) ~= 0)
                axis(prob.plotpareto.range);
            end
            ax = gca; ax.FontSize = prob.plotpareto.fontsize;
            xlabel('$f_1$: objective-1', 'FontSize', prob.plotpareto.fontsize);
            ylabel('$f_2$: objective-2', 'FontSize', prob.plotpareto.fontsize);
            title(['[iteration ', num2str(k), ']'], ...
                'FontSize', (prob.plotpareto.fontsize - 2));
            legend([ph1, ph2, ph3], {'predicted Pareto frontier', ...
                'validation sample point', 'high fidelity model result'}, ...
                'Location', 'best');
            %-------------------------------------------------------------------
        case 'pareto3d'
            %-------------------------------------------------------------------
        case 'parallel-line'
            %-------------------------------------------------------------------
        otherwise
    end
    %---------------------------------------------------------------------------
    if (prob.control.plotexport)
        figure(fg3);
        eval(['export_fig ''', ...
            fullfile( ...
                prob.control.plotpath, [ ...
                    prob.control.case, '_fig03_iter', num2str(k,'%04d')] ...
            ), ''' -pdf']);
    end
end
%===============================================================================
