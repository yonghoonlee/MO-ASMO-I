%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Plot Figure 08: Pareto set of High-fidelity function results
%===============================================================================
try figure(fg8); % Open figure window
catch, fg8 = figure('Color',[1 1 1]); end; fg8.Position = [1320 520 560 220];
hold off;
%-------------------------------------------------------------------------------
figure(fg8);
ph1xd = cell2mat(R.data.c07_PoolXFea(k,1));
ph1fd = cell2mat(R.data.c08_PoolHffFFea(k,1));
if size(ph1xd,1)>0
    [ph1xdSort,ph1fdSort,ph1idSort] = ndSort(ph1xd,ph1fd);
    maxi = max(ph1idSort);
    cm = gray(ceil(1.6*maxi));
    for idx = 1:maxi
        ph1fdR = ph1fdSort(ph1idSort==idx,:);
        plot(ph1fdR(:,1),ph1fdR(:,2),'.','Color',cm(idx,:),'MarkerSize',10);
        hold on;
    end
    clear ph1xd ph1fd ph1xdSort ph1fdSort ph1idSort maxi idx ph1fdR;
else
    clear ph1xd ph1fd;
end
ph1fd = cell2mat(R.data.c04_smpHffFFea(k,1));
if size(ph1fd,1)>0
    plot(ph1fd(:,1),ph1fd(:,2),'ro','MarkerSize',6,'LineWidth',2);
    hold on;
end
clear ph1fd;
%-------------------------------------------------------------------------------
ax = gca; ax.FontSize = prob.plotpareto.fontsize;
xlabel('$f_1$: objective-1', 'FontSize', prob.plotpareto.fontsize);
ylabel('$f_2$: objective-2', 'FontSize', prob.plotpareto.fontsize);
title(['[iteration ', num2str(k), ']'], ...
    'FontSize', (prob.plotpareto.fontsize - 2));
%-------------------------------------------------------------------------------
if (prob.control.plotexport)
    figure(fg8);
    eval(['export_fig ''', ...
        fullfile( ...
            prob.control.plotpath, [ ...
                prob.control.case, '_fig08_iter', num2str(k,'%04d')] ...
        ), ''' -pdf']);
end
%===============================================================================
