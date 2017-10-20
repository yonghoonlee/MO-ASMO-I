%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Plot Figure 04: Convergence
%===============================================================================
try % Open figure window
    figure(fg4);
catch
    fg4 = figure('Color',[1 1 1]);
end
%-------------------------------------------------------------------------------
figure(fg4);
hold off;
if k>1
    for idx = 2:k
        itr = [(idx-1), idx];
        errmax = [max(cell2mat(R.data.c33_valErrorVec(itr(1)))), ...
                  max(cell2mat(R.data.c33_valErrorVec(itr(2))))];
        erravg = [R.data.c34_valErrorAvg(itr(1)), ...
                  R.data.c34_valErrorAvg(itr(2))];
        errmin = [min(cell2mat(R.data.c33_valErrorVec(itr(1)))), ...
                  min(cell2mat(R.data.c33_valErrorVec(itr(2))))];
        ph1 = semilogy(itr, errmax, '-', 'Color', [0.2 0.2 1], 'LineWidth', 1);
        hold on;
        ph2 = semilogy(itr, erravg, '-', 'Color', [0 0 0], 'LineWidth', 2);
        hold on;
        ph3 = semilogy(itr, errmin, '-', 'Color', [1 0.2 0.2], 'LineWidth', 1);
        hold on;
    end
else
    errmax = max(cell2mat(R.data.c33_valErrorVec(k)));
    erravg = R.data.c34_valErrorAvg(k);
    errmin = min(cell2mat(R.data.c33_valErrorVec(k)));
    ph1 = semilogy(k, errmax, '.', 'Color', [0.2 0.2 1], 'LineWidth', 1);
    hold on;
    ph2 = semilogy(k, erravg, '.', 'Color', [0 0 0], 'LineWidth', 2);
    hold on;
    ph3 = semilogy(k, errmin, '.', 'Color', [1 0.2 0.2], 'LineWidth', 1);
    hold on;
end
%-------------------------------------------------------------------------------
ax = gca; ax.FontSize = prob.plotpareto.fontsize;
xlabel('iteration', 'FontSize', prob.plotpareto.fontsize);
ylabel('error: $||\bf{f}_{\rm{P}}-\bf{f}_{\rm{hff}}||$', ...
    'FontSize', prob.plotpareto.fontsize);
legend([ph1, ph2, ph3], ...
    {'largest error', 'average error', 'smallest error'}, ...
    'Location', 'northeast');
%-------------------------------------------------------------------------------
figure(fg4);
if (prob.control.plotexport)
    eval(['export_fig ''', ...
        fullfile( ...
            prob.control.plotpath, [ ...
                prob.control.case, '_fig04_iter', num2str(k,'%04d')] ...
        ), ''' -pdf']);
end
%===============================================================================
