%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Plot Figure 01 for ECDO: Pareto frontier in the objective space
%===============================================================================
try % Open figure window
    figure(fg1);
catch
    fg1 = figure('Color',[1 1 1]);
end
cm = plasma(5);
%-------------------------------------------------------------------------------
% Design exploration plot
ph1 = plot(initScr(:,1),initScr(:,2),'.','Color',cm(4,:),'MarkerSize',9); hold on;
%-------------------------------------------------------------------------------
% Initial population plot
ph2 = plot(R.fval(:,1),R.fval(:,2),'o','Color',cm(2,:),...
        'MarkerFaceColor',cm(2,:),'MarkerSize',5); hold on;
%-------------------------------------------------------------------------------
% Plot labels
ax = gca; ax.FontSize = prob.plotpareto.fontsize;
xlabel('$f_1$: objective-1', 'FontSize', prob.plotpareto.fontsize);
ylabel('$f_2$: objective-2', 'FontSize', prob.plotpareto.fontsize);
legend([ph1, ph2], {'design explored by MO-ASMO', ...
    'direct optimization result'}, ...
    'Location', 'best');
%-------------------------------------------------------------------------------
figure(fg1);
if (prob.control.plotexport)
    eval(['export_fig ''', ...
        fullfile( ...
            prob.control.plotpath, [ ...
                prob.control.case, '_fig01'] ...
        ), ''' -pdf']);
end
%===============================================================================
