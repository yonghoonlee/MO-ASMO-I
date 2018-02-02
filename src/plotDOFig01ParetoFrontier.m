%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Plot Figure 01 for DO: Pareto frontier in the objective space
%===============================================================================
try % Open figure window
    figure(fg1);
catch
    fg1 = figure('Color',[1 1 1]);
end
cm = plasma(5);
%-------------------------------------------------------------------------------
% Design exploration plot
[~,~,nGen] = size(result.scrhistory);
for idx1 = 1:nGen
    scr = result.scrhistory(:,:,idx1);
    ph2 = plot(scr(:,1),scr(:,2),'.','Color',cm(5,:),'MarkerSize',9); hold on;
end
%-------------------------------------------------------------------------------
% Initial population plot
try
    size(initScr,1);
catch
    initScr = [];
end
if (size(initScr,1) ~= 0)
    ph1 = plot(initScr(:,1),initScr(:,2),'o','Color',cm(4,:),...
        'MarkerFaceColor',cm(4,:),'MarkerSize',5); hold on;
end
%-------------------------------------------------------------------------------
% DO result plot
ph3 = plot(fopt(:,1),fopt(:,2),'s','Color',cm(1,:),'MarkerSize',6,...
    'MarkerFaceColor',cm(1,:)); hold on;
%-------------------------------------------------------------------------------
% Plot labels
ax = gca; ax.FontSize = prob.plotpareto.fontsize;
xlabel('$f_1$: objective-1', 'FontSize', prob.plotpareto.fontsize);
ylabel('$f_2$: objective-2', 'FontSize', prob.plotpareto.fontsize);
try
    legend([ph1, ph2, ph3], {'initial population from MO-ASMO', ...
        'explored design in direct optimization', 'direct optimization result'}, ...
        'Location', 'best');
catch
    legend([ph2, ph3], {'explored design in direct optimization', ...
        'direct optimization result'}, ...
        'Location', 'best');
end
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
