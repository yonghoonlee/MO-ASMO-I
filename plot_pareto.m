% The PLOT_PARETO script is called by main script to plot Pareto set.
% Usage: PLOT_PARETO

% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.

% Pareto set plot in design space
figure(fg0); clf;
cmap = hcparula(size(DATA{k,5},1));
subplot(2,1,1);
for i = 1:size(DATA{k,5},1)
    tmpdat = DATA{k,5}(i,:);
    tmpdat = (tmpdat - problem.xlb')./(problem.xub' - problem.xlb');
    plot(1:(problem.nxvar),tmpdat,'-','Color',cmap(i,:));
    hold on;
end
clear tmpdat;
subplot(2,1,2);
for i = 1:size(DATA{k,6},1)
    tmpdat = DATA{k,6}(i,:);
    tmpdat = (tmpdat - problem.flb')./(problem.fub' - problem.flb');
    plot(1:(problem.nfvar),tmpdat,'-','Color',cmap(i,:));
    hold on;
end
clear tmpdat;

% Pareto set plot in objective function space
figure(fg1); clf;
p1 = plot(DATA{k,6}(:,1),DATA{k,6}(:,2),'x',...
    'MarkerEdgeColor',[0 0 0]); hold on;
if (isfield(problem,'plotrange'))
    if (isfield(problem.plotrange,'xmin') ...
            && isfield(problem.plotrange,'xmax') ...
            && isfield(problem.plotrange,'ymin') ...
            && isfield(problem.plotrange,'ymax'))
        axis([min(problem.plotrange.xmin, min(DATA{k,6}(:,1))) ...
            max(problem.plotrange.xmax, max(DATA{k,6}(:,1))) ...
            min(problem.plotrange.ymin, min(DATA{k,6}(:,2))) ...
            max(problem.plotrange.ymax, max(DATA{k,6}(:,2)))]);
    end
end
set(gca,'FontSize',16);
xlabel('$f_1$'); ylabel('$f_2$');
l1 = legend([p1],{'predicted Pareto set'},...
    'Location','southwest','Interpreter','latex','FontSize',16);
legend('boxoff');
if (problem.control.plotexport ~= 0)
    eval(['export_fig ',fullfile(problem.probpath,...
        ['fig_pareto_',num2str(k),'.pdf']), ' -pdf']);
end

% Pareto set plot for inexpensive functions in objective function space
if (problem.highfidelity.expensive == 0)        % If fn eval is not costly
    figure(fg2); clf;
    p3 = plot(DATA{k,7}(:,1),DATA{k,7}(:,2),'o',...
        'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',6,'LineWidth',1.5);
    hold on;
    p2 = plot(DATA{k,6}(:,1),DATA{k,6}(:,2),'x',...
        'MarkerEdgeColor',[0 0 0]);
    for i = 1:size(DATA{k,6},1)
        p4 = plot([DATA{k,6}(i,1);DATA{k,7}(i,1)],...
            [DATA{k,6}(i,2);DATA{k,7}(i,2)],':','Color',[0.5 0.5 0.5]);
    end
    if (isfield(problem,'plotrange'))
        if (isfield(problem.plotrange,'xmin') ...
                && isfield(problem.plotrange,'xmax') ...
                && isfield(problem.plotrange,'ymin') ...
                && isfield(problem.plotrange,'ymax'))
            axis([min(problem.plotrange.xmin, ...
                    min([DATA{k,6}(:,1); DATA{k,7}(:,1)])) ...
                max(problem.plotrange.xmax, ...
                    max([DATA{k,6}(:,1); DATA{k,7}(:,1)])) ...
                min(problem.plotrange.ymin, ...
                    min([DATA{k,6}(:,2); DATA{k,7}(:,2)])) ...
                max(problem.plotrange.ymax, ...
                    max([DATA{k,6}(:,2); DATA{k,7}(:,2)]))]);
        end
    end
    set(gca,'FontSize',16);
    xlabel('$f_1$'); ylabel('$f_2$');
    l1 = legend([p2,p3],{'predicted Pareto set','high-fidelity function'},...
        'Location','southwest','Interpreter','latex','FontSize',16);
    legend('boxoff');
    if (problem.control.plotexport ~= 0)
        eval(['export_fig ',fullfile(problem.probpath,...
            ['fig_pareto_compare_',num2str(k),'.pdf']), ' -pdf']);
    end
end