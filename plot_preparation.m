% The PLOT_PREPARATION script is called by main script to prepare figures.
% Usage: PLOT_PREPARATION

% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.

% Prepare plotting
if (problem.control.plot ~= 0)
    fg0 = figure('Color',[1 1 1]);              % fg0 with white background
    set(fg0,'Position',[50 100 560 420]);       % set position of fg0
    fg1 = figure('Color',[1 1 1]);              % fg1 with white background
    set(fg1,'Position',[630 100 560 420]);       % set position of fg1
    if (problem.highfidelity.expensive == 0)    % If fn eval is not costly
        fg2 = figure('Color',[1 1 1]);          % fg2 with white background
        set(fg2,'Position',[1210 100 560 420]);  % set position of fg2
    end
    if isfield(problem,'plotcustom')
        for i = 1:length(problem.plotcustom)
            fgcustom{i} = figure('Color',[1 1 1]);
            set(fgcustom{i},'Position',[(1000+40*i) (520-40*i) 560 420]);
        end
    end
    set(0, 'defaultTextInterpreter', 'latex');  % Text interpreter: LaTeX
    set(0, 'defaultLegendInterpreter', 'latex');% Legend interpreter: LaTeX
end