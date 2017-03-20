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
    fg1 = figure('Color',[1 1 1]);              % fg1 with white background
    if (problem.highfidelity.expensive == 0)    % If fn eval is not costly
        fg2 = figure('Color',[1 1 1]);          % fg2 with white background
    end
    set(0, 'defaultTextInterpreter', 'latex');  % Text interpreter: LaTeX
    set(0, 'defaultLegendInterpreter', 'latex');% Legend interpreter: LaTeX
end