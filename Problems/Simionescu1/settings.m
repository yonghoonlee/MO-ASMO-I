% Setting override function of Osyczka and Kundu 1995 problem
% Osyczka and Kundu, "A new method to solve generalized multicriteria
% optimization problems using the simple genetic algorithm," Structural
% Optimization, 10(2), 1995, pp. 94-99.
% Usage: problem = SETTINGS(problem)
% Input: problem
% Output: problem
%   problem: Problem definition structure

% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.

function problem = settings(problem)
    problem.highfidelity.expensive = 0;     % Not expensive
    problem.highfidelity.vectorized = 1;    % Function evaluation in vector
    problem.plotrange.xmin = -0.08;
    problem.plotrange.xmax = 0;
    problem.plotrange.ymin = -0.12;
    problem.plotrange.ymax = -0.08;
end