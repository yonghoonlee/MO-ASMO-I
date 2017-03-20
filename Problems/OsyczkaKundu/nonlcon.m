% Nonlinear constraint functions of Osyczka and Kundu 1995 problem
% Osyczka and Kundu, "A new method to solve generalized multicriteria
% optimization problems using the simple genetic algorithm," Structural
% Optimization, 10(2), 1995, pp. 94-99.
% Usage: [c,ceq] = NONLCON(x)
% Input: x(n,1:6)
% Output: c(n,1:2),ceq(n,0)
%   x: Points in design space
%   c: Inequality constraint function values of given design points
%   ceq: Equality constraint function values of given design points

% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.

function [c,ceq] = nonlcon(x)
    c = [(x(:,3) - 3).^2 + x(:,4) - 4, ...
        -(x(:,5) - 3).^2 - x(:,6) + 4];
    ceq = [];
end
