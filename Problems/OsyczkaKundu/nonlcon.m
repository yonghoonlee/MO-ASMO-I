%==========================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.
%==========================================================================
% NONLCON FUNCTION
% Osyczka and Kundu, "A new method to solve generalized multicriteria
% optimization problems using the simple genetic algorithm," Structural
% Optimization, 10(2), 1995, pp. 94-99.
%==========================================================================
% Nonlinear constraint function of Osyczka and Kundu 1995 problem
% Input: x(n,1:6)
% Output: c(n,1:2), ceq(n,0)
%==========================================================================

function [c,ceq] = nonlcon(x)
    c = [(x(:,3) - 3).^2 + x(:,4) - 4, ...
        -(x(:,5) - 3).^2 - x(:,6) + 4];
    ceq = [];
end
