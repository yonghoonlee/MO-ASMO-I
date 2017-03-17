%==========================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.
%==========================================================================
% OBJ FUNCTION
% Osyczka and Kundu, "A new method to solve generalized multicriteria
% optimization problems using the simple genetic algorithm," Structural
% Optimization, 10(2), 1995, pp. 94-99.
%==========================================================================
% Objective function of Osyczka and Kundu 1995 problem
% Input: x(n,1:6)
% Output: f(n,1:2)
%==========================================================================

function f = obj(x)
    f = [-25*(x(:,1) - 2).^2 - (x(:,2) - 2).^2 - (x(:,3) - 1).^2 ...
            - (x(:,4) - 4).^2 - (x(:,5) - 1).^2, ...
        x(:,1).^2 + x(:,2).^2 + x(:,3).^2 + x(:,4).^2 + x(:,5).^2 + ...
            x(:,6).^2];
end
