%==========================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.
%==========================================================================
% SETUP_BOUNDS FUNCTION
% Osyczka and Kundu, "A new method to solve generalized multicriteria
% optimization problems using the simple genetic algorithm," Structural
% Optimization, 10(2), 1995, pp. 94-99.
%==========================================================================
% Linear constraint setup function of Osyczka and Kundu 1995 problem
% Input: none
% Output: A, b, Aeq, beq
%==========================================================================

function [xlb,xub,flb,fub] = setup_bounds()
    xlb = [0;0;1;0;1;0];
    xub = [10;10;5;6;5;10];
    flb = [-500;0];
    fub = [0;100];
end
