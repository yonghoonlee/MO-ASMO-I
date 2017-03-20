% Bounds setup function of Osyczka and Kundu 1995 problem
% Osyczka and Kundu, "A new method to solve generalized multicriteria
% optimization problems using the simple genetic algorithm," Structural
% Optimization, 10(2), 1995, pp. 94-99.
% Usage: [xlb,xub,flb,fub] = SETUP_BOUNDS()
% Input:
% Output: xlb,xub,flb,fub
%   xlb: Lower bound for design variable
%   xub: Upper bound for design variable
%   flb: Predicted lower bound for objective function variable
%   fub: Predicted upper bound for objective function variable

% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.

function [xlb,xub,flb,fub] = setup_bounds()
    xlb = [0;0;1;0;1;0];
    xub = [10;10;5;6;5;10];
    flb = [-500;0];
    fub = [0;100];
end
