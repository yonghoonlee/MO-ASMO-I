% Linear constraint setup function of Osyczka and Kundu 1995 problem
% Osyczka and Kundu, "A new method to solve generalized multicriteria
% optimization problems using the simple genetic algorithm," Structural
% Optimization, 10(2), 1995, pp. 94-99.
% Usage: [A,b,Aeq,beq] = SETUP_LINCON()
% Input:
% Output: A,b,Aeq,beq
% $$A \mathbf{x} \le \mathbf{b}$$
% $$A_{eq} \mathbf{x} = \mathbf{b}_{eq}$$
%   A: Linear inequality constraints matrix
%   b: Linear inequality constraints vector
%   Aeq: Linear equality constraints matrix
%   beq: Linear equality constraints vector

% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.

function [A,b,Aeq,beq] = setup_lincon()
    A = [];
    b = [];
    Aeq = [];
    beq = [];
end