%==========================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.
%==========================================================================
% SETUP_LINCON FUNCTION
% Osyczka and Kundu, "A new method to solve generalized multicriteria
% optimization problems using the simple genetic algorithm," Structural
% Optimization, 10(2), 1995, pp. 94-99.
%==========================================================================
% Linear constraint setup function of Osyczka and Kundu 1995 problem
% Input: none
% Output: A, b, Aeq, beq
%==========================================================================

function [A,b,Aeq,beq] = setup_lincon()
    A = [-1,-1,0,0,0,0;
        1,1,0,0,0,0;
        -1,1,0,0,0,0;
        1,-3,0,0,0,0];
    b = [-2;6;2;2];
    Aeq = [];
    beq = [];
end
