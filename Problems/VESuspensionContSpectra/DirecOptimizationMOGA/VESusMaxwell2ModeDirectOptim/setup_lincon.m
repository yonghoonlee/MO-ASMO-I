% Objective function of Lee-Corman-Allison viscoelastic suspension problem.
% * Lee, Corman, Ewoldt, and Allison, "A Multiobjective Adaptive Surrogate
% Modeling-Based Optimization (MO-ASMO) Framework Using Efficient Sampling
% Strategies," ASME 2017 IDETC/CIE, Cleveland, OH, 2017, to appear.
% * Corman, Rao, Bharadwaj, Allison, and Ewoldt, "Setting Material Function
% Design Targets for Linear Viscoelastic Materials and Structures," Journal
% of Mechanical Design, 138(5), p.051402(12pp), doi: 10.1115/1.4032698.
% * Allison, Guo, Han, "Co-Design of an Active Suspension Using
% Simultaneous Dynamic Optimization," Journal of Mechanical Design, 136(8),
% p.081003(14pp), doi: 10.1115/1.4027335.
%
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
    nmode = 2;
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    
    if (nmode > 1)
        ncol = 2*nmode + 1;
        for i = 1:(nmode - 1)
            Aadd = zeros(1,ncol);
            Aadd(nmode+i) = 1;
            Aadd(nmode+i+1) = -1;
            badd = -0.5;
            A = [A; Aadd];
            b = [b; badd];
        end
    end
end
