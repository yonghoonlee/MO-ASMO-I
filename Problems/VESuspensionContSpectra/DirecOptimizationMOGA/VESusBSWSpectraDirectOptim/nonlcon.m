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
% Usage: [c,ceq] = NONLCON(x)
% Input: x
% Output: c,ceq
%   x: Points in design space
%   c: Inequality constraint function values of given design points
%   ceq: Equality constraint function values of given design points

% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.

function [c,ceq] = nonlcon(x,p)
    c = [];
    ceq = [];

    % Decompose design variables
    n_e = x(1);             % 0.01 - 10
    n_g = x(2);             % 0.01 - 10
    G_N0 = 10^x(3);         % 0 - 6
    lambda_c = 10^x(4);     % -4 - 0
    lambda_max = 10^x(5);   % -3 - 1
    Ge = x(6);
    k1 = 10^x(7);           % Suspension stiffness [N/m]
        
    H = @(tau) n_e.*G_N0.*((tau./lambda_c).^(-n_g) + (tau./lambda_max).^n_e);

    Kmax = Ge + real(integral(@(s) H(exp(s)),-36,log(lambda_max)));
    if abs(Kmax) > 1e16
    	c = [NaN];
    end
    c = Kmax - 10^(p.xub(end));

end
