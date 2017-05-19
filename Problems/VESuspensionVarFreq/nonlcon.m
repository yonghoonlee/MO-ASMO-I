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

    if (size(x,1) == 1)
        if (size(x,2) > 1)
            xin = x;
        else
        	error('Number of design variable does not match');
        end
    else
        if (size(x,2) == 1)
            xin = x';
        else
            xin = x;
        end
    end

    nmode = floor(size(x,2)/2);
    xub = repmat(p.xub',size(xin,1),1);

    K = [10.^xin(:,1:nmode),10.^xin(:,end)];
    K = sum(K,2);
    c = [c; K - (10.^xub(:,1) + 10.^xub(:,end))];

end
