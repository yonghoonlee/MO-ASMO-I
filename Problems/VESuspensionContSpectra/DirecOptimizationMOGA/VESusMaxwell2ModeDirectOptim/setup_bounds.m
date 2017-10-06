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
    nmode = 2;
    xlb = [0.0*ones(nmode,1); -3.5*ones(nmode,1); 2.5];
    xub = [3.3*ones(nmode,1); 1*ones(nmode,1); 2.8];
    flb = [0;1.1];
    fub = [4;1.5];
end
