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

    problem = settings();

    Nr = problem.p.Nr;
    Nth = problem.p.Nth;

    x_geom_lb = 269e-6*ones((Nr+1)*Nth,1);
    x_geom_ub = 1000e-6*ones((Nr+1)*Nth,1) + x_geom_lb;

    x_fluid_lb = [0.010; 0.003; 0.0010; 0.00010; 0.01; 0.01];
    x_fluid_ub = [0.020; 0.010; 0.0020; 0.00020; 0.10; 0.10];

    xlb = [x_geom_lb; x_fluid_lb];
    xub = [x_geom_ub; x_fluid_ub];
    flb = [0.5e-4;-0.1];
    fub = [2.0e-4;0.01];
end
