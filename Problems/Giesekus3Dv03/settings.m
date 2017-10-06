% Setting override function for texture design problem with
% rotational tribo-rheometer, Giesekus model, and 3D pseudospectral solver.
% Lee et al., AIAA SciTech 2018.
% Usage: problem = SETTINGS(problem)
% Input: problem
% Output: problem

% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.

function problem = settings(problem)
    problem.highfidelity.expensive = 1;     % Expensive
    problem.highfidelity.vectorized = 0;    % Fn.eval. cannot be vectorized
    %problem.plotrange.xmin = 0;
    %problem.plotrange.xmax = 20;
    %problem.plotrange.ymin = -15;
    %problem.plotrange.ymax = 5;
    problem.p.Omega = 10;
    problem.p.Ntex = 10;
    problem.p.ri = 0.5e-3;
    problem.p.ro = 20e-3;
    problem.p.Nr = 10;
    problem.p.Nth = 10;
    problem.p.Nz = 5;
    problem.p.eta = 9.624e-3;
    problem.p.rho = 873.4;
    problem.sampling.initnumber = 60;
    problem.sampling.upnumber = 15;
    problem.sampling.upexpnumber = 15;
end