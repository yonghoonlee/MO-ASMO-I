% Setting override function for texture design problem with
% rotational tribo-rheometer, Giesekus model, and 3D pseudospectral solver.
% Lee et al., AIAA SciTech 2018.
% Usage: [c,ceq] = NONLCON(x)
% Input: x(n,1:6)
% Output: c(n,1:2),ceq(n,0)
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
end
