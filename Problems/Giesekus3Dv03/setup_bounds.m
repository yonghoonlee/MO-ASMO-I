% Setting override function for texture design problem with
% rotational tribo-rheometer, Giesekus model, and 3D pseudospectral solver.
% Lee et al., AIAA SciTech 2018.
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
    Nr = 10;
    Nth = 10;
    hmin =  269e-6;
    hmax = 1000e-6;
    
    xlb_geom = hmin*ones(((Nr+1)*Nth),1);
    xub_geom = hmax*ones(((Nr+1)*Nth),1);
    xtyp_fld = [0.013250667;
                0.006625333;
                0.001654582;
                1.42e-4;
                0.05;
                0.05];
    xlb_fld = 0.1*xtyp_fld;
    xub_fld = 10.0*xtyp_fld;
    xlb = [xlb_geom; xlb_fld];
    xub = [xub_geom; xub_fld];
    flb = [0;-15];
    fub = [20;5];
end
