%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Osyczka-Kundu 1995 Problem
% Osyczka and Kundu, "A new method to solve generalized multicriteria
% optimization problems using the simple genetic algorithm," Structural
% Optimization, 10(2), 1995, pp. 94-99.
%===============================================================================
function result = OK_case001
    casefile = mfilename('fullpath');
    result = runMOASMO(@settings, @obj, @nonlcon, casefile);
end
%===============================================================================
function prob = settings(prob)
    prob.bound.xlb = [0;0;1;0;1;0];
    prob.bound.xub = [10;10;5;6;5;10];
    prob.bound.flb = [-300;0];
    prob.bound.fub = [0;100];
    prob.highfidelity.expensive = false;
    prob.highfidelity.vectorized = true;
    prob.lincon.A = [-1,-1,0,0,0,0;
                        1,1,0,0,0,0;
                        -1,1,0,0,0,0;
                        1,-3,0,0,0,0];
    prob.lincon.b = [-2;6;2;2];
    prob.lincon.Aeq = [];
    prob.lincon.beq = [];
    prob.param = [];
    prob.plotpareto.type = 'pareto2d';
    prob.plotpareto.range = [-300, 0, 0, 100];
    prob.control.maxiter = 20;      % maximum number of iterations
    prob.control.maxerror = 1e-6;   % maximum allowable error
    prob.sampling.initnumber = 6;   % initial
    prob.sampling.valnumber = 6;    % validation
    prob.sampling.upnumber = 4;     % exploitation
    prob.sampling.upexpnumber = 2;  % exploration
    prob.surrogate.method = 'GPR';  % regression model
end
%===============================================================================
function f = obj(x,param)
    f = [-25*(x(:,1) - 2).^2 - (x(:,2) - 2).^2 - (x(:,3) - 1).^2 ...
            - (x(:,4) - 4).^2 - (x(:,5) - 1).^2, ...
        x(:,1).^2 + x(:,2).^2 + x(:,3).^2 + x(:,4).^2 + x(:,5).^2 + ...
            x(:,6).^2];
end
%===============================================================================
function [c,ceq] = nonlcon(x,param)
    c = [(x(:,3) - 3).^2 + x(:,4) - 4, ...
        -(x(:,5) - 3).^2 - x(:,6) + 4];
    ceq = [];
end
%===============================================================================
