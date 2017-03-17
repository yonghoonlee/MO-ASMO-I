%==========================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.
%==========================================================================
% SETUP_PROBLEM FUNCTION
%==========================================================================
% This function setup the problem and and paths
% Input: problem structure containing 'problem.name'
% Output: problem structure
%==========================================================================

function problem = setup_problem(problem)
    fprintf('%s','Setup problem...');
    
    if (hff_getpoolsize() > 0)
        problem.parallel = true;
    else
        problem.parallel = false;
    end
    
    % Root path and problem path
    problem.rootpath = pwd;
    problem.probpath = fullfile(problem.rootpath,'Problems',problem.name);
    
    % List of toolbox paths
    problem.toolpathlist = string({...
        fullfile('Library','hcparula');
        fullfile('Library','export_fig') ...
    });
    
    % Adding paths
    for i = 1:length(problem.toolpathlist)
        path(path,fullfile(problem.rootpath,problem.toolpathlist{i}));
    end
    
    % Setup problem functions and constraints
    cd(problem.probpath);
    objfun = str2func('obj');           % Objective functions
    nonlconfun = str2func('nonlcon');   % Nonlinear constraint functions
    if (problem.parallel == true)
        poolobj = gcp;
        addAttachedFiles(poolobj,{'obj.m','nonlcon.m'});
    end
    settingsfun = str2func('settings'); % Settings function
    [A,b,Aeq,beq] = setup_lincon();     % Setup linear constraints
    [xlb,xub,flb,fub] = setup_bounds(); % Setup bounds
    problem.nxvar = length(xlb);        % Number of design variables
    problem.nfvar = length(flb);        % Number of obj function variables
    cd(problem.rootpath);
    
    problem.objfun = objfun;
    problem.nonlconfun = nonlconfun;
    problem.settingsfun = settingsfun;
    problem.A = A;
    problem.b = b;
    problem.Aeq = Aeq;
    problem.beq = beq;
    problem.xlb = xlb;
    problem.xub = xub;
    problem.flb = flb;
    problem.fub = fub;

    fprintf('%s\n','done');
end