% This function initializes global parameters to the default preset values.
% Problem-specific parameters can be overrided by "settings.m" function in
% the problem directory.
% Usage: problem = SETUP_DEFAULT_PARAMETER(problem)
% Input: problem
% Output: problem
%   problem: problem definition structure

% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.

function problem = setup_default_parameter(problem)
    fprintf('%s','Setup parameter...');
    
    problem.control.verbose = 2;            % 0:no, 1:yes, [2:debug]
    problem.control.maxiter = 20;           % maximum number of iterations
    problem.control.plot = 1;               % 0:no, [1:yes]
    problem.control.plotexport = 1;         % 0:no, [1:yes]
    problem.control.savedata = 1;           % 0:no, [1:yes]
    
    problem.highfidelity.expensive = 1;     % 0:no, [1:yes] (expensive)
    problem.highfidelity.vectorized = 0;    % [0:no], 1:yes
    
    problem.sampling.initmethod = 'LHS';
    problem.sampling.initnumber = 5;
    problem.sampling.initconopt = optimoptions('fmincon', ...
        'Algorithm','sqp','Display','none','MaxIterations',1000, ...
        'MaxFunctionEvaluations',1000*problem.nxvar, ...
        'TolCon',1e-6,'TolFun',1e-8,'TolX',1e-8);
%        'TolCon',1e-12,'TolFun',1e-16,'TolX',1e-16);
    problem.sampling.initconobjw = 1e-3;    % Weight factor for obj fn
    problem.sampling.initcondispw = 1e-9;   % Weight factor for disperse fn
    
    problem.sampling.valnumber = 10;
    
    problem.sampling.upmethod = 'FDL';      % ['FDL'], 'CDD'
    problem.sampling.upnumber = 5;
    problem.sampling.upexpmethod = 'LHS';   % ['LHS']
    problem.sampling.upexpnumber = 5;
    
    problem.surrogate.method = 'RBF';
    problem.surrogate.epsilon = 1;
    problem.surrogate.basisfn = 'TPS';
        % Basis function for RBF: Linear, Cubic, [TPS], Gaussian,
        % MQ (multiquadric), InvMQ (inverse multiquadric), user-defined
    
    opt = gaoptimset(@gamultiobj);
    opt.PopulationSize = min(max(1000,500*problem.nxvar),10000); % [1k,10k]
    opt.CrossoverFraction = 0.25;
    opt.ParetoFraction = 0.1;
    opt.Generations = 200;
    opt.StallGenLimit = 20;
    opt.TolFun = 5e-5;
    opt.TolCon = 5e-5;
    opt.Vectorized = 'on';
    opt.UseParallel = false;
    if (problem.control.verbose > 0)
        opt.PlotFcns = [@gaplotpareto];
        opt.Display = 'final';
        if (problem.control.verbose > 1)
            opt.Display = 'iter';
        end
    else
        opt.PlotFcns = [];
        opt.Display = 'off';
    end
    problem.gamultiobj.opt = opt;           % Save options generated
    problem.gamultiobj.parallel = 0;        % [0:vectorize], 1:parallel
                                            % In most cases, vectorization
                                            % will be faster than parallel
    problem.p = [];
    
    fprintf('%s\n','done');
end