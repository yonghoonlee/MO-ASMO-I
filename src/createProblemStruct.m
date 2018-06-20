%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Creates a default problem structure.
%===============================================================================
function prob = createProblemStruct(settingsfun,objfun,nonlconfun,casename)
    prob = defaultProblemStructure(); % retrieve default values
    [casepath, casefile] = fileparts(which(casename));
    prob.control.path = casepath;
    prob.control.case = casefile; % set case name as calling example name
    if ((prob.control.plot) && (prob.control.plotexport))
        prob.control.plotpath = mfoldername(casename,'plot');
    end
    prob.control.solpath = mfoldername(casename,'solution');
    %---------------------------------------------------------------------------
    % check if settingsfun is provided
    try
        finfo = functions(settingsfun); % if the fn exists, finfo is a struct
    catch
        finfo = []; % otherwise, finfo is an empty variable
    end
    if isstruct(finfo) % if settingsfun is provided, execute it
        prob = feval(settingsfun,prob);
    end
    %---------------------------------------------------------------------------
    % check if objfun is provided
    try
        finfo = functions(objfun); % if the fn exists, finfo is a struct
    catch
        finfo = []; % otherwise, finfo is an empty variable
    end
    if isstruct(finfo) % if objfun is provided, store the function handle
        prob.function.objfun = objfun;
    else
        error('no objective function provided'); % objfun is always required
    end
    %---------------------------------------------------------------------------
    % check if nonlconfun is provided
    try
        finfo = functions(nonlconfun); % if the fn exists, finfo is a struct
    catch
        finfo = []; % otherwise, finfo is an empty variable
    end
    if isstruct(finfo) % if objfun is provided, store the function handle
        prob.function.nonlconfun = nonlconfun;
    end
    %---------------------------------------------------------------------------
    % number of design variables, lower/upper bounds
    if ~(size(prob.bound.xlb) == [0 0])
        prob.nxvar = length(prob.bound.xlb);
    elseif (size(prob.nxvar) == [1 1])
        prob.bound.xlb = -Inf*ones(prob.nxvar,1);
        prob.bound.xub = Inf*ones(prob.nxvar,1);
    else
        error('number of design variables or bounds are required');
    end
    %---------------------------------------------------------------------------
    % number of objective function variables, lower/upper bounds
    if ~(size(prob.bound.flb) == [0 0])
        prob.nfvar = length(prob.bound.flb);
    elseif (size(prob.nfvar) == [1 1])
        prob.bound.flb = -Inf*ones(prob.nfvar,1);
        prob.bound.fub = Inf*ones(prob.nfvar,1);
    else
        error('number of objective function variables or bounds are required');
    end
    %---------------------------------------------------------------------------
    % check if parallel pool is active
    if (hffPoolSize() > 0)
        prob.highfidelity.parallel = true;
        prob.sampling.initconopt.UseParallel = true;
    else
        prob.highfidelity.parallel = false;
    end
    %---------------------------------------------------------------------------
    %prob.gamultiobj.optDO.PopulationSize = min(max(100,5*prob.nxvar),1000); % [100,1000]
end
%===============================================================================
function prob = defaultProblemStructure()
    prob.bound.xlb = [];
    prob.bound.xub = [];
    prob.bound.flb = [];
    prob.bound.fub = [];
    prob.bound.adaptive = true;
    prob.control.case = [];
    prob.control.path = [];
    prob.control.plot = true;           % false, [true]
    prob.control.plotexport = true;     % false, [true]
    prob.control.plotpath = [];
    prob.control.solpath = [];
    prob.control.verbose = 1;           % 0:no, [1:yes], 2:debug
    prob.control.maxiter = 20;          % maximum number of iterations
    prob.control.maxerror = 1e-6;       % maximum allowable error
    prob.control.savedata = true;       % false, [true]
    prob.function.objfun = [];
    prob.function.nonlconfun = [];
    prob.highfidelity.parallel = false;
    prob.highfidelity.expensive = true;
    prob.highfidelity.vectorized = false;
    prob.highfidelity.infeaparam.C = 0.4;
    prob.highfidelity.infeaparam.q = 0.4;
    prob.highfidelity.infeaparam.epsilon = 1e-6;
    prob.lincon.A = [];
    prob.lincon.b = [];
    prob.lincon.Aeq = [];
    prob.lincon.beq = [];
    prob.nxvar = 0;
    prob.nfvar = 0;
    prob.random.seed = 0;
    prob.random.generator = 'twister';
    prob.plotpareto.type = 'pareto2d';   % none, [pareto2d],
                                         % pareto3d, parallel-line
    prob.plotpareto.range = [];
    prob.plotpareto.fontsize = 16;
    prob.sampling.initmethod = 'LHS';
    prob.sampling.initnumber = 5;
    prob.sampling.initconopt = optimoptions('fmincon', ...
        'Algorithm','sqp','Display','iter-detailed','MaxIterations',20, ...
        'MaxFunctionEvaluations',Inf, ...
        'ConstraintTolerance',1e-9,'OptimalityTolerance',1e-4,...
        'StepTolerance',1e-4,'UseParallel',false);
    prob.sampling.initconobjw = 1e-3;    % Weight factor for obj fn
    prob.sampling.initcondispw = 1e-9;   % Weight factor for disperse fn
    prob.sampling.valnumber = 10;
    prob.sampling.upmethod = 'FDL';      % ['FDL'], 'CDD'
    prob.sampling.upnumber = 5;
    prob.sampling.upexpmethod = 'LHS';   % ['LHS']
    prob.sampling.upexpnumber = 5;
    prob.surrogate.method = 'GPR'; % 'RBF', 'RBN', 'SNN', ['GPR'] etc.
                                   % RBF: Radial-basis function network (custom)
                                   % RBN: Radial-basis function network (MATLAB)
                                   % SNN: Shallow ANN (MATLAB)
                                   % GPR: Gaussian Process Regression (MATLAB)
                                   % - kernel function: squared exponential fn.
    prob.surrogate.epsilon = 1;
    prob.surrogate.basisfn = 'TPS';
        % Basis function for RBF: Linear, Cubic, [TPS], Gaussian,
        % MQ (multiquadric), InvMQ (inverse multiquadric), user-defined
    prob.surrogate.snntrainfnc = 'trainlm';
        % ['trainlm']: Levenberg-Marquardt backpropagation
        % 'trainbr': Bayesian regularization backpropagation
    prob.surrogate.snnmaxfail = 30;
    opt = gaoptimset(@gamultiobj);
    opt.PopulationSize = min(max(1000,500*prob.nxvar),10000); % [1k,10k]
    opt.CrossoverFraction = 0.25;
    opt.ParetoFraction = 0.1;
    opt.Generations = 100;
    opt.StallGenLimit = 20;
    opt.TolFun = 5e-5;
    opt.TolCon = 5e-5;
    opt.Vectorized = 'on';
    opt.UseParallel = false;
    if (prob.control.verbose > 0)
        opt.PlotFcns = [@gaplotpareto];
        opt.Display = 'final';
        if (prob.control.verbose > 1)
            opt.Display = 'iter';
        end
    else
        opt.PlotFcns = [];
        opt.Display = 'off';
    end
    optDO = gaoptimset(@gamultiobj);
    optDO.PopulationSize = min(max(1000,5*prob.nxvar),10000); % [1k,10k]
    optDO.CrossoverFraction = 0.25;
    optDO.ParetoFraction = 0.3;
    optDO.Generations = 1000;
    optDO.StallGenLimit = 20;
    optDO.TolFun = 5e-5;
    optDO.TolCon = 5e-5;
    optDO.Vectorized = 'off';
    optDO.UseParallel = true;
    optDO.PlotFcns = [@gaplotpareto];
    optDO.Display = 'iter';
    optDO.OutputFcns = [@runDOFnOut];
    prob.gamultiobj.opt = opt;           % Save options generated
    prob.gamultiobj.optDO = optDO;       % For direct optimization
    prob.gamultiobj.parallel = 0;        % [0:vectorize], 1:parallel
                                         % In most cases, vectorization
                                         % will be faster than parallel
    prob.param = [];
end
%===============================================================================
