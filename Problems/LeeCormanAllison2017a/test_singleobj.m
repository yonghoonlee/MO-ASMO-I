clear; clc;

objfun = str2func('obj');           % Objective functions
nonlconfun = str2func('nonlcon');   % Nonlinear constraint functions
problem.probpath = pwd;

% Parallel pool
poolobj = gcp('nocreate');  % If no pool, do not create new one.
if isempty(poolobj)
    npool = 0;
    problem.parallel = false;
else
    npool = poolobj.NumWorkers;
    problem.parallel = true;
    poolobj = gcp;
    addAttachedFiles(poolobj,{'obj.m','nonlcon.m'});
end

[A,b,Aeq,beq] = setup_lincon();     % Setup linear constraints
[xlb,xub,flb,fub] = setup_bounds(); % Setup bounds
problem.nxvar = length(xlb);        % Number of design variables
problem.nfvar = length(flb);        % Number of obj function variables
problem.objfun = objfun;
problem.nonlconfun = nonlconfun;
problem.A = A;
problem.b = b;
problem.Aeq = Aeq;
problem.beq = beq;
problem.xlb = xlb;
problem.xub = xub;
problem.flb = flb;
problem.fub = fub;
problem = settings(problem);        % Problem-specific params

%% SINGLE-OBJECTIVE OPTIMIZATIONS

gopt = gaoptimset('ga');
gopt.PopulationSize = 100;
gopt.CrossoverFraction = 0.35;
gopt.Generations = 100;
gopt.TolFun = 1e-9;
gopt.TolCon = 1e-6;
gopt.Display = 'iter';
gopt.MutationFcn = @mutationadaptfeasible;
gopt.UseParallel = true;

p = problem.p;
p.single = 1;
[xopt1,fopt1] = ga(@(x)obj(x,p),5,[],[],[],[],xlb,xub,[],gopt);

p.single = 2;
[xopt2,fopt2] = ga(@(x)obj(x,p),5,[],[],[],[],xlb,xub,[],gopt);

save('test_singleobj.mat');

fprintf('Optimization for Objective-1\n');
fprintf('xopt: %20.16f \n',xopt1);
fprintf('fopt: %20.16f \n',fopt1);

fprintf('Optimization for Objective-2\n');
fprintf('xopt: %20.16f \n',xopt2);
fprintf('fopt: %20.16f \n',fopt2);
