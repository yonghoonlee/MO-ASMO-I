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

% Full Factorial Design
N = 4;
x = (fullfact(N*ones(1,5)) - 1)./(N-1);
Xlb = repmat(xlb',size(x,1),1);
Xub = repmat(xub',size(x,1),1);
X = Xlb + (Xub - Xlb).*x;
F = zeros(size(X,1),2);

% Full Factorial Analysis
p = problem.p;
if (problem.parallel == true)
    parfor i = 1:size(X,1)
        F(i,:) = transpose(obj(X(i,:),p));
    end
else
    for i = 1:size(X,1)
        F(i,:) = transpose(obj(X(i,:),p));
    end
end