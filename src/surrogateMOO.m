%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Runs multiobjective optimization on the surrogate model using gamultiobj
%===============================================================================
function [xopt, fopt, out] = surrogateMOO(surrogate, infeaRegion, initpop, prob)
    if (prob.control.verbose > 0)
        fprintf('%s','Solving multiobjective optimization...');
    end
    optGa = prob.gamultiobj.opt;
    A = prob.lincon.A;
    b = prob.lincon.b;
    Aeq = prob.lincon.Aeq;
    beq = prob.lincon.beq;
    lb = prob.bound.xlb;
    ub = prob.bound.xub;
    nxvar = prob.nxvar;
    nonlconfun = prob.function.nonlconfun;
    param = prob.param;
    %---------------------------------------------------------------------------
    % Modify options
    % Put initial population
    if size(initpop,1) > 1
        optGa.InitialPopulation = initpop;
    end
    % If parallel GA is enabled (not recommended)
    if ((prob.gamultiobj.parallel ~= 0) && (hffPoolSize ~= 0))
        opt.Vectorized = 'off';
        opt.UseParallel = true;
    end
    %---------------------------------------------------------------------------
    % Run multiobjective optimization using gamultiobj
    [xga, fga, exitflag, output, population, score] = gamultiobj( ...
        @(x) surrogateEval(x, surrogate), nxvar, A, b, Aeq, beq, lb, ub, ...
        @(x) nonlconEval(x, param, nonlconfun, infeaRegion), optGa);
    %---------------------------------------------------------------------------
    % Make solution unique
    xsize = size(xga,2);
    xcomb = [xga, fga];
    xcomb = unique(xcomb, 'rows');
    xopt = xcomb(:,1:xsize);
    fopt = xcomb(:,(xsize + 1):end);
    %---------------------------------------------------------------------------
    % Return output from gamultiobj
    out.exitflag = exitflag;
    out.output = output;
    out.population = population;
    out.score = score;
    %---------------------------------------------------------------------------
    if (prob.control.verbose > 0)
        fprintf('%s\n','done');
    end
    %===========================================================================
    function [c, ceq] = nonlconEval(x, param, nonlconfun, infeaRegion)
        [c1, ceq] = feval(nonlconfun, x, param);
        c2 = infeaRegionCheck(infeaRegion, x);
        c = [c1, c2];
    end
    %===========================================================================
end
%===============================================================================