%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Main routine for running Direct Optimization using NSGA-II algorithm.
%===============================================================================
function R = runDO(settingsfun, objfun, nonlconfun, casefile, varargin)
    prob = createProblemStruct(settingsfun,objfun,nonlconfun,casefile);
    rng(prob.random.seed,prob.random.generator); % set random seed
    %---------------------------------------------------------------------------
    % load from MO-ASMO result
    initPop = [];
    if (nargin == 5)
        try
            load(fullfile(prob.control.solpath, ...
                [cell2mat(varargin(1))]), 'result');
            try
                initPop = cell2mat(result.data.c07_PoolXFea(end));
            catch
                try
                    initPop = cell2mat(result.c07_PoolXFea);
                catch
                    initPop = [];
                end
            end
            try
                initScr = cell2mat(result.data.c08_PoolHffFFea(end));
            catch
                try
                    initScr = cell2mat(result.c08_PoolHffFFea);
                catch
                    initScr = [];
                end
            end
        catch
            initPop = [];
            initScr = [];
        end
    end
    %---------------------------------------------------------------------------
    % option setting
    optGa = prob.gamultiobj.optDO;
    A = prob.lincon.A;
    b = prob.lincon.b;
    Aeq = prob.lincon.Aeq;
    beq = prob.lincon.beq;
    lb = prob.bound.xlb;
    ub = prob.bound.xub;
    nxvar = prob.nxvar;
    objfun = prob.function.objfun;
    nonlconfun = prob.function.nonlconfun;
    param = prob.param;
    %---------------------------------------------------------------------------
    % Initial population and parallel settings
    if size(initPop,1) > 1
        if optGa.PopulationSize < size(initPop,1)
            [initPop, initScr, ~] = ndSort(initPop, initScr);
            initPop = initPop(1:(optGa.PopulationSize), :);
            initScr = initScr(1:(optGa.PopulationSize), :);
        end
        optGa.InitialPopulation = initPop;
    end
    if (hffPoolSize == 0)
        optGa.UseParallel = false;
    end
    % multiobjective optimization
    tic;
    [xopt, fopt, exitflag, output, population, score] = gamultiobj( ...
        @(x) objfun(x, param), nxvar, A, b, Aeq, beq, lb, ub, ...
        @(x) nonlconfun(x, param), optGa);
    timeDO = toc;
    %---------------------------------------------------------------------------
    % Save final result
    clear result;
    result = struct('xopt',xopt,'fopt',fopt,'exitflag',exitflag,...
        'output',output,'population',population,'score',score,...
        'timeDO',timeDO);
    try
        result.pophistory = evalin('base', 'pophistory');
    catch
        result.pophistory = [];
    end
    try
        result.scrhistory = evalin('base', 'scrhistory');
    catch
        result.scrhistory = [];
    end
    R = result;
    save(fullfile(prob.control.solpath, ...
        [prob.control.case, '_DOfinal.mat']), 'result');
    %---------------------------------------------------------------------------
    % Plot
    %plotPreparation;
    plotDOFig01ParetoFrontier;
end
%===============================================================================
