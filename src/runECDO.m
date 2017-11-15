%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Main routine for running Epsilon-constraint-based Direct Optimization.
%===============================================================================
function R = runECDO(settingsfun, objfun, nonlconfun, casefile, varargin)
    prob = createProblemStruct(settingsfun,objfun,nonlconfun,casefile);
    R = createResultStruct(prob); % create result struct
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
    R.initPop = initPop;
    R.initScr = initScr;
    %---------------------------------------------------------------------------
    % option setting
    optEC = optimoptions('fmincon','Algorithm','sqp');
    optEC.Display = 'iter-detailed';
    optEC.MaxFunctionEvaluations = Inf;
    optEC.MaxIterations = 40;
    optEC.ScaleProblem = true;
    optEC.OptimalityTolerance = 1e-4;
    optEC.TypicalX = sum(initPop,1)/size(initPop,1);
    A = prob.lincon.A;
    b = prob.lincon.b;
    Aeq = prob.lincon.Aeq;
    beq = prob.lincon.beq;
    lb = prob.bound.xlb;
    ub = prob.bound.xub;
    objfun = prob.function.objfun;
    nonlconfun = prob.function.nonlconfun;
    param = prob.param;
    %---------------------------------------------------------------------------
    % parallel settings
    if (hffPoolSize == 0)
        optEC.UseParallel = false;
    else
        optEC.UseParallel = true;
    end
    %---------------------------------------------------------------------------
    % multiobjective optimization 1: finding anchor point for each
    % objective function
    X = [initPop, initScr];
    for idx = 1:prob.nfvar
        R.data.c01_anchorX(idx) = {[]};
        R.data.c02_anchorF(idx) = {[]};
        R.time.c01_02(idx) = {[]};
    end
    for idx = 1:prob.nfvar
        tic;
        Xtmp = sortrows(X,(prob.nxvar+idx));
        x0 = Xtmp(1,1:(prob.nxvar));
        [xo{idx},~] = fmincon(@(x)objecdo_anchor(objfun,x,param,idx),...
            x0,A,b,Aeq,beq,lb',ub',@(x)nonlconfun(x,param),optEC);
        fo{idx} = objfun(xo{idx},param);
        to{idx} = toc;
    end
    utopiapt = [];
    for idx = 1:prob.nfvar
        R.data.c01_anchorX(idx) = xo(idx);
        R.data.c02_anchorF(idx) = fo(idx);
        R.time.c01_02(idx) = to(idx);
        utopiapt = [utopiapt, cell2mat(R.data.c02_anchorF(idx))];
    end
    utopiapt1 = min(utopiapt,[],2);
    utopiapt2 = max(utopiapt,[],2);
    R.utopiapt = utopiapt1;
    R.antiutopiapt = utopiapt2;
    
    save(fullfile(prob.control.solpath, ...
        [prob.control.case, '_01_anchorpoints.mat']), 'R');
    
    %---------------------------------------------------------------------------
    % multiobjective optimization 2: epsilon-constraint methods
    clear xo;
    clear fo;
    clear to;
    ofun = @hobjfun;
    cfun = @hconfun;
    if (prob.nfvar == 2)
        nPdiv = 24;
        fPcon = utopiapt1(2) + transpose(1:(nPdiv-1)).*((utopiapt2(2)-utopiapt1(2))/nPdiv);
        optEC.UseParallel = false;
        [sortPop, sortScr, sortInd] = ndSort(initPop, initScr);
        
        funcs.objfun = objfun;
        funcs.nonlconfun = nonlconfun;
        opt.A = A; opt.b = b; opt.Aeq = Aeq; opt.beq = beq;
        opt.lb = lb'; opt.ub = ub'; opt.optEC = optEC;
        
        parfor idx = 1:(nPdiv-1)
            tic;
            % Get initial point for epsilon-constraint optimization
            fPconstr = fPcon(idx);
            dfPcon = abs(fPconstr - sortScr(:,2)).*exp(sortInd.^2);
            [~,idx2] = min(dfPcon);
            x0 = sortPop(idx2,:);
%             [xo{idx},~] = fmincon(@(x)ofun(objfun,x,param), ...
%                 x0, A, b, Aeq, beq, lb', ub', ...
%                 @(x)cfun(nonlconfun,x,param,fPconstr), optEC);
            [xo{idx},~,~,~] = runECDOobjconstr(x0,opt,param,fPconstr,funcs);
            fo{idx} = objfun(xo{idx},param);
            to{idx} = toc;
        end
        xoo = [reshape(R.data.c01_anchorX{1},1,prob.nxvar);
            reshape(R.data.c01_anchorX{2},1,prob.nxvar)];
        foo = [reshape(R.data.c02_anchorF{1},1,prob.nfvar);
            reshape(R.data.c02_anchorF{2},1,prob.nfvar)];
        for idx = 1:(nPdiv-1)
            xoo = [xoo; reshape(xo{idx},1,prob.nxvar)];
            foo = [foo; reshape(fo{idx},1,prob.nfvar)];
        end
        [xoo,foo,ioo] = ndSort(xoo,foo);
        xPo = xoo(find(ioo==1),:);
        fPo = foo(find(ioo==1),:);
        R.data.c03_ParetoX(1) = {xPo};
        R.data.c04_ParetoF(1) = {fPo};
        R.x = xPo;
        R.fval = fPo;
        save(fullfile(prob.control.solpath, ...
            [prob.control.case, '_02_ParetoSet.mat']), 'R');
    else
        error('3D or more dimensions in obj space is not supported yet');
    end
    %---------------------------------------------------------------------------
    % Plot
    plotPreparation;
    plotECDOFig01ParetoFrontier;
end
%===============================================================================
function f = objecdo_anchor(objfun,x,param,idx)
    fout = feval(objfun,x,param);
    f = fout(idx);
end
%===============================================================================
