function [x,f,eflag,outpt] = runECDOobjconstr(x0,opt,param,fPconstr,funcs)

    xLast = []; % Last place computeall was called
    myf = []; % Use for objective at xLast
    myc = []; % Use for nonlinear inequality constraint
    myceq = []; % Use for nonlinear equality constraint

    fun = @objfn; % the objective function, nested below
    cfun = @constr; % the constraint function, nested below

    % Call fmincon
    [x,f,eflag,outpt] = fmincon(@(x)fun(x,param,fPconstr,funcs), ...
        x0, opt.A, opt.b, opt.Aeq, opt.beq, ...
        opt.lb, opt.ub, @(x)cfun(x,param,fPconstr,funcs), opt.optEC);

    function y = objfn(x,param,fPconstr,funcs)
        if ~isequal(x,xLast) % Check if computation is necessary
            [myf,myc,myceq] = computeall(x,param,fPconstr,funcs);
            xLast = x;
        end
        % Now compute objective function
        y = myf(1);
    end
    function [c,ceq] = constr(x,param,fPconstr,funcs)
        if ~isequal(x,xLast) % Check if computation is necessary
            [myf,myc,myceq] = computeall(x,param,fPconstr,funcs);
            xLast = x;
        end
        % Now compute constraint functions
        c = myc; % In this case, the computation is trivial
        ceq = myceq;
    end
end

function [F,C,Ceq] = computeall(x,param,fPconstr,funcs)
    objfun = funcs.objfun;
    nonlconfun = funcs.nonlconfun;
    
    % objective function and constraint evaluation result
    F = feval(objfun,x,param);
    [C,Ceq] = feval(nonlconfun,x,param);
    
    % smooth scaling of epsilon-constraint function
    ecdoConst = exp(F(2)-fPconstr);
    if (ecdoConst > exp(3))
        ecdoConst = exp(3) + exp(3)*2*sqrt(3)*sqrt(F(2)-fPconstr) - 6*exp(3);
    end
    ecdoConst = ecdoConst - 1;
    
    % combining inequality constraints with epsilon constraint
    C = [reshape(C,numel(C),1); ecdoConst];
    
end