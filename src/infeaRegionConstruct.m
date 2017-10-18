%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Infeasible region domain model construction using SVDD
% Reference: Malak Jr., R.J. and Paredis, C.J.J., "Using Support Vector Machines
% to Formalize the Valid Input Domain of Predictive Models in Systems Design
% Problems," Journal of Mechanical Design, 132, Oct 2010, pp.101001,
% doi: 10.1115/1.4002151
%===============================================================================
function infeaRegion = infeaRegionConstruct(infeaX, prob)
    p.C = prob.highfidelity.infeaparam.C;
    p.q = prob.highfidelity.infeaparam.q;
    p.epsilon = prob.highfidelity.infeaparam.epsilon;
    p.x = infeaX;
    nx = size(p.x,1);
    if (nx < 5)
        infeaRegion = [];
        return;
    end
    %---------------------------------------------------------------------------
    d0 = ones(nx,1);
    lb = zeros(size(d0));
    ub = p.C*ones(size(d0));
    Aeq = ones(1,nx);
    beq = 1;
    npool = hffPoolSize();
    if npool > 0
        par = true;
    else
        par = false;
    end
    opt = optimoptions('fmincon','Algorithm','sqp',...
        'ConstraintTolerance',1e-9,'Display','iter-detailed',...
        'FiniteDifferenceType','central','MaxFunctionEvaluations',Inf,...
        'MaxIterations',Inf,'OptimalityTolerance',1e-9,...
        'StepTolerance',1e-12,'UseParallel',par);
    [dopt,~] = fmincon(@(d)obj(d,p),d0,[],[],Aeq,beq,lb,ub,[],opt);
    %---------------------------------------------------------------------------
    infeaRegion.dopt = dopt;
    infeaRegion.infeaX = infeaX;
    infeaRegion.C = p.C;
    infeaRegion.q = p.q;
    infeaRegion.epsilon = p.epsilon;
    infeaRegion.x = p.x;
    %===========================================================================
    function f = obj(d,p)
        x = p.x;
        npx = size(x,1);
        f1 = 0;
        for m = 1:npx
            f1 = f1 + d(m)*exp(0);
        end
        f2 = 0;
        for n = 1:npx
            for m = 1:npx
                % Use Gaussian kernel as Mercer kernel of nonlinear transform
                f2 = f2 + d(m)*d(n)*exp(-p.q*sum((x(m,:)-x(n,:)).^2));
            end
        end
        f = -(f1-f2);
    end
    %===========================================================================
end
%===============================================================================
