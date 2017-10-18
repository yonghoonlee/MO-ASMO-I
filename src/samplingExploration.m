%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Generates samples to explore unexplored region of the design space.
%===============================================================================
function xf = samplingExploration(prob,k,varargin)
    if (prob.control.verbose > 0)
        fprintf('Generating samples to explore unexplored design space...');
    end
    %---------------------------------------------------------------------------
    if (k == 0)
        method = prob.sampling.initmethod;
        number = prob.sampling.initnumber;
        prevPoints = [];
    else
        method = prob.sampling.upexpmethod;
        number = prob.sampling.upexpnumber;
        if (nargin > 2)
            prevPoints = varargin{1};
        else
            prevPoints = [];
        end
    end
    %---------------------------------------------------------------------------
    dimension = prob.nxvar;
    A = prob.lincon.A;
    b = prob.lincon.b;
    Aeq = prob.lincon.Aeq;
    beq = prob.lincon.beq;
    lb = prob.bound.xlb;
    ub = prob.bound.xub;
    nonlconfun = prob.function.nonlconfun;
    p = prob.param;
    opt = prob.sampling.initconopt;
    wei1 = prob.sampling.initconobjw;
    wei2 = prob.sampling.initcondispw;
    %---------------------------------------------------------------------------
    switch lower(method)
        case 'lhs'
            ex = -1; % To enter the while loop
            while (ex < 0)
                ex = 2;
                % Sampling and descaling
                xt = samplingLHS(number,dimension);
                xt = samplingDeScale(xt,lb,ub);
                % Sample adjustment (to satisfy constraints)
                fprintf('%s','adjust samples to satisfy constraints...');
                xf = xt;
                if (prob.control.verbose == 2); fprintf('\n'); end
                for i = 1:size(xt,1)
                    xtmp = xt(i,:);
                    [xs,f,e,o] = fmincon( ...
                        @(x)samplingCObj(x,xtmp,wei1,wei2,prevPoints),...
                        xtmp,A,b,Aeq,beq,lb,ub,@(x)nonlconfun(x,p),opt);
                    xf(i,:) = reshape(xs,1,numel(xs));
                    if (size(prevPoints,1) ~= 0)
                        prevPoints = [prevPoints; reshape(xs,1,numel(xs))];
                    end
                    ex = min(ex,e);
                    if (prob.control.verbose == 2)
                        disp(strcat('exitflag:',num2str(e),...
                            '/objfn:',num2str(f),...
                            '/constrviolation:',num2str(o.constrviolation)));
                    end
                end
            end
        otherwise
            error(strcat(method,'::not supported.'));
    end
    %---------------------------------------------------------------------------
    if (prob.control.verbose == 2); fprintf('...'); end
    if (prob.control.verbose > 0)
        fprintf('%s\n','done');
    end
    %===========================================================================
    function f = samplingCObj(x,x0,wei1,wei2,prevPoints)
        xm = repmat(x,size(prevPoints,1),1);
        f1 = wei1*sum((x - x0).^2);
        if (size(prevPoints,1) == 0)
            f2 = 0;
        else
            f2 = wei2*sum(-log(max(sum((xm - prevPoints).^2,2),1e-12))) ...
                /size(prevPoints,1);
        end
        f = f1 + f2;
    end
    %===========================================================================
end
%===============================================================================
