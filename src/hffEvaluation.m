%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% High fidelity function evaluation
%===============================================================================
function f = hffEvaluation(x, prob)
    if (prob.control.verbose > 0)
        fprintf('%s','High fidelity function evaluation...');
    end
    nx = size(x,1);                 % Number of designs
    mf = prob.nfvar;                % Number of obj function variables
    f = zeros(nx,mf);               % Allocate result matrix
    npool = hffPoolSize();          % Number of parallel pools
    if (prob.control.verbose > 0)
        fprintf('%s',strcat(num2str(size(x,1)),'-points...'));
    end
    %---------------------------------------------------------------------------
    if (npool == 0)
        if (prob.highfidelity.vectorized == 0)
            if (prob.control.verbose > 0)
                fprintf('%s','serial...');
            end
            for i = 1:nx
                result = feval(prob.function.objfun, x(i,:), prob.param);
                f(i,:) = reshape(result,1,mf);
            end
        else
            if (prob.control.verbose > 0)
                fprintf('%s','serial-vectorized...');
            end
            f = feval(prob.function.objfun, x, prob.param);
        end
    else
        if (prob.control.verbose > 0)
            fprintf('%s','parallel...');
        end
        currentpool = gcp('nocreate');
        for i = 1:nx
            fevFuture(i) = parfeval( ...
            	currentpool, prob.function.objfun, 1, x(i,:), prob.param);
        end
        for i = 1:nx
            [completedIdx,value] = fetchNext(fevFuture);
            f(completedIdx,:) = reshape(value,1,mf);
            fprintf('%d,',completedIdx);
        end
    end
    %---------------------------------------------------------------------------
    if (prob.control.verbose == 2)
        fprintf('\n');
        for i = 1:nx
            fprintf('%12.4f',x(i,:));
            fprintf('  /');
            fprintf('%12.4f',f(i,:));
            fprintf('\n');
        end
        fprintf('...');
    end
    if (prob.control.verbose > 0); fprintf('%s\n','done'); end
end
%===============================================================================
