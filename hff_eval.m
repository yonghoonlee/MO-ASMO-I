%==========================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.
%==========================================================================
% HFF_EVAL FUNCTION
%==========================================================================
% This function runs high fidelity function to evaluate designs
% Input: x, problem
%   x: designs for high fidelity function evaluation
% Output: f
%==========================================================================

function f = hff_eval(x,problem)
    if (problem.control.verbose > 0)
        fprintf('%s','High fidelity function evaluation...');
    end
    nx = size(x,1);                 % Number of designs
    mf = problem.nfvar;             % Number of obj function variables
    f = zeros(nx,mf);               % Allocate result matrix
    npool = hff_getpoolsize();      % Number of parallel pools
    if (problem.control.verbose > 0)
        fprintf('%s',strcat(num2str(size(x,1)),'-points...'));
    end
    if (npool == 0)
        if (problem.highfidelity.vectorized == 0)
            if (problem.control.verbose > 0)
                fprintf('%s','serial...');
            end
            for i = 1:nx
                result = feval(problem.objfun,x(i,:));
                f(i,:) = reshape(result,1,mf);
            end
        else
            if (problem.control.verbose > 0)
                fprintf('%s','serial-vectorized...');
            end
            f = feval(problem.objfun,x);
        end
    else
        if (problem.control.verbose > 0)
            fprintf('%s','parallel...');
        end
        parfor i = 1:nx
            result = feval(problem.objfun,x(i,:));
            f(i,:) = reshape(result,1,mf);
        end
    end
    if (problem.control.verbose > 1)
        fprintf('\n');
        for i = 1:nx
            fprintf('%12.4f',x(i,:));
            fprintf('  /');
            fprintf('%12.4f',f(i,:));
            fprintf('\n');
        end
        fprintf('...');
    end
    if (problem.control.verbose > 0); fprintf('%s\n','done'); end
end