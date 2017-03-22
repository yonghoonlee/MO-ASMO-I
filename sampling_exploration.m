% This function generates samples for exploration of design space.
% Usage: xf = SAMPLING_EXPLORATION(problem,k,(prevPoints))
% Input: problem,k,(prevPoints)
% Output: xf
%   problem: problem definition structure
%   k: if initial sampling, k=0
%   prevPoints: if k>0, prevPoints matrix is required
%   xf: sample points for exploration of design space

% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.

function xf = sampling_exploration(problem,k,varargin)
    if (problem.control.verbose > 0)
        fprintf('Sampling for exploration...');
    end
    if (k == 0)
        method = problem.sampling.initmethod;
        number = problem.sampling.initnumber;
        prevPoints = [];
    else
        method = problem.sampling.upexpmethod;
        number = problem.sampling.upexpnumber;
        if (nargin > 2)
            prevPoints = varargin{1};
        else
            prevPoints = [];
        end
    end
    dimension = problem.nxvar;
    A = problem.A;
    b = problem.b;
    Aeq = problem.Aeq;
    beq = problem.beq;
    lb = problem.xlb;
    ub = problem.xub;
    nonlcon = problem.nonlconfun;
    opt = problem.sampling.initconopt;
    wei1 = problem.sampling.initconobjw;
    wei2 = problem.sampling.initcondispw;
    
    switch lower(method)
        case 'lhs'
            ex = 0; % To enter the while loop
            while (ex < 1)
                ex = 2;
                % Sampling and descaling
                xt = sampling_LHS(number,dimension);
                xt = sampling_descaling(xt,lb,ub);
                % Sample adjustment (to satisfy constraints)
                fprintf('%s','adjust samples to satisfy constraints...');
                xf = xt;
                if (problem.control.verbose == 2); fprintf('\n'); end
                for i = 1:size(xt,1)
                    xtmp = xt(i,:);
                    [xs,f,e] = fmincon( ...
                        @(x)sampling_cobj(x,xtmp,wei1,wei2,prevPoints),...
                        xtmp,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x),opt);
                    xf(i,:) = reshape(xs,1,numel(xs));
                    if (size(prevPoints,1) ~= 0)
                        prevPoints = [prevPoints; reshape(xs,1,numel(xs))];
                    end
                    ex = min(ex,e);
                    if (problem.control.verbose == 2)
                        disp(strcat('exitflag:',num2str(e),...
                                    '/objfn:',num2str(f)));
                    end
                end
            end
        otherwise
            error(strcat(method,'::not supported.'));
    end
    if (problem.control.verbose == 2); fprintf('...'); end
    if (problem.control.verbose > 0)
        fprintf('%s\n','done');
    end
end
