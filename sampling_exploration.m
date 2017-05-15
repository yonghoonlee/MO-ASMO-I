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
    p = problem.p;
    opt = problem.sampling.initconopt;
    wei1 = problem.sampling.initconobjw;
    wei2 = problem.sampling.initcondispw;
    
    switch lower(method)
        case 'lhs'
            ex = -1; % To enter the while loop
            while (ex < 0)
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
                        xtmp,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,p),opt);
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
    
%     %==
%     % For Simionescu Problem Only
%     fgtest = figure('Color',[1 1 1]);
%     angle = linspace(0,2*pi,361);
%     for iangle = 1:length(angle)-1
%         tangle1 = angle(iangle);
%         rangle1 = (1 + 0.2*cos(8*atan(cot(tangle1))));
%         tangle2 = angle(iangle+1);
%         rangle2 = (1 + 0.2*cos(8*atan(cot(tangle2))));
%         plot([rangle1*cos(tangle1);rangle2*cos(tangle2)],...
%              [rangle1*sin(tangle1);rangle2*sin(tangle2)],'r-');
%         hold on;
%     end
%     axis 'square';
%     plot(xt(:,1),xt(:,2),'bo');
%     plot(xf(:,1),xf(:,2),'kx');
%     axis([problem.xlb(1),problem.xub(1),problem.xlb(2),problem.xub(2)]);
%     if (problem.control.plotexport ~= 0)
%         eval(['export_fig ',fullfile(problem.probpath,...
%             ['fig_smpexplore_',num2str(k),'.pdf']), ' -pdf']);
%     end
%     close(fgtest);
%     %==
end
