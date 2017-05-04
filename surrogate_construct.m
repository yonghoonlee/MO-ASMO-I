% This function constructs a surrogate model based on training samples.
% Usage: SM = SURROGATE_CONSTRUCT(x,f,problem)
% Input: x,SM
% Output: f
%   x: Training sample points in design space
%   f: Training sample points in objective function space
%   problem: problem definition structure
%   SM: Surrogate model data structure

% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.

function SM = surrogate_construct(x,f,problem)
    if (problem.control.verbose > 0)
        fprintf('%s','Constructing surrogate model...');
    end
    method = problem.surrogate.method;
    switch lower(method)
        case 'rbf'                          % Radial-basis function
            basisfn = problem.surrogate.basisfn;
            epsilon = problem.surrogate.epsilon;
            %SM = RBF(x,f,[],[],'construct',basisfn,epsilon,problem);
            [nx,~] = size(x);
            [~,mf] = size(f);
            if (problem.control.verbose > 0)
                fprintf('%s',['using RBF with ', basisfn ,'...']);
            end
            c = x;                          % Set centers using known pts
            w = zeros(nx,mf);               % Pre-allocate weight matrix
            for i = 1:mf                    % For each obj fn variable
                phi = zeros(nx,nx);         % Gram matrix
                for j = 1:nx                % For each training pts
                    % Radius (distance) computation in vector
                    r = sqrt(sum((repmat(x(j,:),nx,1) - x).^2,2));
                    switch lower(basisfn)
                        case 'linear'       % Linear basis fn
                            phi(:,j) = r;
                        case 'cubic'        % Cubic basis fn
                            phi(:,j) = r.^3;
                        case 'tps'          % Thin plate spline basis fn
                            phi(:,j) = r.^2.*log(r);
                            phi(r<eps,j) = 0;
                        case 'gaussian'     % Gaussian basis fn
                            phi(:,j) = exp(-(epsilon.*r).^2);
                        case 'mq'           % Multiquadric basis fn
                            phi(:,j) = sqrt(1+(epsilon.*r).^2);
                        case 'invmq'        % Inverse multiquadric basis fn
                            phi(:,j) = 1./sqrt(1+(epsilon.*r).^2);
                        otherwise
                            error(strcat(basisfn,'::not supported.'));
                    end
                end
                % w(:,i) = pinv(phi)*f(:,i);
                w(:,i) = phi\f(:,i);
            end
            SM = [];                        % Clear structure
            SM.w = w;                       % Save weights
            SM.c = c;                       % Save centers
            SM.method = method;
            SM.basisfn = basisfn;
            SM.epsilon = epsilon;
        otherwise
            error(strcat(method,'::not supported.'));
    end
    if (problem.control.verbose > 0)
        fprintf('%s\n','done');
    end
end