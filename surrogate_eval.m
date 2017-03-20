% This function evaluates surrogate model for given input.
% Usage: f = SURROGATE_EVAL(x,SM)
% Input: x,SM
% Output: f
%   x: Input design points for evaluating surrogate model
%   SM: Surrogate model data structure
%   f: Predicted objective function value for given x using surrogate model

% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.

function f = surrogate_eval(x,SM)
    switch lower(SM.method)
        case 'rbf'
            w = SM.w;
            c = SM.c;
            basisfn = SM.basisfn;
            epsilon = SM.epsilon;
            [nx,~] = size(x);
            [nc,~] = size(c);
            [~,mw] = size(w);
            f = zeros(nx,mw);
            for i = 1:mw
                phi = zeros(nx,nc);
                for j = 1:nc
                    r = sqrt(sum((repmat(c(j,:),nx,1) - x).^2,2));
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
                f(:,i) = phi*w(:,i);
            end
        otherwise
            error(strcat(SM.method,'::not supported.'));
    end
end