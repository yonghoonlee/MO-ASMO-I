%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Evaluate function values on the surrogate model
%===============================================================================
function fout = surrogateEval(xin, surrogate)
    method = surrogate.method;
    nxin = size(xin,1);
    xlb = surrogate.scale.xlb;
    xub = surrogate.scale.xub;
    flb = surrogate.scale.flb;
    fub = surrogate.scale.fub;
    %---------------------------------------------------------------------------
    % Scale input
    x = (xin - repmat(reshape(xlb,1,numel(xlb)),nxin,1)) ...
        ./(repmat(reshape(xub,1,numel(xub)),nxin,1) ...
            - repmat(reshape(xlb,1,numel(xlb)),nxin,1));
    %---------------------------------------------------------------------------
    switch lower(method)
        case 'rbf'                          % Radial-basis function (deprecated)
            basisfn = surrogate.basisfn;
            epsilon = surrogate.epsilon;
            w = surrogate.w;
            c = surrogate.c;
            nx = size(x,1);
            nc = size(c,1);
            mw = size(w,2);
            f = zeros(nx,mw);
            for idx1 = 1:mw
                phi = zeros(nx,nc);
                for idx2 = 1:nc
                    r = sqrt(sum((repmat(c(idx2,:),nx,1) - x).^2,2));
                    switch lower(basisfn)
                        case 'linear'       % Linear basis fn
                            phi(:,idx2) = r;
                        case 'cubic'        % Cubic basis fn
                            phi(:,idx2) = r.^3;
                        case 'tps'          % Thin plate spline basis fn
                            phi(:,idx2) = r.^2.*log(r);
                            phi(r<eps,idx2) = 0;
                        case 'gaussian'     % Gaussian basis fn
                            phi(:,idx2) = exp(-(epsilon.*r).^2);
                        case 'mq'           % Multiquadric basis fn
                            phi(:,idx2) = sqrt(1+(epsilon.*r).^2);
                        case 'invmq'        % Inverse multiquadric basis fn
                            phi(:,idx2) = 1./sqrt(1+(epsilon.*r).^2);
                        otherwise
                            error(strcat(basisfn,'::not supported.'));
                    end
                end
                f(:,idx1) = phi*w(:,idx1);
            end
            %-------------------------------------------------------------------
        case 'rbn'                          % Radial-basis function network
            rbmodel = surrogate.rbmodel; % Gaussian with epsilon=1
            f = rbmodel(x')';
            %-------------------------------------------------------------------
        case 'snn'                          % Shallow neural network
            snnmodel = surrogate.nnmodel;
            f = snnmodel(x')';
            %-------------------------------------------------------------------
        case 'gpr'
            mf = size(surrogate.gpm,2);
            f = zeros(size(x,1),mf);
            for idx = 1:mf
                f(:,idx) = predict(surrogate.gpm{idx}, x);
            end
            %-------------------------------------------------------------------
        otherwise
            error(strcat(method,'::not supported.'));
    end
    %---------------------------------------------------------------------------
    % Descale
    nf = size(f,1);
    fout = repmat(reshape(flb,1,numel(flb)),nf,1) ...
        + (repmat(reshape(fub,1,numel(fub)),nf,1) ...
            - repmat(reshape(flb,1,numel(flb)),nf,1)).*f;
end
%===============================================================================
