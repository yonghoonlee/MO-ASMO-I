%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Surrogate model construction
%===============================================================================
function surrogate = surrogateConstruct(x, f, prob)
    if (prob.control.verbose > 0)
        fprintf('%s','Constructing surrogate model...');
    end
    %---------------------------------------------------------------------------
    method = prob.surrogate.method;
    switch lower(method)
        case 'rbf'                          % Radial-basis function (deprecated)
            basisfn = prob.surrogate.basisfn;
            epsilon = prob.surrogate.epsilon;
            nx = size(x,1);
            mf = prob.nfvar;
            if (prob.control.verbose > 0)
                fprintf('%s',['using RBF with ', basisfn ,'...']);
            end
            c = x;                          % Set centers using known pts
            w = zeros(nx,mf);               % Pre-allocate weight matrix
            for k = 1:mf                    % For each obj fn variable
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
                w(:,k) = pinv(phi)*f(:,k);
                % w(:,k) = phi\f(:,k);
            end
            surrogate = [];                        % Clear structure
            surrogate.w = w;                       % Save weights
            surrogate.c = c;                       % Save centers
            surrogate.method = method;
            surrogate.basisfn = basisfn;
            surrogate.epsilon = epsilon;
            %-------------------------------------------------------------------
        case 'rbn'
            surrogate.method = method;
            surrogate.rbmodel = newrb(x', f');
            %-------------------------------------------------------------------
        case 'snn'
            hiddenLayerSize = ceil(0.5*prob.nxvar + prob.nfvar);
            net = fitnet(hiddenLayerSize);
            net.divideParam.trainRatio = 8.0;
            net.divideParam.valRatio = 2.0;
            net.divideParam.testRatio = 0.0;
            net.trainFcn = prob.surrogate.snntrainfnc;
            net.trainParam.max_fail = prob.surrogate.snnmaxfail;
            [tr,~] = train(net,x',f');
            surrogate.method = method;
            surrogate.nnmodel = tr;
            %-------------------------------------------------------------------
        case 'gpr'
            mf = prob.nfvar;
            for k = 1:mf
                surrogate.gpm{k} = fitrgp(x,f(:,k));
            end
            surrogate.method = method;
            %-------------------------------------------------------------------
        otherwise
            error(strcat(method,'::not supported.'));
    end
    %---------------------------------------------------------------------------
    if (prob.control.verbose > 0)
        fprintf('%s\n','done');
    end
end
%===============================================================================
