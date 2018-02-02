%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Quarter car suspension multiobjective optimization problem
%===============================================================================
function result = QCsus_CriticalGel651
    casefile = mfilename('fullpath');
    result = runMOASMO(@settings, @obj, @nonlcon, casefile);
end
%===============================================================================
function prob = settings(prob)
    prob.random.seed = 1;
	% x = [log10(S); log10(lambda_0); log10(lambda_max); n]
	prob.bound.xlb = [+0; -9;  0; 1e-9];
    prob.bound.xub = [+6; +2; 20; 1e+0];
    prob.bound.flb = [ 0; 0];
    prob.bound.fub = [50; 1.5];
    prob.gamultiobj.opt.Generations = 100;
    prob.highfidelity.expensive = true;
    prob.highfidelity.vectorized = false;
    prob.lincon.A = [0 1 -1 0];
    prob.lincon.b = [-1e-3];
    prob.lincon.Aeq = [];
    prob.lincon.beq = [];
    [mpath,~] = fileparts(mfilename('fullpath'));
    load(fullfile(mpath,'DoubleSine651.mat'), 'road*');
    % road_x = 0:0.01:100; road_z = sin(t/100*pi)*0.01.*sin(t*pi);
    prob.param.m1 = 325;
    prob.param.m2 = 65;
    prob.param.k1 = 500;
    prob.param.k2 = 232500;
    prob.param.g = 9.806;
    prob.param.v = [5,10]; % Lower velocity values are tested
    prob.param.xlb = prob.bound.xlb;
    prob.param.xub = prob.bound.xub;
    prob.param.road_x = road_x;
    prob.param.road_z = road_z;
    prob.plotpareto.type = 'pareto2d';
    prob.plotpareto.range = [];
    prob.control.plotexport = false;
    prob.control.maxiter = 100;     % maximum number of iterations
    prob.control.maxerror = 1e-4;   % maximum allowable error
    prob.sampling.initnumber = 20;  % initial
    prob.sampling.valnumber = 10;    % validation
    prob.sampling.upnumber = 5;     % exploitation
    prob.sampling.upexpnumber = 5;  % exploration
    prob.surrogate.method = 'GPR';  % regression model
end
%===============================================================================
function f = obj(x,p)
    % Parameters
    m1 = p.m1;                              % 1/4 sprung mass [kg]
    m2 = p.m2;                              % 1/4 unsprung mass [kg]
    k1 = p.k1;                              % Suspension spring constant [N/m]
    k2 = p.k2;                              % Tire stiffness [N/m]
    g = p.g;                                % Gravitational acc [m/s^2]
    vlist = p.v;                            % Vehicle velocity list [m/s]
    road_x = p.road_x;                      % Road profile in x [m]
    road_z = p.road_z;                      % Road profile in z [m]
    %---------------------------------------------------------------------------
    % Decompose design variables
    S = 10^x(1);            % Strength of the critical gel
    lambda_0 = 10^x(2);     % Integrating lower bound
    lambda_max = 10^x(3);   % Integrating upper bound
    n = x(4);               % Deviation of the continuous spectra
    %---------------------------------------------------------------------------
    % Gamma function can be approximated with following formula,
    % gamma_approximation = @(x) 1./x.*(1-0.1138.*(1-4.*(x-0.5).^2));
    % Otherwise, use MATLAB intrinsic gamma function
    %---------------------------------------------------------------------------
    f = [0;0];
    %---------------------------------------------------------------------------
    H = @(tau) S/gamma(n)*tau.^(-n);
    %---------------------------------------------------------------------------
    xlbnd = lambda_0;
    xrbnd = lambda_max;
    if xlbnd >= xrbnd
        f = [NaN; NaN];
        return;
    end
    %---------------------------------------------------------------------------
    for k = 1:length(vlist)
        v = vlist(k);
        % Time span
        t = [road_x/v];                     % t [s]
        dt = diff(t);                       % dt = t_{k+1} - t_{k} [s]
        dt = [dt, dt(end)];
        
        K = zeros(1,size(t,2));
        for i = 1:size(t,2)
            Kinteg = real( ...
                integral( ...
                    @(s) H(exp(s)).*exp(-t(i)./exp(s)),log(xlbnd),log(xrbnd)));
            if abs(Kinteg) > 1e16
                f = [NaN; NaN];
                return;
            end
            K(i) = Kinteg;
        end
        %-----------------------------------------------------------------------
        % Simulation
        xi = zeros(4,size(t,2));                % Create state
        xidot = zeros(4,size(t,2));             % Create derivative of state
        %-----------------------------------------------------------------------
        % Euler predictor-corrector method for simulation
        for i = 1:(size(t,2)-1)
            convint = 0;
            for j = 1:i
                convint = convint ...
                    + K(j)*(xi(3,(i-j+1)) - xi(4,(i-j+1)))*dt(j);
            end
            
            xidot(:,i) = [xi(3,i); xi(4,i);
                (-convint - k1*(xi(1,i) - xi(2,i)))/m1;% - g;
                (convint + k1*(xi(1,i) - xi(2,i)) ...
                - k2*(xi(2,i) - road_z(i)))/m2 - g];
            xist = xi(:,i) + dt(i)*xidot(:,i); % Predictor step (\xi_{*})
            xitp = [xi(:,1:i), xist];
            convintst = 0;
            for j = 1:i+1
                convintst = convintst ...
                    + K(j)*(xitp(3,(i-j+2)) - xitp(4,(i-j+2)))*dt(j);
            end
            fst = [xist(3); xist(4);
                (-convintst - k1*(xist(1) - xist(2)))/m1;% - g;
                (convintst + k1*(xist(1) - xist(2)) ...
                - k2*(xist(2) - road_z(i+1)))/m2 - g];
            xi(:,i+1) = xi(:,i) + dt(i)/2*(xidot(:,i) + fst);
        end
        %-----------------------------------------------------------------------
        % Objective functions
        f = f + [ max(abs(xidot(3,:)));                             % Comfort
                  max(xi(2,:) - road_z) - min(xi(2,:) - road_z) ];  % Handling
        if ~(isreal(f(1)) && isreal(f(2)))
            f = [NaN; NaN];
        end
        if (f(1) < 0.01) && (f(2) < 0.01)
            f = [NaN; NaN];
        end
    end
end
%===============================================================================
function [c,ceq] = nonlcon(x,p)
    c = [];
    ceq = [];
end
%===============================================================================
