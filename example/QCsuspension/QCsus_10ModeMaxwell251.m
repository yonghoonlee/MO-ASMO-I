%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Quarter car suspension multiobjective optimization problem
%===============================================================================
function result = QCsus_10ModeMaxwell251
    casefile = mfilename('fullpath');
    result = runMOASMO(@settings, @obj, @nonlcon, casefile);
end
%===============================================================================
function prob = settings(prob)
    prob.random.seed = 1;
	% x = [Ge; k1; lambda_m; G_m]
	prob.bound.xlb = [   0; 300; -5; -5; -5; -5; -5; -5; -5; -5; -5; -5; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
    prob.bound.xub = [1000; 700;  2;  2;  2;  2;  2;  2;  2;  2;  2;  2; 6; 6; 6; 6; 6; 6; 6; 6; 6; 6];
    prob.bound.flb = [0; 0];
    prob.bound.fub = [50; 1.5];
    prob.gamultiobj.opt.Generations = 100;
    prob.highfidelity.expensive = true;
    prob.highfidelity.vectorized = false;
    prob.lincon.A = [0 0 1 -1  0  0  0  0  0  0  0  0 0 0 0 0 0 0 0 0 0 0;  % lambda(1) <= lambda(2)
                     0 0 0  1 -1  0  0  0  0  0  0  0 0 0 0 0 0 0 0 0 0 0;  % lambda(2) <= lambda(3)
                     0 0 0  0  1 -1  0  0  0  0  0  0 0 0 0 0 0 0 0 0 0 0;  % lambda(3) <= lambda(4)
                     0 0 0  0  0  1 -1  0  0  0  0  0 0 0 0 0 0 0 0 0 0 0;  % lambda(4) <= lambda(5)
                     0 0 0  0  0  0  1 -1  0  0  0  0 0 0 0 0 0 0 0 0 0 0;  % lambda(5) <= lambda(6)
                     0 0 0  0  0  0  0  1 -1  0  0  0 0 0 0 0 0 0 0 0 0 0;  % lambda(6) <= lambda(7)
                     0 0 0  0  0  0  0  0  1 -1  0  0 0 0 0 0 0 0 0 0 0 0;  % lambda(7) <= lambda(8)
                     0 0 0  0  0  0  0  0  0  1 -1  0 0 0 0 0 0 0 0 0 0 0;  % lambda(8) <= lambda(9)
                     0 0 0  0  0  0  0  0  0  0  1 -1 0 0 0 0 0 0 0 0 0 0]; % lambda(9) <= lambda(10)
    prob.lincon.b = [0; 0; 0; 0; 0; 0; 0; 0; 0];
    prob.lincon.Aeq = [];
    prob.lincon.beq = [];
    [mpath,~] = fileparts(mfilename('fullpath'));
    load(fullfile(mpath,'IRI_737b.mat'), 'road*');
    prob.param = struct('m1', 325, 'm2', 65, 'k2', 232500, ...
        'g', 9.806, 'v', [5,10,15,30,45,60,80], ...
        'xlb', prob.bound.xlb, 'xub', prob.bound.xub, ...
        'road_x', road_x, 'road_z', road_z);
    prob.plotpareto.type = 'pareto2d';
    prob.plotpareto.range = [];
    prob.control.maxiter = 30;     % maximum number of iterations
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
    k2 = p.k2;                              % Tire stiffness [N/m]
    g = p.g;                                % Gravitational acc [m/s^2]
    vlist = p.v;                            % Vehicle velocity list [m/s]
    road_x = p.road_x;                      % Road profile in x [m]
    road_z = p.road_z;                      % Road profile in z [m]
    %---------------------------------------------------------------------------
    % Decompose design variables
    %Hmax = 10^x(1);     % Max value of log-normal continuous spectra
    %taumax = 10^x(2);   % Time scale of max in the continuous spectra
    %sigma = x(3);       % Deviation of the continuous spectra
    Ge = x(1);              % Base value in the relaxation kernel
    k1 = x(2);              % Suspension stiffness [N/m]
    lambda_m = 10.^x(3:12);  % Maxwell model relaxation time scale parameter
    K_m = 10.^x(13:22);       % Maxwell model relaxation stress parameter
    %---------------------------------------------------------------------------
    f = [0;0];
    %---------------------------------------------------------------------------
    for k = 1:length(vlist)
        v = vlist(k);
        % Time span
        t = [road_x/v];                     % t [s]
        dt = diff(t);                       % dt = t_{k+1} - t_{k} [s]
        dt = [dt, dt(end)];
        
        K = zeros(1,size(t,2));
        for j = 1:10
            K = K + K_m(j)*exp(-t/lambda_m(j));
        end
        K = K + Ge;
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
