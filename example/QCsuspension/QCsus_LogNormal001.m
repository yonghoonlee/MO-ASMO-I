%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Quarter car suspension multiobjective optimization problem
%===============================================================================
function result = QCsus_LogNormal001
    casefile = mfilename('fullpath');
    result = runMOASMO(@settings, @obj, @nonlcon, casefile);
end
%===============================================================================
function prob = settings(prob)
	% x = [log10(Hmax); log10(taumax); sigma; Ge; k1]
	prob.bound.xlb = [0; -5; 1e-5;    0; 300];
    prob.bound.xub = [6;  2;    5; 1000; 700];
    prob.bound.flb = [0; 1.1];
    prob.bound.fub = [4; 1.5];
    prob.gamultiobj.opt.Generations = 50;
    prob.highfidelity.expensive = true;
    prob.highfidelity.vectorized = false;
    prob.lincon.A = [];
    prob.lincon.b = [];
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
    prob.control.maxiter = 100;     % maximum number of iterations
    prob.control.maxerror = 1e-6;   % maximum allowable error
    prob.sampling.initnumber = 10;  % initial
    prob.sampling.valnumber = 5;    % validation
    prob.sampling.upnumber = 3;     % exploitation
    prob.sampling.upexpnumber = 2;  % exploration
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
    Hmax = 10^x(1);     % Max value of log-normal continuous spectra
    taumax = 10^x(2);   % Time scale of max in the continuous spectra
    sigma = x(3);       % Deviation of the continuous spectra
    Ge = x(4);          % Base value in the relaxation kernel
    k1 = x(5);          % Suspension stiffness [N/m]
    %---------------------------------------------------------------------------
    f = [0;0];
    %---------------------------------------------------------------------------
    H = @(tau) Hmax*exp(-0.5*(log(tau/taumax)/sigma).^2);
    Hcutoff = @(tau) log10(H(tau)) + 12; % cutoff at 1e-12
    %---------------------------------------------------------------------------
    osol = optimoptions('fsolve');
    osol.Display = 'off';
    osol.FiniteDifferenceType = 'central';
    osol.FunctionTolerance = 1e-3;
    osol.OptimalityTolerance = 1e-3;
    osol.StepTolerance = 1e-12;
    osol.TypicalX = taumax;
    xlbnd = real(fsolve(Hcutoff,taumax/exp(1)^sigma,osol));
    xrbnd = real(fsolve(Hcutoff,taumax*exp(1)^sigma,osol));
    if (xlbnd < 1e-8) && (xrbnd > 1e-8)
        xlbnd = 1e-8;
    elseif (xlbnd < 1e-8) && (xrbnd > 1e-10)
        xlbnd = 1e-16;
    elseif (xlbnd < 1e-8) && (xrbnd <= 1e-10)
        xlbnd = 1e-32;
    end
    %figure(); fplot(H,[xlbnd,xrbnd]);
    if (xlbnd >= xrbnd)
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
            K(i) = Ge + Kinteg;
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
    %---------------------------------------------------------------------------
    % Decompose design variables
    HmaxNLC = 10^x(1);     % Max value of log-normal continuous spectra
    taumaxNLC = 10^x(2);   % Time scale of max in the continuous spectra
    sigmaNLC = x(3);       % Deviation of the continuous spectra
    GeNLC = x(4);          % Base value in the relaxation kernel
    %k1NLC = x(5);          % Suspension stiffness [N/m]
    %---------------------------------------------------------------------------
    HNLC = @(tau) HmaxNLC*exp(-0.5*(log(tau/taumaxNLC)/sigmaNLC).^2);
    HcutoffNLC = @(tau) log10(HNLC(tau)) + 12; % cutoff at 1e-12
    %---------------------------------------------------------------------------
    osolNLC = optimoptions('fsolve');
    osolNLC.Display = 'off';
    osolNLC.FiniteDifferenceType = 'central';
    osolNLC.FunctionTolerance = 1e-3;
    osolNLC.OptimalityTolerance = 1e-3;
    osolNLC.StepTolerance = 1e-12;
    osolNLC.TypicalX = taumaxNLC;
    xlbndNLC = real(fsolve(HcutoffNLC,taumaxNLC/exp(1)^sigmaNLC,osolNLC));
    xrbndNLC = real(fsolve(HcutoffNLC,taumaxNLC*exp(1)^sigmaNLC,osolNLC));
    if (xlbndNLC < 1e-8) && (xrbndNLC > 1e-8)
        xlbndNLC = 1e-8;
    elseif (xlbndNLC < 1e-8) && (xrbndNLC > 1e-10)
        xlbndNLC = 1e-16;
    elseif (xlbndNLC < 1e-8) && (xrbndNLC <= 1e-10)
        xlbndNLC = 1e-32;
    end
    %figure(); fplot(H,[xlbnd,xrbnd]);
    if (xlbndNLC >= xrbndNLC)
        c = [NaN];
        return;
    end
    %---------------------------------------------------------------------------
    Kmax = GeNLC + real(integral(@(s) HNLC(exp(s)),log(xlbndNLC),log(xrbndNLC)));
    if abs(Kmax) > 1e16
    	c = [Kmax/1e16];
    end
    c = Kmax - p.xub(end);
end
%===============================================================================
