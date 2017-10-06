function f = obj(x,p)

    % Parameters
    m1 = p.m1;                              % 1/4 sprung mass [kg]
    m2 = p.m2;                              % 1/4 unsprung mass [kg]
    k2 = p.k2;                              % Tire stiffness [N/m]
    g = p.g;                                % Gravitational acc [m/s^2]
    vlist = p.v;                            % Vehicle velocity list [m/s]
    road_x = p.road_x;                      % Road profile in x [m]
    road_z = p.road_z;                      % Road profile in z [m]

    % Decompose design variables
    Hmax = 10^x(1);
    taumax = 10^x(2);
    sigma = 10^x(3);
    Ge = x(4);
    k1 = 10^x(5);                           % Suspension stiffness [N/m]

    f = [0;0];

    H = @(tau) Hmax*exp(-0.5*(log(tau/taumax)/sigma).^2);
    Hcutoff = @(tau) log10(H(tau)) + 12; % cutoff at 1e-12
    deltar = 1e6;
    noinf = false;
    while ~noinf
        test = real(Hcutoff(taumax+deltar));
        if ((test > 1e6) || (test < -1e6))
            deltar = deltar/2;
        else
            noinf = true;
        end
        deltar = deltar/10;
    end
    deltal = 1e6;
    noinf = false;
    while ~noinf
        test = real(Hcutoff(taumax-deltal));
        if ((test > 1e6) || (test < -1e6))
            deltal = deltal/2;
        else
            noinf = true;
        end
        deltal = deltal/10;
    end
    resultr = Hcutoff(taumax+deltar);
    resultl = Hcutoff(taumax-deltal);
    if (min(resultr,resultl) < -1e6) || (max(resultr,resultl) > 1e6)
        f = [NaN; NaN];
        return;
    end
    osol = optimoptions('fsolve');
    osol.Display = 'off';
    osol.FiniteDifferenceType = 'central';
    osol.FunctionTolerance = 1e-8;
    osol.OptimalityTolerance = 1e-8;
    osol.StepTolerance = 1e-8;
    osol.TypicalX = taumax;
    xlbnd = real(fsolve(Hcutoff,taumax-deltal,osol));
    xrbnd = real(fsolve(Hcutoff,taumax+deltar,osol));
    if xlbnd<1e-8
        xlbnd = 1e-8;
    end
    if (xlbnd >= xrbnd)
        f = [NaN; NaN];
        return;
    end

    for k = 1:length(vlist)
        v = vlist(k);
        % Time span
        t = [road_x/v];                     % t [s]
        dt = diff(t);                       % dt = t_{k+1} - t_{k} [s]
        dt = [dt, dt(end)];
        
        K = zeros(1,size(t,2));
        for i = 1:size(t,2)
            Kinteg = real(integral(@(s) H(exp(s)).*exp(-t(i)./exp(s)),log(xlbnd),log(xrbnd)));
            if abs(Kinteg) > 1e16
                f = [NaN; NaN];
                return;
            end
            K(i) = Ge + Kinteg;
        end
        %subplot(2,1,2); loglog(t+1e-5,K);
        
        % Simulation
        xi = zeros(4,size(t,2));                % Create state
        xidot = zeros(4,size(t,2));             % Create derivative of state

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

        % Objective functions
        f = f + [   max(abs(xidot(3,:)));                               % Comfort
                    max(xi(2,:) - road_z) - min(xi(2,:) - road_z)   ];  % Handling
        if ~(isreal(f(1)) && isreal(f(2)))
            f = [NaN; NaN];
        end
    end
end