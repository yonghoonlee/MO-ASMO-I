%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Texture and fluid design problem using 3D Cauchy momentum equation solver with
% 2-mode Giesekus model.
%===============================================================================
function result = CEF2_case2531
    casefile = mfilename('fullpath');
    result = runMOASMO(@settings, @obj, @nonlcon, casefile);
end
%===============================================================================
function prob = settings(prob)
    N       =      5;
    hmin    =    269.0e-6;
    hmax    =    800.0e-6;
    eta     =      9.624e-3;
    rho     =    873.4;
    Omega   =     10.0;
    Ntex    =     10;
    ri      =      0.5e-3;
    ro      =     20.0e-3;
    tan_angle = sqrt(3); % tan(60deg) = sqrt(3), mesh control point slope
    xlb_geom = hmin*ones(((N+1)*N),1);
    xub_geom = hmax*ones(((N+1)*N),1);
    nmode = 2;
    xlb_fld = [0*ones(nmode,1);       % eta_p_i
               0.00001*ones(nmode,1); % lambda_i
               0.01*ones(nmode,1)];   % alpha_i
    xub_fld = [5/2*eta*ones(nmode,1); % eta_p_i <= 5/2*eta_s (Explained below)
               0.01*ones(nmode,1);    % lambda_i
               0.5*ones(nmode,1)];    % alpha_i <= 0.5 (Based on limits of Giesekus model)
    % limiting upper limit of eta_p_i <= 5/2*eta_s
    % : This has to do with a limit that we?ve been looking at where we can
    % relate the polymer viscosity to the concentration of the polymer in
    % solution. In order to do that, we need the solution to be dilute to
    % semi-dilute, which limits the concentration, and hence limits the max
    % value of eta_p.
    %---------------------------------------------------------------------------
    prob.bound.xlb = [xlb_geom; xlb_fld];
    prob.bound.xub = [xub_geom; xub_fld];
    prob.bound.flb = [1e-4; -0.1];
    prob.bound.fub = [1e-2;  0.0];
    prob.control.plotexport = false;
    prob.control.maxiter = 100;     % maximum number of iterations
    prob.control.maxerror = 1e-8;   % maximum allowable error
    prob.gamultiobj.opt.Generations = 50;
    prob.highfidelity.expensive = true;
    prob.highfidelity.vectorized = false;
    prob.highfidelity.infeaparam.C = 0.5;
    prob.highfidelity.infeaparam.q = 0.1;
    prob.plotpareto.type = 'pareto2d';
    prob.plotpareto.range = [];
    prob.sampling.initnumber = 36;  % initial
    prob.sampling.valnumber = 12;    % validation
    prob.sampling.upnumber = 12;    % exploitation
    prob.sampling.upexpnumber = 12;  % exploration
    prob.surrogate.method = 'GPR';  % regression model
    %---------------------------------------------------------------------------
    % Linear Constraints
    A1 = []; b1 = [];
    % Specify lowest point at the outer-most arc.
    for idx = 2:(N) % (Nth-1) rows in sparse constraint vector
        A1 = [A1; idx-1, N+1, -1; idx-1, idx*(N+1), 1];
        b1 = [b1; 0];
    end
    % Specify slope constraints
    [~,~,~,~,z,~] = semhat(N);
    r = (1+z)/2*ro+(1-z)/2*ri;
    theta = (2*pi/Ntex)/2*z;
    [Rmat,Theta] = ndgrid(r,theta);
    X2d = Rmat.*cos(Theta);
    Y2d = Rmat.*sin(Theta);
    cdx = N-1;
    for jdx = 1:(N) % 4*(Nth*Nr) rows in sparse constraint vector
        for idx = 1:(N)
            kdx = (jdx-1)*(N+1)+idx;
            kdxpr = kdx + 1;
            kdxpth = kdx + (N+1);
            if kdxpth > (N+1)*N
                kdxpth = kdxpth - (N+1)*N;
            end
            dr = sqrt((X2d(idx+1,jdx) - X2d(idx,jdx))^2 + (Y2d(idx+1,jdx) - Y2d(idx,jdx))^2);
            dth = sqrt((X2d(idx,jdx+1) - X2d(idx,jdx))^2 + (Y2d(idx,jdx+1) - Y2d(idx,jdx))^2);
            cdx = cdx + 1;
            A1 = [A1; cdx, kdx, -1; cdx, kdxpr, 1];
            b1 = [b1; tan_angle*dr];
            cdx = cdx + 1;
            A1 = [A1; cdx, kdx, 1; cdx, kdxpr, -1];
            b1 = [b1; tan_angle*dr];
            cdx = cdx + 1;
            A1 = [A1; cdx, kdx, -1; cdx, kdxpth, 1];
            b1 = [b1; tan_angle*dth];
            cdx = cdx + 1;
            A1 = [A1; cdx, kdx, 1; cdx, kdxpth, -1];
            b1 = [b1; tan_angle*dth];
        end
    end
    A1 = [A1; cdx+1, length(xlb_geom)+1, -1; cdx+1, length(xlb_geom)+2, 1];
    A1 = [A1; cdx+2, length(xlb_geom)+3, -1; cdx+2, length(xlb_geom)+4, 1];
    A1 = [A1; cdx+3, length(xlb_geom)+5, -1; cdx+3, length(xlb_geom)+6, 1];
    b1 = [b1; zeros(3,1)];
    prob.lincon.A = full(sparse(A1(:,1),A1(:,2),A1(:,3),...
        (N-1)+(4*N*N)+3,length(prob.bound.xlb)));
    prob.lincon.b = full(sparse(b1));
    prob.lincon.Aeq = [];
    prob.lincon.beq = [];
    %---------------------------------------------------------------------------
    prob.param = struct('N', N, 'hmin', hmin, 'hmax', hmax, 'eta', eta, 'rho', rho, ...
        'Omega', Omega, 'Ntex', Ntex, 'ri', ri, 'ro', ro, 'nmode', nmode);
end
%===============================================================================
function f = obj(x,param)
    % Operating condition
    Omega = param.Omega;
    %---------------------------------------------------------------------------
    % Geometry
    Ntex = param.Ntex;
    phi = 2*pi/Ntex;
    R1 = param.ri;
    R2 = param.ro;
    hmin = param.hmin;
    %---------------------------------------------------------------------------
    % Mesh
    N = param.N;
    x_geom = x(1:((N+1)*N));
    x_fluid = x(((N+1)*N+1):end);
    %---------------------------------------------------------------------------
    % Surface Geometry
    % 2.69e-6 <= H(i,j) <= Hmax
    H = reshape(x_geom,N+1,N);
    H = [H, H(:,1)];
    %---------------------------------------------------------------------------
    % 1-D matrices for the pseudospectral method
    [Kh,Mh,Ch,Dh,z,w] = semhat(N);
    q = zeros(1,N);
    q(N) = 1;
    Q = [q;eye(N)];
    I = eye(N+1);
    Res = I(2:N,:);
    Pro = Res';
    Nuem = I(1:N,:)';
    r = (R2-R1)/2*(z)+(R2+R1)/2;
    theta = phi/2*(z);
    [Rmat,Theta] = ndgrid(r,theta);
    Rdiff=R2-R1;
    R1d = diag(r);
    R = sparse(kron(I,R1d));
    Rinv1d = diag(1./r);
    Rinv = sparse(kron(I,Rinv1d));
    r2d = diag(R);
    M2d = sparse(kron(Mh,Mh));
    B = sparse(kron(Q,Nuem));
    B2 = sparse(kron(Q,Nuem));
    Hmat = H;
    dh = zeros((N+1)^2,1);
    h = zeros((N+1)^2,1);
    for j = 1:(N+1)
        h((j-1)*(N+1)+1:j*(N+1)) = H(:,j);
    end
    Dr = kron(I,Dh);
    Dt = kron(Dh,I);
    dhdr = (2/Rdiff)*Dr*h;
    dhdt = (2/phi)*Dt*h;
    %---------------------------------------------------------------------------
    % Fluid properties
    gdot = r2d.*Omega./h;
    nmode = param.nmode;
    etas = param.eta;
    rho = param.rho;
    psi1 = zeros((N+1).^2,1);
    eta = etas*ones((N+1).^2,1);
    for j = 1:nmode
        etap = x_fluid(j);
        lam = x_fluid(nmode+j);
        al = x_fluid(2*nmode+j);
        Wi = lam*gdot;
        fx = 8.0*al*(1-al)*Wi.^2;
        chis = zeros((N+1).^2,1);
        for l = 1:((N+1)^2)
            if fx(l) <= 0.35
                chis(l) = 1.0-0.5*fx(l)+0.5*fx(l).^2-(5.0/8.0)*fx(l).^3 ...
                    +(7.0/8.0)*fx(l).^4-(21.0/16.0)*fx(l).^5 ...
                    +(33.0/16.0)*fx(l).^6-(429.0/128.0)*fx(l).^7 ...
                    +(715.0/128.0)*fx(l).^8;
            else
                chis(l) = (sqrt(1+2*fx(l))-1)./(fx(l));
            end
        end
        chi = sqrt(chis);
        f = (1-chi)./(1+(1-2*al)*chi);
        eta = eta + etap*((1.0-f).^2./(1.0+(1.0-2.0*al)*f));
        psi1 = psi1 + 2*etap*lam*((f.*(1.0-al*f))./(al*(1.0-f).*(Wi.^2)));
    end
    %---------------------------------------------------------------------------
    % Calculation
    ETA=sparse(diag(1./(12*eta)));
    H=sparse(diag(h.^3));
    Kr=-sparse((phi/Rdiff)*(Dr'*(M2d*R*H*ETA)*Dr));
    Kt=-sparse((Rdiff/phi)*(Dt'*(M2d*Rinv*H*ETA)*Dt));
    Kb=(Kr+Kt);
    fb=((Rdiff*phi)/4)*M2d*R*(Omega/2*dhdt);

    K=B'*Kb*B; f=B'*fb;
    p=K\f;

    P=B*p;

    P0=zeros((N+1),1);
    for j=1:N+1,
        P0(j)=P((N+1)+(j-1)*(N+1));
    end
    Pref=-0.5*(w'*P0);
    P=P+Pref;
    P1=P;

    dPdth=(2.0/phi)*Dt*P;
    dPdr=(2.0/Rdiff)*Dr*P;
    d2Pdrdth=(2.0/Rdiff)*Dr*(dPdth);
    d2Pdt2=(2.0/phi)*Dt*dPdth;
    d2Pdr2=(2.0/Rdiff)*Dr*dPdr;

    detadr=(2.0/Rdiff)*Dr*eta;
    detadt=(2.0/phi)*Dt*eta;

    br=1.0+(1.0/28.0)*(1.0./(eta*Omega).*dPdth.*(h./r2d).^2).^2 ...
          -(1.0/3.0)*(1.0./(eta*Omega).*dPdth.*(h./r2d).^2);
    brp=-((psi1).*Omega./eta)...
        -(1/20)*((psi1)*Omega./eta).*(1./(eta*Omega).*dPdth.*(h./r2d).^2).^2 ...
        +(1/20)*((psi1)*Omega./eta).*(1./(eta*Omega).*r2d.*dPdr.*(h./r2d).^2).^2;
    gr=(rho*Omega.*h.^2./eta).*(r2d*Omega.*h./40.0).*br...
        +(r2d*Omega.*h./12.0).*brp;
    corr=(1.0./r2d).*((2.0/Rdiff).*(Dr*(r2d.*gr)));

    bt=(1./(eta*Omega).*r2d.*dPdr).*(h./r2d).^2 ...
      -(3./14)*(1./(eta*Omega).^2).*(r2d.*dPdr.*dPdth).*(h./r2d).^4;
    btp=((psi1)*Omega./eta).*(1./(eta*Omega).^2.*(r2d.*dPdr.*dPdth).*(h./r2d).^4);
    gt=(rho*Omega*h.^2./eta).*(r2d*Omega.*h./240.0).*bt ...
      +(r2d*Omega.*h./120).*btp;
    corrt=(1.0./r2d).*((2.0/phi)*(Dt*(gt)));

    fb=(Rdiff*phi/4.0)*M2d*R*((corr+corrt));

    K=transpose(B2)*Kb*B2;
    f=transpose(B2)*fb;

    p2=K\f;
    P2=B2*p2;

    P0=zeros((N+1),1);
    for j=1:N+1,
        P0(j)=P2((N+1)+(j-1)*(N+1));
    end
    Pref=-0.5*(w'*P0);
    P2=P2+Pref;

    P=P1+P2;
    P0=zeros((N+1),1);
    for j=1:N+1,
        P0(j)=P((N+1)+(j-1)*(N+1));
    end
    Pref=-0.5*(w'*P0);
    P=P+Pref;

    dP2dr=(2/Rdiff)*Dr*P2;
    dP2dt=(2/phi)*Dt*P2;
    dP1dt=(2.0/phi)*Dt*P1;
    dP1dr=(2.0/Rdiff)*Dr*P1;
    dPdth=(2.0/phi)*Dt*P;
    dPdr=(2.0/Rdiff)*Dr*P;

    dutdznd=1.0+(1./(2.0*eta*Omega).*dPdth.*(h./r2d).^2)...
               -(1/40).*(rho*Omega*h.^2./eta).*(1./(eta*Omega).*(r2d.*dP1dr).*(h./r2d).^2)...
               +(1/240).*(rho*Omega*h.^2./eta).*(1./(eta*Omega).^2.*(r2d.*dP1dr.*dP1dt).*(h./r2d).^4)...
               -(1/6)*(Omega.*(psi1)./eta).*((1./(eta*Omega).*r2d.*dP1dr).*(h./r2d).^2)...
               -(1/12)*(Omega.*(psi1)./eta).*(1./(eta*Omega).^2.*(r2d.*dP1dr.*dP1dt).*(h./r2d).^4);
    tau=eta.*(r2d*Omega./h).*dutdznd;

    Fn = Ntex*(Rdiff*phi/4.0)*((kron(w,w))'*R*(P));
    M = Ntex*((Rdiff*phi)/4)*((kron(w,w))'*R*(R*tau));
    
    
    % Objective functions
    PowerInput = M*Omega; % PowerInput = Torque * Omega
    %TauStar = (2/(pi*ro^3)*M)/(eta*(Omega*ro)/hmin); % Normalized Apparent Viscosity
    InertiaEffect = -3/40*pi*R2^2*rho*(R2*Omega)^2; % Inertia effect should be considered for computing normal force generation
    NormalForce = Fn - InertiaEffect;
    
    % Objective functions
    f = [PowerInput; -NormalForce];

end
%===============================================================================
function [c,ceq] = nonlcon(x,param)
    c = [];
    ceq = [];
end
%===============================================================================
function [semhatAh,semhatBh,semhatCh,semhatDh,semhatz,semhatw] = semhat(semhatN)
    [semhatz,semhatw] = zwgll(semhatN); %  z = [z0,z1,...,zN] Gauss-Lobatto points 
                                        %  w = [w0,w1,...,wN] Gauss-Lobatto weights 
    semhatBh    = diag(semhatw);        %  1-D mass matrix
    semhatDh    = dhat(semhatz);        %  1-D derivative matrix
    semhatAh    = semhatDh'*semhatBh*semhatDh;  %  1-D stiffness matrix
    semhatCh    = semhatBh*semhatDh;            %  1-D convection operator
end
%===============================================================================
function [zwgllz,zwgllw] = zwgll(zwgllp)
    zwglln = zwgllp+1;
    zwgllz = zeros(1,zwglln);
    zwgllw = zeros(1,zwglln);
    zwgllz(1)=-1;
    zwgllz(zwglln)= 1;
    if zwgllp>1
        if zwgllp==2 
            zwgllz(2)=0; 
        else
            zwgllM=zeros(zwgllp-1,zwgllp-1);
            for zwglli=1:zwgllp-2
                zwgllM(zwglli,zwglli+1)=(1/2)*sqrt((zwglli*(zwglli+2))/((zwglli+1/2)*(zwglli+3/2)));
                zwgllM(zwglli+1,zwglli)=zwgllM(zwglli,zwglli+1);
            end
            zwgllz(2:zwgllp)=sort(real(eig(zwgllM)));
        end
    end
    %compute the weights w
    zwgllw(1)=2/(zwgllp*(zwglln));
    zwgllw(zwglln)=zwgllw(1);
    for zwglli=2:zwgllp

        zwgllx=zwgllz(zwglli);
        zwgllz0=1;
        zwgllz1=zwgllx;
        zwgllz2=0;
        for zwgllj=1:zwgllp-1
            zwgllz2=zwgllx.*zwgllz1*(2*zwgllj+1)/(zwgllj+1)-zwgllz0*zwgllj/(zwgllj+1);
            zwgllz0=zwgllz1;
            zwgllz1=zwgllz2;
        end
        zwgllw(zwglli)=2/(zwgllp*(zwglln)*zwgllz2*zwgllz2);
    end
    zwgllz=zwgllz';
    zwgllw=zwgllw';
end
%===============================================================================
function[fdwfc] = fd_weights_full(fdwfxx,fdwfx,fdwfm)
    fdwfn1 = length(fdwfx);
    fdwfn  = fdwfn1-1;
    fdwfm1 = fdwfm+1;
    fdwfc1       = 1.;
    fdwfc4       = fdwfx(1) - fdwfxx;
    fdwfc = zeros(fdwfn1,fdwfm1);
    fdwfc(1,1) = 1.;
    for fdwfi=1:fdwfn
        fdwfi1  = fdwfi+1;
        fdwfmn = min(fdwfi,fdwfm);
        fdwfc2 = 1.;
        fdwfc5 = fdwfc4;
        fdwfc4 = fdwfx(fdwfi1)-fdwfxx;
        for fdwfj=0:fdwfi-1;              fdwfj1  = fdwfj+1;
            fdwfc3 = fdwfx(fdwfi1)-fdwfx(fdwfj1);
            fdwfc2 = fdwfc2*fdwfc3;
            for fdwfk=fdwfmn:-1:1;         fdwfk1  = fdwfk+1;
                fdwfc(fdwfi1,fdwfk1) = fdwfc1*(fdwfk*fdwfc(fdwfi1-1,fdwfk1-1)-fdwfc5*fdwfc(fdwfi1-1,fdwfk1))/fdwfc2;
            end
            fdwfc(fdwfi1,1) = -fdwfc1*fdwfc5*fdwfc(fdwfi1-1,1)/fdwfc2;
            for fdwfk=fdwfmn:-1:1;         fdwfk1  = fdwfk+1;
                fdwfc(fdwfj1,fdwfk1) = (fdwfc4*fdwfc(fdwfj1,fdwfk1)-fdwfk*fdwfc(fdwfj1,fdwfk1-1))/fdwfc3;
            end
            fdwfc(fdwfj1,1) = fdwfc4*fdwfc(fdwfj1,1)/fdwfc3;
        end
        fdwfc1 = fdwfc2;
    end
end
%===============================================================================
function [dhatDh] = dhat(dhatx)
    dhatn1 = length(dhatx);
    %dhatw  = zeros(dhatn1,2);
    dhatDh = zeros(dhatn1,dhatn1);
    for dhati=1:dhatn1
        dhatw = fd_weights_full(dhatx(dhati),dhatx,1);  % Bengt Fornberg's interpolation algorithm
        dhatDh(:,dhati) = dhatw(:,2);                   % Derivative of pn is stored in column 2.
    end
    dhatDh = dhatDh';
end
%===============================================================================
