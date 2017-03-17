function z = multimode_maxwell_multiobjective_obj(x,p)
    % z = [maxAmp, maxAcc, AUC_Amp, AUC_Acc],
    % where
    % maxAmp - maximum amplitude of mass, m1
    % maxAcc - maximum acceleration of mass, m1
    % AUC_Amp - integral of position of mass, m1
    % AUC_Acc - integral of acceleration of mass, m1

    % Set parameters
    nx = size(x,1);
    n = p.n; % number of modes in the multi-mode maxwell model
    m = p.m; % mass of m1
    k = p.k; %  spring constant
    zeta = p.zeta;
    omega1 = p.omega1; % natural frequency of base system
    omega = p.omega;
    Y = p.Y; % amplitude of forcing frequency

    % Define Maxwell element parameters
    eta = zeros(nx,n);
    lambda = zeros(nx,n);
    for i=1:n
        eta(:,i) = 10.^x(:,2*i-1);
        lambda(:,i) = 10.^x(:,2*i);
    end

    % Necessary calculated parameters
    c = zeta*2*m*omega1; % absorber damping coefficient
    C = zeros(nx,length(omega)); %pre-allocating
    S = zeros(nx,length(omega)); %pre-allocating
    X_m = zeros(nx,length(omega)); %pre-allocating
    A_m = zeros(nx,length(omega)); %pre-allocating

    for i = 1:length(omega)
        for q = 1:n
            C(:,i) = C(:,i) + eta(:,q) ./ (1+(lambda(:,q) .* repmat(omega(i),nx,1)).^2);
            S(:,i) = S(:,i) + eta(:,q) .* lambda(:,q) .* repmat(omega(i),nx,1) ./ (1 + (lambda(:,q) .* repmat(omega(i),nx,1)).^2);
        end
        %%%%SPRING AND DASHPOT IN PARALLEL
        a = repmat(m,nx,1).*(repmat(omega(i),nx,1)).^2 - repmat(k,nx,1) - repmat(omega(i),nx,1).*S(:,i);
        b = repmat(omega(i),nx,1).*(repmat(c,nx,1) + C(:,i));

        for j = 1:nx
            A = [a(j) b(j); -b(j) a(j)];
            B = Y*[-omega(i)*S(j,i)-k; omega(i)*(-C(j,i)-c)];
            X = A\B;
            X_m(j,i) = sqrt(X(1)^2 + X(2)^2); % Amplitude of mass, m
            A_m(j,i) = X_m(j,i)*omega(i).^2; % Acceleration function
        end
    end
    
    % Function outputs
    X = X_m/Y; % nondimensional displacement vector of m
    
    % Finding a relevant natural frequency for the system
    %[~, Index] = max(X,[],2); % Finds the index of the maximum amplitude
    %omega_max = omega(Index); %Selects the omega at this maximum value from the original omega vector
    %omega_shifted = omega/omega_max;
    
    A = A_m; % Nondimensional acceleration
    
    maxAmp = max(X,[],2); % maximum amplitude of m
    maxAcc = max(A,[],2); % maximum accleration of m

    AUC_Amp = trapz(X,2); % area under displacement curve
    AUC_Acc = trapz(A,2); % area under acceleartion curve
    
    z = [maxAmp, maxAcc, AUC_Amp, AUC_Acc];

end