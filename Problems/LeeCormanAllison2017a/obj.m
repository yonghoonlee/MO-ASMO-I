% Objective function of Lee-Corman-Allison viscoelastic suspension problem.
% * Lee, Corman, Ewoldt, and Allison, "A Multiobjective Adaptive Surrogate
% Modeling-Based Optimization (MO-ASMO) Framework Using Efficient Sampling
% Strategies," ASME 2017 IDETC/CIE, Cleveland, OH, 2017, to appear.
% * Corman, Rao, Bharadwaj, Allison, and Ewoldt, "Setting Material Function
% Design Targets for Linear Viscoelastic Materials and Structures," Journal
% of Mechanical Design, 138(5), p.051402(12pp), doi: 10.1115/1.4032698.
% * Allison, Guo, Han, "Co-Design of an Active Suspension Using
% Simultaneous Dynamic Optimization," Journal of Mechanical Design, 136(8),
% p.081003(14pp), doi: 10.1115/1.4027335.
%
% Usage: f = OBJ(x)
% Input: x
% Output: f
%   x: Points in design space
%   f: Objective function value of given design points

% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.

%% Quarter car suspension model with viscoelastic damper
%  +++++++++++++++++++++++++++
%  +                 ^ x_1   +
%  +   +---------+   |       +
%  +   |   m_1   |---+       +
%  +   +-+-----+-+           +
%  +     |     |             +
%  +    k_1   K(s)           +
%  +     |     |     ^ x_2   +
%  +   +-+-----+-+   |       +
%  +   |   m_2   |---+       +
%  +   +----+----+           +
%  +        |                +
%  +       k_2       ^ z     +
%  +        |        |       +
%  +   -----+--------+       +
%  +++++++++++++++++++++++++++
%
% $$ m_1 \ddot{x}_1 + \int_{0}^{t} K\left(s\right) \left[\dot{x}_1 \left(t
% - s\right) - \dot{x}_2 \left(t - s\right) \right] ds + k_1 \left(x_1 -
% x_2\right) + m_1 g = 0 $$
%
% $$ m_2 \ddot{x}_2 - \int_{0}^{t} K\left(s\right) \left[\dot{x}_1 \left(t
% - s\right) - \dot{x}_2 \left(t - s\right) \right] ds - k_1 \left(x_1 -
% x_2\right) + k_2 \left(x_2 - z\right) + m_2 g = 0 $$
%
% $$ K \left(t\right) 
% = \sum_{m=1}^{M} K_{m} \exp\left(-t / \lambda_{m}\right) $$
%
% $$ \underline{\xi} =
% \left[\xi_1,\; \xi_2,\; \xi_3,\; \xi_4\right]^{T} =
% \left[x_1,\; x_2,\; \dot{x}_1,\; \dot{x}_2\right]^{T} $$
%
% $$ \dot{\underline{\xi}} =
% \left[\dot{\xi}_1,\; \dot{\xi}_2,\; \dot{\xi}_3,\; \dot{\xi}_4\right]^{T}
% = \left[\dot{x}_1,\; \dot{x}_2,\; \ddot{x}_1,\; \ddot{x}_2\right]^{T} $$
%
% $$ \ddot{x}_1 = \left[
% -\int_{0}^{t}{ K\left(s\right) \left\{ \dot{x}_1 \left(t - s\right)
% - \dot{x}_2 \left(t - s\right) \right\} } ds - k_1 \left(x_1 - x_2\right)
% \right] / m_1 - g $$
%
% $$ \ddot{x}_2 = \left[
% \int_{0}^{t}{ K\left(s\right) \left\{ \dot{x}_1 \left(t - s\right)
% - \dot{x}_2 \left(t - s\right) \right\} } ds + k_1 \left(x_1 - x_2\right)
% - k_2 \left(x_2 - z\right) \right] / m_2 - g $$
%
% $$ \underline{\xi}^{k+1} = \underline{\xi}^{k} + \Delta t^{k}/2
% \left(  \right) $$

%% Objective function
function f = obj(x,p)

    % Log scale to linear scale
    xlinear = 10.^x;

    % Parameters
    m1 = p.m1;                              % 1/4 sprung mass [kg]
    m2 = p.m2;                              % 1/4 unsprung mass [kg]
    k2 = p.k2;                              % Tire stiffness [N/m]
    g = p.g;                                % Gravitational acc [m/s^2]
    v = p.v;                                % Vehicle velocity=60mph [m/s]
    road_x = p.road_x;                      % Road profile in x [m]
    road_z = p.road_z;                      % Road profile in z [m]
    
    % Time span
    t = [road_x/v];                         % t [s]
    dt = diff(t);                           % dt = t_{k+1} - t_{k} [s]
    dt = [dt, dt(end)];
    
    % Decompose design variables
    n_m = floor(length(xlinear)/2);         % Number of modes
    K_m = xlinear(1:n_m);                   % K_m
    lambda_m = xlinear((n_m + 1):(2*n_m));  % lambda_m
    k1 = xlinear(end);                      % Suspension stiffness [N/m]
    K = zeros(1,size(t,2));                 % Initialize kernel
    for i = 1:n_m                           % Multimode Maxwell model
        K = K  + K_m(i)*exp(-t/lambda_m(i));
    end
    
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
        
%         %%%% TEST
%         
%         convint = 500*(xi(3,i) - xi(4,i));
%         
%         %%%% END TEST
        
        
        
        xidot(:,i) = [xi(3,i); xi(4,i);
            (-convint - k1*(xi(1,i) - xi(2,i)))/m1;% - g;
            (convint + k1*(xi(1,i) - xi(2,i)) ...
            - k2*(xi(2,i) - road_z(i)))/m2 - g];
        %xi(:,i+1) = xi(:,i) + dt(i)*xidot(:,i);
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
    
    % Debug purpose
    simplot = figure();
    plot(t,road_z,'k-'); hold on; plot(t,xi(1,:),'r-'); plot(t,xi(2,:),'b-');
    legend('road','sprung mass','unsprung mass'); hold off;
    
    % Objective functions
    f = [   max(xi(1,:)) - min(xi(1,:));                        % Comfort
            max(xi(2,:) - road_z) - min(xi(2,:) - road_z)   ];  % Handling
    
    
end
