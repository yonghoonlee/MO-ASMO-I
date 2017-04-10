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
% Usage: [c,ceq] = NONLCON(x)
% Input: x
% Output: c,ceq
%   x: Points in design space
%   c: Inequality constraint function values of given design points
%   ceq: Equality constraint function values of given design points

% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.

function [c,ceq] = nonlcon(x,p)
    c = [];
    ceq = [];
%     for li = 1:size(x,1)
%         xlinear = 10.^x(li,:);
%         
%         % Parameters
%         m1 = p.m1;                              % 1/4 sprung mass [kg]
%         m2 = p.m2;                              % 1/4 unsprung mass [kg]
%         k2 = p.k2;                              % Tire stiffness [N/m]
%         g = p.g;                                % Gravitational acc [m/s^2]
%         v = p.v;                                % Vehicle velocity=60mph [m/s]
%         road_x = p.road_x;                      % Road profile in x [m]
%         road_z = p.road_z;                      % Road profile in z [m]
% 
%         % Time span
%         t = [road_x/v];                         % t [s]
%         dt = diff(t);                           % dt = t_{k+1} - t_{k} [s]
%         dt = [dt, dt(end)];
% 
%         % Decompose design variables
%         n_m = floor(length(xlinear)/2);         % Number of modes
%         K_m = xlinear(1:n_m);                   % K_m
%         lambda_m = xlinear((n_m + 1):(2*n_m));  % lambda_m
%         k1 = xlinear(end);                      % Suspension stiffness [N/m]
%         K = zeros(1,size(t,2));                 % Initialize kernel
%         for i = 1:n_m                           % Multimode Maxwell model
%             K = K  + K_m(i)*exp(-t/lambda_m(i));
%         end
% 
%         % Simulation
%         xi = zeros(4,size(t,2));                % Create state
%         xidot = zeros(4,size(t,2));             % Create derivative of state
% 
%         % Euler predictor-corrector method for simulation
%         for i = 1:(size(t,2)-1)
%             convint = 0;
%             for j = 1:i
%                 convint = convint ...
%                     + K(j)*(xi(3,(i-j+1)) - xi(4,(i-j+1)))*dt(j);
%             end
%             xidot(:,i) = [xi(3,i); xi(4,i);
%                 (-convint - k1*(xi(1,i) - xi(2,i)))/m1 - g;
%                 (convint + k1*(xi(1,i) - xi(2,i)) ...
%                 - k2*(xi(2,i) - road_z(i)))/m2 - g];
%             %xi(:,i+1) = xi(:,i) + dt(i)*xidot(:,i);
%             xist = xi(:,i) + dt(i)*xidot(:,i); % Predictor step (\xi_{*})
%             xitp = [xi(:,1:i), xist];
%             convintst = 0;
%             for j = 1:i+1
%                 convintst = convintst ...
%                     + K(j)*(xitp(3,(i-j+2)) - xitp(4,(i-j+2)))*dt(j);
%             end
%             fst = [xist(3); xist(4);
%                 (-convintst - k1*(xist(1) - xist(2)))/m1 - g;
%                 (convintst + k1*(xist(1) - xist(2)) ...
%                 - k2*(xist(2) - road_z(i+1)))/m2 - g];
%             xi(:,i+1) = xi(:,i) + dt(i)/2*(xidot(:,i) + fst);
%         end
% 
%         % Debug purpose
%         close all;
%         plot(t,road_z,'k-'); hold on; plot(t,xi(1,:),'r-'); plot(t,xi(2,:),'b-');
%         legend('road','sprung mass','unsprung mass'); hold off;
% 
%         % Objective functions
%         f1 = max(xi(1,:)) - min(xi(1,:)); % Comfort
%         f2 = max(xi(2,:) - road_z) - min(xi(2,:) - road_z); % Handling
%         
%         % Inequality
%         c = [c; max(f1 - 0.25, 0), max(f2 - 0.1, 0)];
%     end
%     
%     ceq = [];
    
    
end
