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
% Usage: problem = SETTINGS(problem)
% Input: problem
% Output: problem
%   problem: Problem definition structure

% Multiobjective Adaptive Surrogate Modeling-based Optimization Toolbox I
% Author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Please refer to LICENSE.TXT for licensing details.
% Some directories may include codes from different author or with
% different license. In this case, please refer to LICENSE file or
% LICENSE.TXT file in each corresponding subdirectories.

function problem = settings(problem)

    % Case 2-1 Log normal continuous spectra
    problem.control.maxiter = 30*4;
    problem.sampling.initnumber = 10*4;
    problem.sampling.valnumber = 12 + 4;
    problem.sampling.upnumber = 2*4;
    problem.sampling.upexpnumber = 2*4;
    
    

    % =======================
    % Setting
    % =======================
    
    problem.highfidelity.expensive = 1;     % Not expensive
    problem.highfidelity.vectorized = 0;    % Function evaluation in scalar
    problem.p.m1 = 325;                     % 1/4 sprung mass [kg]
    problem.p.m2 = 65;                      % 1/4 unsprung mass [kg]
    problem.p.k2 = 232500;                  % Tire stiffness [N/m]
    problem.p.g = 9.806;                    % Gravitational acc [m/s^2]
    problem.p.v = [20,30,35,40,45,55,60,70,80]*1600/3600; % Vehicle velocity [m/s]
    problem.p.xlb = problem.xlb;
    problem.p.xub = problem.xub;
    
    load(fullfile(problem.probpath,'IRI_737b.mat'));
    problem.p.road_x = road_x;              % Road profile in x [m]
    problem.p.road_z = road_z;              % Road profile in z [m]
    
    problem.plotrange.xmin = 0;
    problem.plotrange.xmax = 4;
    problem.plotrange.ymin = 1.1;
    problem.plotrange.ymax = 1.5;
    
%    problem.plotcustom{1} = 'plotscript1.m';
%    problem.plotcustom{2} = 'plotscript2.m';
%    problem.plotcustom{3} = 'plotscript3.m';
    
end