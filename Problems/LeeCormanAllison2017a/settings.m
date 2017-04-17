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

    problem.control.maxiter = 40;
    problem.sampling.initnumber = 20;
    
    problem.highfidelity.expensive = 0;     % Not expensive
    problem.highfidelity.vectorized = 0;    % Function evaluation in scalar
    problem.p.m1 = 325;                     % 1/4 sprung mass [kg]
    problem.p.m2 = 65;                      % 1/4 unsprung mass [kg]
    problem.p.k2 = 232500;                  % Tire stiffness [N/m]
    problem.p.g = 9.806;                    % Gravitational acc [m/s^2]
    problem.p.v = 60*1600/3600;             % Vehicle velocity=60mph [m/s]
    
    load(fullfile(problem.probpath,'IRI_737b.mat'));
    problem.p.road_x = road_x;              % Road profile in x [m]
    problem.p.road_z = road_z;              % Road profile in z [m]
    
    problem.plotrange.xmin = 0.0045;
    problem.plotrange.xmax = 0.006;
    problem.plotrange.ymin = 0.15;
    problem.plotrange.ymax = 0.18;
    
    % Approximated anchor point solutions obtained from other method
    problem.approximated.anchor1_x ...
        = [ 1.601440240667433068111336069705430418252945e+00,...
            3.287396192937855676774461244349367916584015e+00,...
            5.148919924139860881950880866497755050659180e-01,...
            2.583896432208809823549700013245455920696259e-01,...
            2.500154854692612005351293191779404878616333e+00 ];
    problem.approximated.anchor1_f ...
        = [ 4.591321895652368881290961155627883272245526e-03,...
            1.749896094060851559071068095363443717360497e-01 ];
    problem.approximated.anchor2_x ...
        = [ 2.299994088071814601903497532475739717483521e+00,...
            3.299996580890698805177407848532311618328094e+00,...
           -9.999818143393159886045395978726446628570557e-01,...
           -9.999813468601501664778652411769144237041473e-01,...
            2.799992852481252292307090101530775427818298e+00 ];
    problem.approximated.anchor2_f ...
        = [ 5.768371395759553169058087007670110324397683e-03,...
            1.575889905018117154167356375182862393558025e-01 ];
end