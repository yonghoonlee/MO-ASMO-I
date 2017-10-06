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

    % Decompose design variables
    Hmax = 10^x(1);
    taumax = 10^x(2);
    sigma = 10^x(3);
    Ge = x(4);
    k1 = 10^x(5);                           % Suspension stiffness [N/m]

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
        c = [NaN];
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
        c = [NaN];
        return;
    end

    Kmax = Ge + real(integral(@(s) H(exp(s)),log(xlbnd),log(xrbnd)));
    if abs(Kmax) > 1e16
    	c = [NaN];
    end
    c = Kmax - 10^(p.xub(end));

end
