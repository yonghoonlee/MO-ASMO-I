function problem = settings(problem)

    problem.control.maxiter = 100;
    problem.sampling.initnumber = 40;
    problem.sampling.valnumber = 20;
    problem.sampling.upnumber = 30;
    problem.sampling.upexpnumber = 10;
    problem.gamultiobj.opt.Generations = 30;
    problem.gamultiobj.opt.PopulationSize = 3000;
    
    problem.highfidelity.expensive = 1;     % Expensive
    problem.highfidelity.vectorized = 0;    % Function evaluation in scalar

    % Operating condition
    problem.p.Omega = 10;

    % Geometry
    problem.p.Ntex = 10;
    problem.p.ri = 0.5e-3;
    problem.p.ro = 20e-3;

    % Mesh
    problem.p.Nr = 5;
    problem.p.Nth = 5;
    problem.p.Nz = 3;

    % Fluid properties
    problem.p.eta = 9.624e-3;
    problem.p.rho = 873.4;
    %problem.p.etap1 = 0.013250667;
    %problem.p.etap2 = 0.006625333;
    %problem.p.lambda1 = 0.001654582;
    %problem.p.lambda2 = 1.42e-4;
    %problem.p.alpha1 = 0.05;
    %problem.p.alpha2 = 0.05;

    %problem.plotrange.xmin = 0;
    %problem.plotrange.xmax = 4;
    %problem.plotrange.ymin = 1.1;
    %problem.plotrange.ymax = 1.5;
    
    %problem.plotcustom{1} = 'plotscript1.m';
    %problem.plotcustom{2} = 'plotscript2.m';
    %problem.plotcustom{3} = 'plotscript3.m';
    
end