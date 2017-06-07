function [x, p] = Giesekus_2M3D_v002_default_parameters()
    
    % Operating condition
    p.Omega = 10;
    
    % Geometry
    p.Ntex = 10;
    p.ri = 0.5e-3;
    p.ro = 20e-3;
    
    % Mesh
    p.Nr = 10;
    p.Nth = 10;
    p.Nz = 5;
    
    % Fluid properties
    p.eta = 9.624e-3;
    p.rho = 873.4;
    p.etap1 = 0.013250667;
    p.etap2 = 0.006625333;
    p.lambda1 = 0.001654582;
    p.lambda2 = 1.42e-4;
    p.alpha1 = 0.05;
    p.alpha2 = 0.05;
    
    xH = 269e-6 * ones(p.Nr+1,p.Nth);
    
    x = [reshape(xH,numel(xH),1); p.etap1; p.etap2; p.lambda1; p.lambda2; p.alpha1; p.alpha2];
    
end