p.Omega = 10;
p.Ntex = 10;
p.ri = 0.5e-3;
p.ro = 20e-3;
p.Nr = 10;
p.Nth = 10;
p.Nz = 5;
p.eta = 9.624e-3;
p.rho = 873.4;
p.hmin =  269e-6;
p.hmax = 1000e-6;

xlb_geom = p.hmin*ones(((p.Nr+1)*p.Nth),1);
xub_geom = p.hmax*ones(((p.Nr+1)*p.Nth),1);
xtyp_geom = xlb_geom + rand(size(xlb_geom)).*xub_geom;
xtyp_fld = [0.013250667;
            0.006625333;
            0.001654582;
            1.42e-4;
            0.05;
            0.05];
xtyp = [xtyp_geom; xtyp_fld];

close all;
p.plot = true;

tic;
ObjectiveFunction = obj(xtyp,p);
subplot(3,1,1)
ylabel('M');
subplot(3,1,2)
ylabel('Fn');
subplot(3,1,3)
axis([0 160 1e-6 1e2])
set(gca,'YTick',[1e-6 1e-4 1e-2 1e0 1e2])
ylabel('residual');
xlabel('iteration');
M = ObjectiveFunction(1)
Fn = -ObjectiveFunction(2)
toc