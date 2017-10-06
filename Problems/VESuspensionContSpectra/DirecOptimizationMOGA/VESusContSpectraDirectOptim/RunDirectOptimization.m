% Direct optimization using MOGA

clear; close all; clc;

[A,b,Aeq,beq] = setup_lincon();
[xlb,xub,flb,fub] = setup_bounds();
problem.probpath = pwd;
problem.A = A;
problem.b = b;
problem.Aeq = Aeq;
problem.beq = beq;
problem.xlb = xlb;
problem.xub = xub;
problem.flb = flb;
problem.fub = fub;
problem = settings(problem);
initpop = [];

opt = gaoptimset('gamultiobj');
opt.PopulationSize = 200;
opt.ParetoFraction = 0.5;
opt.Generations = 200;
opt.StallGenLimit = 10;
opt.TolFun = 1e-5;
opt.TolCon = 1e-3;
opt.InitialPopulation = initpop;
opt.Display = 'iter';
opt.PlotFcns = @gaplotpareto;
opt.OutputFcns = @RunDO_Output;
opt.Vectorized = 'off';
opt.UseParallel = true;

[xopt,fopt,eopt,oopt,popt,sopt] = gamultiobj(@(x)obj(x,problem.p),...
    numel(problem.xlb),problem.A,problem.b,problem.Aeq,problem.beq,...
    problem.xlb,problem.xub,@(x)nonlcon(x,problem.p),opt);

save('RunDO_LogNormal');

scores = [];
for i = 1:size(scrhistory,3)
    scores = [scores; scrhistory(:,:,i)];
end

hf1 = figure();
plot(scores(:,1),scores(:,2),'r.'); hold on;
xlabel('f1: comfort metric');
ylabel('f2: handling metric');
title('objective function space');
legend('Log-normal model');




