clear; close all;

% Frequency range
omega = logspace(-1,5,1000);
n = 3; % number of modes

% Fixed parameters
m = 1; % mass
k = 1; % spring constant
Y = 1; % magnitude of base amplitude
zeta = 0;

% Calculated parameters
omega1 = sqrt(k/m); % natural frequency of mass m1 and spring k1
r = omega./omega1;

% Setting the fixed parameter vector
p.r = r;
p.m = m;
p.k = k;
p.omega1 = omega1;
p.omega = omega;
p.Y = Y;
p.n = n;
p.zeta = zeta;

% Maxwell (constant zeta, logeta)
% Design Variables (initial conditions):
%for i = 1:n
%    logeta(i) = -0.5;
%    loglambda(i) = -0.5;
%end
%x_mm = zeros(1,2*n);
%for i = 1:n
%    x_mm(1,2*i-1) = logeta(i);
%    x_mm(1,2*i) = loglambda(i);
%end
%x_mm = [-0.5 -0.5 -0.5 -0.5 -0.5 -0.5]; [omega_shifted_mm,X_mm,A_mm,maxAmp_mm,maxAcc_mm,AUC_Amp_mm,AUC_Acc_mm,halfAmp_mm]= dvi_maxwell_multi(x_mm,p); disp([maxAmp_mm,maxAcc_mm,AUC_Amp_mm,AUC_Acc_mm]);
%x_mm = [-0.5 -0.6 -0.7 -0.8 -0.9 -1.0]; [omega_shifted_mm,X_mm,A_mm,maxAmp_mm,maxAcc_mm,AUC_Amp_mm,AUC_Acc_mm,halfAmp_mm]= dvi_maxwell_multi(x_mm,p); disp([maxAmp_mm,maxAcc_mm,AUC_Amp_mm,AUC_Acc_mm]);
%x_mm = [-0.5 -0.5 -0.5 -0.5 -0.5 -0.5; -0.5 -0.6 -0.7 -0.8 -0.9 -1.0]; z = multimode_maxwell_multiobjective_obj(x_mm,p); disp(z);

for i = 1:n
    lb_x_mm(1,2*i-1) = -100; % log10(eta) lower bound
    lb_x_mm(1,2*i) = -100; % log10(lambda) lower bound
    ub_x_mm(1,2*i-1) = 5; % log10(eta) upper bound
    ub_x_mm(1,2*i) = 5; % log10(lambda) upper bound
end

opt = optimoptions('gamultiobj');
opt.CrossoverFraction = 0.5;
opt.ConstraintTolerance = 1e-6;
opt.FunctionTolerance = 1e-6;
opt.PopulationSize = 1000;
opt.MaxGenerations = 100;
opt.ParetoFraction = 0.1;
opt.UseVectorized = true;
opt.UseParallel = false;
opt.Display = 'iter';
opt.PlotFcn = @gaplotpareto;

[xopt,fopt,eflag,output,pop,score] = gamultiobj(...
    @(x) multimode_maxwell_multiobjective_obj(x,p), ...
    length(lb_x_mm),[],[],[],[],lb_x_mm,ub_x_mm,[],opt);

save('multimode_maxwell_multiobjective_result');

m1 = min(fopt(:,3));
m2 = max(fopt(:,3));
cgray = 1 - (fopt(:,3)-m1)/(m2-m1);
crgb = zeros(length(fopt(:,1)),3);
for i = 1:length(fopt(:,1))
    crgb(i,:) = interp1(linspace(0,1,64)',viridis(64),cgray(i));
end

%=============================
set(0,'DefaultTextInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter', 'latex');
fg1 = figure('Color',[1 1 1]);
%
subplot(4,6,1);
for i = 1:length(fopt(:,1))
    plot(xopt(i,1),fopt(i,1),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12); hold on;
end
set(gca,'FontSize',16); axis([-100 5 0 30]);
subplot(4,6,2);
for i = 1:length(fopt(:,1))
    plot(xopt(i,2),fopt(i,1),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12); hold on;
end
set(gca,'FontSize',16); axis([-100 5 0 30]);
subplot(4,6,3);
for i = 1:length(fopt(:,1))
    plot(xopt(i,3),fopt(i,1),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12); hold on;
end
set(gca,'FontSize',16); axis([-100 5 0 30]);
subplot(4,6,4);
for i = 1:length(fopt(:,1))
    plot(xopt(i,4),fopt(i,1),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12); hold on;
end
set(gca,'FontSize',16); axis([-100 5 0 30]);
subplot(4,6,5);
for i = 1:length(fopt(:,1))
    plot(xopt(i,5),fopt(i,1),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12); hold on;
end
set(gca,'FontSize',16); axis([-100 5 0 30]);
subplot(4,6,6);
for i = 1:length(fopt(:,1))
    plot(xopt(i,6),fopt(i,1),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12); hold on;
end
set(gca,'FontSize',16); axis([-100 5 0 30]);
%
subplot(4,6,7);
for i = 1:length(fopt(:,1))
    plot(xopt(i,1),fopt(i,2),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12); hold on;
end
set(gca,'FontSize',16); axis([-100 5 0 8e9]);
subplot(4,6,8);
for i = 1:length(fopt(:,1))
    plot(xopt(i,2),fopt(i,2),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12); hold on;
end
set(gca,'FontSize',16); axis([-100 5 0 8e9]);
subplot(4,6,9);
for i = 1:length(fopt(:,1))
    plot(xopt(i,3),fopt(i,2),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12); hold on;
end
set(gca,'FontSize',16); axis([-100 5 0 8e9]);
subplot(4,6,10);
for i = 1:length(fopt(:,1))
    plot(xopt(i,4),fopt(i,2),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12); hold on;
end
set(gca,'FontSize',16); axis([-100 5 0 8e9]);
subplot(4,6,11);
for i = 1:length(fopt(:,1))
    plot(xopt(i,5),fopt(i,2),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12); hold on;
end
set(gca,'FontSize',16); axis([-100 5 0 8e9]);
subplot(4,6,12);
for i = 1:length(fopt(:,1))
    plot(xopt(i,6),fopt(i,2),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12); hold on;
end
set(gca,'FontSize',16); axis([-100 5 0 8e9]);
%
subplot(4,6,13);
for i = 1:length(fopt(:,1))
    plot(xopt(i,1),fopt(i,3),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12); hold on;
end
set(gca,'FontSize',16); axis([-100 5 200 1000]);
subplot(4,6,14);
for i = 1:length(fopt(:,1))
    plot(xopt(i,2),fopt(i,3),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12); hold on;
end
set(gca,'FontSize',16); axis([-100 5 200 1000]);
subplot(4,6,15);
for i = 1:length(fopt(:,1))
    plot(xopt(i,3),fopt(i,3),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12); hold on;
end
set(gca,'FontSize',16); axis([-100 5 200 1000]);
subplot(4,6,16);
for i = 1:length(fopt(:,1))
    plot(xopt(i,4),fopt(i,3),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12); hold on;
end
set(gca,'FontSize',16); axis([-100 5 200 1000]);
subplot(4,6,17);
for i = 1:length(fopt(:,1))
    plot(xopt(i,5),fopt(i,3),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12); hold on;
end
set(gca,'FontSize',16); axis([-100 5 200 1000]);
subplot(4,6,18);
for i = 1:length(fopt(:,1))
    plot(xopt(i,6),fopt(i,3),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12); hold on;
end
set(gca,'FontSize',16); axis([-100 5 200 1000]);
%
subplot(4,6,19);
for i = 1:length(fopt(:,1))
    plot(xopt(i,1),fopt(i,4),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12); hold on;
end
set(gca,'FontSize',16); axis([-100 5 0 3e11]);
subplot(4,6,20);
for i = 1:length(fopt(:,1))
    plot(xopt(i,2),fopt(i,4),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12); hold on;
end
set(gca,'FontSize',16); axis([-100 5 0 3e11]);
subplot(4,6,21);
for i = 1:length(fopt(:,1))
    plot(xopt(i,3),fopt(i,4),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12); hold on;
end
set(gca,'FontSize',16); axis([-100 5 0 3e11]);
subplot(4,6,22);
for i = 1:length(fopt(:,1))
    plot(xopt(i,4),fopt(i,4),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12); hold on;
end
set(gca,'FontSize',16); axis([-100 5 0 3e11]);
subplot(4,6,23);
for i = 1:length(fopt(:,1))
    plot(xopt(i,5),fopt(i,4),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12); hold on;
end
set(gca,'FontSize',16); axis([-100 5 0 3e11]);
subplot(4,6,24);
for i = 1:length(fopt(:,1))
    plot(xopt(i,6),fopt(i,4),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12); hold on;
end
set(gca,'FontSize',16); axis([-100 5 0 3e11]);
%
ax1 = axes('Position',[0 0 1 1],'Visible','off');
text(0.175,0.05,'$x_1$','FontSize',18);
text(0.310,0.05,'$x_2$','FontSize',18);
text(0.445,0.05,'$x_3$','FontSize',18);
text(0.580,0.05,'$x_4$','FontSize',18);
text(0.715,0.05,'$x_5$','FontSize',18);
text(0.850,0.05,'$x_6$','FontSize',18);
text(0.085,0.190,'$f_4$','FontSize',18);
text(0.085,0.405,'$f_3$','FontSize',18);
text(0.085,0.620,'$f_2$','FontSize',18);
text(0.085,0.835,'$f_1$','FontSize',18);
set(gcf,'Position',[10,200,1280,720]);
%
eval(['export_fig ', 'multimode_maxwell_multiobjective_moga_plot_x_f', ' -png -r200']);

fg2 = figure('Color',[1 1 1]);
subplot(3,3,1);
for i = 1:length(fopt(:,1))
    plot(fopt(i,1),fopt(i,2),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12);
    hold on;
end
set(gca,'FontSize',16);
subplot(3,3,4);
for i = 1:length(fopt(:,1))
    plot(fopt(i,1),fopt(i,3),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12);
    hold on;
end
set(gca,'FontSize',16);
subplot(3,3,7);
for i = 1:length(fopt(:,1))
    plot(fopt(i,1),fopt(i,4),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12);
    hold on;
end
set(gca,'FontSize',16);
subplot(3,3,5);
for i = 1:length(fopt(:,1))
    plot(fopt(i,2),fopt(i,3),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12);
    hold on;
end
set(gca,'FontSize',16);
subplot(3,3,8);
for i = 1:length(fopt(:,1))
    plot(fopt(i,2),fopt(i,4),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12);
    hold on;
end
set(gca,'FontSize',16);
subplot(3,3,9);
for i = 1:length(fopt(:,1))
    plot(fopt(i,3),fopt(i,4),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',crgb(i,:),'MarkerSize',12);
    hold on;
end
set(gca,'FontSize',16);
%
ax1 = axes('Position',[0 0 1 1],'Visible','off');
text(0.225,0.05,'$f_1$','FontSize',18);
text(0.500,0.05,'$f_2$','FontSize',18);
text(0.775,0.05,'$f_3$','FontSize',18);
text(0.05,0.225,'$f_4$','FontSize',18);
text(0.05,0.525,'$f_3$','FontSize',18);
text(0.05,0.825,'$f_2$','FontSize',18);
set(gcf,'Position',[10,200,720,720]);
%
eval(['export_fig ', 'multimode_maxwell_multiobjective_moga_plot_f_f', ' -png -r200']);

fg3 = figure('Color',[1 1 1]);
foptmin = min(fopt,[],1);
foptmax = max(fopt,[],1);
foptnr = (fopt - repmat(foptmin,size(fopt,1),1))./repmat((foptmax-foptmin),size(fopt,1),1);
[~,idx] = sort(foptnr(:,1));
foptn = foptnr(idx,:);
crgbn = crgb(idx,:);
for i = 1:size(foptn,1)
    plot([1 2 3 4], foptn(i,:), '-','Color',crgbn(i,:),'LineWidth',2); hold on;
    %drawnow;
    %pause(0.1);
end
xlabel('objectives'); ylabel('normalized objective function value');
axis([0.8 4.2 -0.1 1.1]);
set(gca,'FontSize',16);
set(gca,'XTick',[1 2 3 4]);
set(gca,'XGrid','on');
%
eval(['export_fig ', 'multimode_maxwell_multiobjective_moga_plot_obj1-4', ' -png -r200']);
