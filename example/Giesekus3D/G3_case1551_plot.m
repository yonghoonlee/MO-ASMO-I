%===============================================================================
% Multiobjective Adaptive Surrogate Modeling-based Optimization Code I
% Main author: Yong Hoon Lee (ylee196@illinois.edu, yonghoonlee@outlook.com)
% Link: https://github.com/yonghoonlee/MO-ASMO-I
%===============================================================================
% Texture and fluid design problem using 3D Cauchy momentum equation solver with
% 2-mode Giesekus model.
% Plot Pareto set
%===============================================================================
function G3_case1551_plot
    close all;
    plotPreparation;
    [mpath,mname] = fileparts(mfilename('fullpath'));
    %fname = 'G3_case1551_final.mat';
    fname = 'G3_case1551_iter0020.mat';
    plotexport = true;
    Nr = 5;
    Nth = 5;
    hmin = 269e-6;
    hmax = 800e-6;
    Ntex = 10;
    ri = 0.5e-3;
    ro = 20e-3;
    xlb_fld = [0.00001; % min(etap1,etap2)
               0.00001; % min(etap1,etap2)
               0.00001; % min(lambda1,lambda2)
               0.00001; % min(lambda1,lambda2)
               0.01;    % min(alpha1,alpha2)
               0.01];   % min(alpha1,alpha2)
    xub_fld = [0.1;     % max(etap1,etap2)
               0.1 ;    % max(etap1,etap2)
               0.01;    % max(lambda1,lambda2)
               0.01;    % max(lambda1,lambda2)
               0.99;    % max(alpha1,alpha2)
               0.99];   % max(alpha1,alpha2)
    %-------------------------------------------------------------------------------
    plotPreparation;
    [mpath,~] = fileparts(mfilename('fullpath'));
    load(fullfile(mpath,'solution',fname),'result');
    cm = parula(64);
    %-------------------------------------------------------------------------------
    % Data preparation
    if ~isstruct(result)
        R.data = result;
        result = R;
    end
    result.population = cell2mat(table2array(result.data(end,7)));
    result.scores = cell2mat(table2array(result.data(end,8)));
    [xsort,fsort,isort] = ndSort(result.population, result.scores);
    xpareto = xsort(isort==1,:);
    fpareto = fsort(isort==1,:);
    mx = size(xpareto,2);
    pareto = [xpareto,fpareto];
    pareto = sortrows(pareto,mx+1);
    xpareto = pareto(:,1:mx);
    fpareto = pareto(:,(mx+1):end);
    result.x = xpareto;
    result.fval = fpareto;
    %-------------------------------------------------------------------------------
    % Mesh preparation
    [~,~,~,~,zth,~] = semhat(Nth);
    [~,~,~,~,zr,~] = semhat(Nr);
    phi = (2*pi/Ntex);
    r = (1+zr)/2*ro+(1-zr)/2*ri;
    theta = phi/2*zth;
    [Rmat,Theta] = ndgrid(r,theta);
    X2d = Rmat.*cos(Theta);
    Y2d = Rmat.*sin(Theta);
    Nplot = 100;
    [~,~,~,~,zplot,~] = semhat(Nplot);
    rplot = (1+zplot)/2*ro+(1-zplot)/2*ri;
    thetaplot = phi/2*zplot;
    [Rplot,Tplot] = ndgrid(rplot,thetaplot);
    X2dplot = Rplot.*cos(Tplot);
    Y2dplot = Rplot.*sin(Tplot);
    %---------------------------------------------------------------------------
    % PLOT1
    for idx = 1:size(fpareto,1)
        Hmat = reshape(xpareto(idx,1:((Nr+1)*Nth)), (Nr+1), Nth);
        Hmat = [Hmat, Hmat(:,1)];
        Fpareto = fpareto(idx,:);
        Hplot = lagrangeR2d(r,theta,Hmat,rplot,thetaplot);
        Hplot(Hplot<hmin) = hmin;
        Hplot(Hplot>hmax) = hmax;
        %-----------------------------------------------------------------------
        try % Open figure window
            figure(fg1);
        catch
            fg1 = figure('Color',[1 1 1]);
        end
        %-----------------------------------------------------------------------
        subplot(1,2,1);
        hold off;
        [foc, foh] = contourf(X2dplot*1e3, Y2dplot*1e3, -Hplot*1e3);
        foh.LineWidth = 0.5;
        foh.LevelListMode = 'manual';
        foh.LevelList = linspace(-hmax*1e3, -hmin*1e3, 12);
        view(-90,90);
        caxis([-hmax*1e3 -hmin*1e3]);
        ax1 = gca;
        ax1.DataAspectRatio = [1 1 1];
        ax1.Visible = 'off';
        colormap(cm);
        text(3, 6, ['$f_1$=$',num2str(Fpareto(1),'%7.4f'),'$'],...
            'FontSize',13);
        text(1.5, 6, ['$f_2$=$',num2str(-Fpareto(2),'%7.4f'),'$'],...
            'FontSize',13);
        %-----------------------------------------------------------------------
        subplot(1,2,2);
        hold off;
        xfld_normalized = (xpareto(idx,(end-5):end) ...
                - reshape(xlb_fld,1,numel(xlb_fld))) ...
            ./ (reshape(xub_fld,1,numel(xub_fld)) ...
                - reshape(xlb_fld,1,numel(xlb_fld)));
        xfld_normalized = fliplr(xfld_normalized);
        barh(1:6,xfld_normalized);
        ax2 = gca;
        ax2.Visible = 'on';
        ax2.Box = 'on';
        ax2.XLim = [0 1];
        ax2.XTick = [0 1];
        ax2.XTickLabelMode = 'manual';
        ax2.XTickLabel = {'LB';'UB'};
        ax2.YAxis.Color = [0 0 0];
        ax2.YAxis.TickLength = [0 0];
        ax2.YAxis.TickLabelsMode = 'manual';
        ax2.YAxis.TickLabels = {'$\alpha_2$'; '$\alpha_1$'; ...
            '$\lambda_2$'; '$\lambda_1$'; ...
            '$\eta_{p2}$'; '$\eta_{p1}$'};
        ax2.FontSize = 13;
        ax1.Position = [0.100 0.110 0.500 0.815];
        ax2.Position = [0.540 0.175 0.140 0.515];
        %-----------------------------------------------------------------------
        if plotexport
            eval(['export_fig ''', ...
                fullfile(mpath,'plot',...
                [mname,'_ParetoDesigns_id',num2str(idx,'%04d')]), ''' -pdf']);
        end
    end
    %---------------------------------------------------------------------------
    % PLOT2
    try % Open figure window
        figure(fg2);
    catch
        fg2 = figure('Color',[1 1 1]);
    end
    %---------------------------------------------------------------------------
    nsort = isort(end);
    cm = plasma(nsort+6);
    cm = cm(4:(end-3),:);
    xsc{1} = xsort(find(isort==1),:);
    fsc{1} = fsort(find(isort==1),:);
    ph1 = plot(fsc{1}(:,1),-fsc{1}(:,2),'.','Color',cm(1,:),'MarkerSize',10);
    hold on;
    for idx = 2:nsort
        xsc{idx} = xsort(find(isort==idx),:);
        fsc{idx} = fsort(find(isort==idx),:);
        plot(fsc{idx}(:,1),-fsc{idx}(:,2),'.',...
            'Color',cm(idx,:),'MarkerSize',10);
        hold on;
    end
    xfsc1 = [xsc{1}, fsc{1}];
    xfsc1 = sortrows(xfsc1,(size(xsort,2)+1));
    xsc1 = xfsc1(:,1:(size(xsort,2)));
    fsc1 = xfsc1(:,(size(xsort,2)+1):end);
    cm = gray(ceil(1.6*size(fsc1,1)));
    ph2 = plot(fsc1(1,1),-fsc1(1,2),'o','Color',[0 0 0],...
        'MarkerSize',5,'LineWidth',2);
    for idx = 2:size(fsc1,1)
        plot(fsc1(idx,1),-fsc1(idx,2),'o','Color',cm(idx,:),...
            'MarkerSize',5,'LineWidth',2);
    end
    ax = gca;
    ax.FontSize = 16;
    xlabel('$f_1$: normalized apparent viscosity');
    ylabel('$f_2$: normal force');
    legend([ph1, ph2], {'Explored designs', 'Non-dominated designs'},...
        'Location','northwest');
    axis([(0 - eps),(15 + eps),(-0.05 - eps),(0.3 + eps)]);
    %---------------------------------------------------------------------------
    if plotexport
        eval(['export_fig ''', ...
            fullfile(mpath,'plot',...
            [mname,'_ParetoDesigns_objective']), ''' -pdf']);
    end
end
%===============================================================================
function fo = lagrangeR2d(ri,ti,fi,ro,to)
    nri = length(ri);
    nro = length(ro);
    nti = length(ti);
    nto = length(to);
    frito = zeros(nri,nto);
    for idx = 1:nri
        frito(idx,:) = reshape(lagrange(reshape(ti,1,nti),fi(idx,:),reshape(to,1,nto)),1,nto);
    end
    fo = zeros(nro,nto);
    for idx = 1:nto
        fo(:,idx) = reshape(lagrange(reshape(ri,nri,1),frito(:,idx),reshape(ro,nro,1)),nro,1);
    end
end
%===============================================================================
function fo = lagrange(xi,fi,xo)
    n = length(xi);
    fo = zeros(size(xo));
    for idx = 1:n
        w = ones(size(xo));
        for jdx = [1:(idx-1), (idx+1):n]
            w = (xo - xi(jdx))./(xi(idx) - xi(jdx)).*w;
        end
        fo = fo + w*fi(idx);
    end
end
%===============================================================================
function [semhatAh,semhatBh,semhatCh,semhatDh,semhatz,semhatw] = semhat(semhatN)
    [semhatz,semhatw] = zwgll(semhatN); %  z = [z0,z1,...,zN] Gauss-Lobatto points 
                                        %  w = [w0,w1,...,wN] Gauss-Lobatto weights 
    semhatBh    = diag(semhatw);        %  1-D mass matrix
    semhatDh    = dhat(semhatz);        %  1-D derivative matrix
    semhatAh    = semhatDh'*semhatBh*semhatDh;  %  1-D stiffness matrix
    semhatCh    = semhatBh*semhatDh;            %  1-D convection operator
end
%===============================================================================
function [zwgllz,zwgllw] = zwgll(zwgllp)
    zwglln = zwgllp+1;
    zwgllz = zeros(1,zwglln);
    zwgllw = zeros(1,zwglln);
    zwgllz(1)=-1;
    zwgllz(zwglln)= 1;
    if zwgllp>1
        if zwgllp==2 
            zwgllz(2)=0; 
        else
            zwgllM=zeros(zwgllp-1,zwgllp-1);
            for zwglli=1:zwgllp-2
                zwgllM(zwglli,zwglli+1)=(1/2)*sqrt((zwglli*(zwglli+2))/((zwglli+1/2)*(zwglli+3/2)));
                zwgllM(zwglli+1,zwglli)=zwgllM(zwglli,zwglli+1);
            end
            zwgllz(2:zwgllp)=sort(real(eig(zwgllM)));
        end
    end
    %compute the weights w
    zwgllw(1)=2/(zwgllp*(zwglln));
    zwgllw(zwglln)=zwgllw(1);
    for zwglli=2:zwgllp

        zwgllx=zwgllz(zwglli);
        zwgllz0=1;
        zwgllz1=zwgllx;
        zwgllz2=0;
        for zwgllj=1:zwgllp-1
            zwgllz2=zwgllx.*zwgllz1*(2*zwgllj+1)/(zwgllj+1)-zwgllz0*zwgllj/(zwgllj+1);
            zwgllz0=zwgllz1;
            zwgllz1=zwgllz2;
        end
        zwgllw(zwglli)=2/(zwgllp*(zwglln)*zwgllz2*zwgllz2);
    end
    zwgllz=zwgllz';
    zwgllw=zwgllw';
end
%===============================================================================
function[fdwfc] = fd_weights_full(fdwfxx,fdwfx,fdwfm)
    fdwfn1 = length(fdwfx);
    fdwfn  = fdwfn1-1;
    fdwfm1 = fdwfm+1;
    fdwfc1       = 1.;
    fdwfc4       = fdwfx(1) - fdwfxx;
    fdwfc = zeros(fdwfn1,fdwfm1);
    fdwfc(1,1) = 1.;
    for fdwfi=1:fdwfn
        fdwfi1  = fdwfi+1;
        fdwfmn = min(fdwfi,fdwfm);
        fdwfc2 = 1.;
        fdwfc5 = fdwfc4;
        fdwfc4 = fdwfx(fdwfi1)-fdwfxx;
        for fdwfj=0:fdwfi-1;              fdwfj1  = fdwfj+1;
            fdwfc3 = fdwfx(fdwfi1)-fdwfx(fdwfj1);
            fdwfc2 = fdwfc2*fdwfc3;
            for fdwfk=fdwfmn:-1:1;         fdwfk1  = fdwfk+1;
                fdwfc(fdwfi1,fdwfk1) = fdwfc1*(fdwfk*fdwfc(fdwfi1-1,fdwfk1-1)-fdwfc5*fdwfc(fdwfi1-1,fdwfk1))/fdwfc2;
            end
            fdwfc(fdwfi1,1) = -fdwfc1*fdwfc5*fdwfc(fdwfi1-1,1)/fdwfc2;
            for fdwfk=fdwfmn:-1:1;         fdwfk1  = fdwfk+1;
                fdwfc(fdwfj1,fdwfk1) = (fdwfc4*fdwfc(fdwfj1,fdwfk1)-fdwfk*fdwfc(fdwfj1,fdwfk1-1))/fdwfc3;
            end
            fdwfc(fdwfj1,1) = fdwfc4*fdwfc(fdwfj1,1)/fdwfc3;
        end
        fdwfc1 = fdwfc2;
    end
end
%===============================================================================
function [dhatDh] = dhat(dhatx)
    dhatn1 = length(dhatx);
    dhatw  = zeros(dhatn1,2);
    dhatDh = zeros(dhatn1,dhatn1);
    for dhati=1:dhatn1
        dhatw = fd_weights_full(dhatx(dhati),dhatx,1);  % Bengt Fornberg's interpolation algorithm
        dhatDh(:,dhati) = dhatw(:,2);                   % Derivative of pn is stored in column 2.
    end
    dhatDh = dhatDh';
end
%===============================================================================

