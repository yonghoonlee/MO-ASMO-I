set(0,'DefaultAxesTickLabelInterpreter','latex');
set(0,'DefaultColorbarTickLabelInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultPolaraxesTickLabelInterpreter','latex');
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultTextarrowshapeInterpreter','latex');
set(0,'DefaultTextboxshapeInterpreter','latex');
set(0,'DefaultFigurePosition',[20 150 560 420]);

[mpath,afile] = fileparts(mfilename('fullpath'));
[~,mfile] = regexp(afile,'_analyze','match','split');
mfile = cell2mat(mfile(1));
load(fullfile(mpath,'solution',[mfile,'_final.mat']));

n = size(result.data,1);
cm = flipud(plasma(n));
fga1 = figure('Color',[1 1 1]);
for idx = n:-1:1
    PoolHffFFea = cell2mat(result.data.c08_PoolHffFFea(idx));
    plot(PoolHffFFea(:,1),PoolHffFFea(:,2),'.','MarkerSize',9,'Color',cm(idx,:));
    hold on;
end
caxis([0.5, n+0.5]);
colormap(cm);
cb = colorbar;
ticks = [1, n];
if n > 6
    ticks = sort(round([ticks, linspace(round(n/6),round(n/6*5),5)]));
end
cb.Ticks = ticks;
