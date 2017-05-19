plotsc1t = [0,logspace(-3,2,101)];
plotsc1xlist = DATA{k,5};
plotsc1xlist = 10.^plotsc1xlist;
plotsc1nmode = round((size(plotsc1xlist,2)-1)/2);
plotsc1Klist = plotsc1xlist(:,1:plotsc1nmode);
plotsc1lambdalist = plotsc1xlist(:,(plotsc1nmode+1):(2*plotsc1nmode));
plotsc1k1list = plotsc1xlist(:,end);

plotsc1K = zeros(size(plotsc1k1list,1),size(plotsc1t,2));

for i = 1:plotsc1nmode
    plotsc1K = plotsc1K + repmat(plotsc1Klist(:,i),1,size(plotsc1t,2)).*exp(-repmat(plotsc1t,size(plotsc1k1list,1),1)./repmat(plotsc1lambdalist(:,i),1,size(plotsc1t,2)));
end

for i = 1:size(plotsc1K,1)
    if mod(i,10) == 0
        loglog(plotsc1t+0.01,plotsc1K(i,:),'LineWidth',1,...
            'Color',cmap(i,:));
        hold on;
    end
end

grid on;
axis([1e-2,1e2,1e0,3e3]);
xlabel('time, $t$','Interpreter','latex','FontSize',18);
ylabel('relaxation kernel, $K(t)$','Interpreter','latex','FontSize',18);
title([num2str(plotsc1nmode),'-Mode Case, Variable Speeds (20-80mph)'],'Interpreter','latex','FontSize',18);
set(gca,'FontSize',18);
set(gca,'TickLabelInterpreter','latex');
set(gca,'YTick',[1e0,1e1,1e2,2e2,3e2,4e2,6e2,1e3,2e3,3e3]);