hyp = categorical({'1/2', '1/3', '1/5', '1/7', '1/9'});
nopt = [16, 180, 321, 572, 144];
yopt = [10.6864, 10.3875, 10.3298, 10.33, 10.3816];
nmax = [19, 506, 1005, 911, 1231];
f=figure();

subplot(3,1,1)

p1 = bar(hyp,nmax);
xtips1 = p1.XEndPoints;
ytips1 = p1.YEndPoints;
labels1 = string(p1.YData);
text(xtips1(1),ytips1(1),labels1(1),'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
grid on
xlabel("Scale")
ylabel("$$\mbox{n}_{max}$$",'Interpreter','latex', 'FontSize', 12)
ax=gca;
ax.GridLineStyle = ':';
ax.TickLength = ax.TickLength / 2;
title(ax,'Naive BO Cov. Hyperparameter Dependency')


subplot(3,1,2)
p2 = bar(hyp,nopt);
xtips1 = p2.XEndPoints;
ytips1 = p2.YEndPoints;
labels1 = string(p2.YData);
text(xtips1(1),ytips1(1),labels1(1),'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
grid on
xlabel("Scale")
ylabel("$$\mbox{n}_{opt}$$",'Interpreter','latex', 'FontSize', 12)
ax=gca;
ax.GridLineStyle = ':';
ax.TickLength = ax.TickLength / 2;


subplot(3,1,3)
p3 = plot(hyp,yopt);
xlabel("Scale")
ylabel("$$\mbox{J}_{opt} \mbox{ [fs]}$$",'Interpreter','latex', 'FontSize', 12)

grid on



%p.XTickLabel = hyp;

