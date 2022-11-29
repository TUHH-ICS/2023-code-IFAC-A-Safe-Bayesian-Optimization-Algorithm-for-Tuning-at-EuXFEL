set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'defaultTextInterpreter','latex')
set(groot, 'defaultAxesFontsize',15)
set(groot, 'defaultTextFontsize',20)
set(groot, 'defaultLegendFontSize',14)
x=linspace(0,30,1000);
A=1;
f=0.1;
var = 0.3;
n = @(x) var*sin(2*pi*f/2*x); 
s=@(x) A*sin(2*pi*f*x);
s_n = @(x) A*sin(2*pi*f*x+n(x));

plot(x,s(x),x,s_n(x),x,n(x),'--')
grid on
legend("$s(t)$","$s_n(t)$","$\phi_n(t)$",'FontSize',14,'Location','Northeast')
ylim([-1.5 3])
ax=gca;
ZoomPlot(ax)