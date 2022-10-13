clear all
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'defaultTextInterpreter','latex')
set(groot, 'defaultAxesFontsize',15)
set(groot, 'defaultTextFontsize',20)
set(groot, 'defaultLegendFontSize',15)

x=[-0.5,0;0,2;0.5,1];
x_m = mean(x(1:end-1,:),1);
rho = 1;
chi = 2;
gamma = 0.5;
sigma = 0.5;
x_r = x_m +rho*(x_m-x(3,:));
x_e = x_m+chi*(x_r-x_m);
x_c = x_m+gamma*(x_r-x_m);
x_cc = x_m-gamma*(x_m-x(3,:));
v = x(1,:)+sigma*(x-x(1,:));

fig=figure(1);
fig.Position(3:4)=[1000,400];
hold on
text(x(3,1),x(3,2),"$x_3$",'Interpreter',"latex","HorizontalAlignment","left",VerticalAlignment="middle")
x(4,:)=x(1,:);
plot(x(:,1),x(:,2),'k-')
text(x(1,1),x(1,2),"$x_1$",'Interpreter',"latex","HorizontalAlignment","center",VerticalAlignment="top")
text(x(2,1),x(2,2),"$x_2$",'Interpreter',"latex","HorizontalAlignment","center",VerticalAlignment="bottom")
x(3,:)=x_r;
plot(x(:,1),x(:,2),'k--')
text(x(3,1),x(3,2),"$r$",'Interpreter',"latex","HorizontalAlignment","right",VerticalAlignment="bottom")
x(3,:)=x_e;
text(x(3,1),x(3,2),"$e$",'Interpreter',"latex","HorizontalAlignment","right",VerticalAlignment="middle")
plot(x(:,1),x(:,2),'k--')
x(3,:)=x_c;
text(x(3,1),x(3,2),"$c$",'Interpreter',"latex","HorizontalAlignment","right",VerticalAlignment="bottom")
plot(x(:,1),x(:,2),'k--')
x(3,:)=x_cc;
text(x(3,1),x(3,2),"$cc$",'Interpreter',"latex","HorizontalAlignment","left",VerticalAlignment="cap")
v(4,:)=v(1,:);
plot(x(:,1),x(:,2),'k--')
plot(v(:,1),v(:,2),'r--')
text(v(3,1),v(3,2),"$v_3$",'Interpreter',"latex","HorizontalAlignment","left",VerticalAlignment="top")
text(v(2,1),v(2,2),"$v_2$",'Interpreter',"latex","HorizontalAlignment","right",VerticalAlignment="bottom")
ylim([-0.5,2.5])
xlim([-2,0.75])
ax=gca;
ax.Box='on';
