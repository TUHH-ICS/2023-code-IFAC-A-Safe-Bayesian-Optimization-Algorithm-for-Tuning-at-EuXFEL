set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'defaultTextInterpreter','latex')
set(groot, 'defaultAxesFontsize',15)
set(groot, 'defaultTextFontsize',20)
set(groot, 'defaultLegendFontSize',15)

sig_n = 1;
X1=-5:0.000001:5;
r=[1;1];
o=[0;0];

fig=figure(1);

[x,y]=ellipse(X1,r,o);
plot([x,flip(x,2)],[y(1,:),flip(y(2,:),2)],'k--')
hold on
[x,y]=ellipse(X1,2*r,o);
plot([x,flip(x,2)],[y(1,:),flip(y(2,:),2)],'k')
plot(o(1),o(2),'rx')
grid on
xlabel("offset $w_1$",Interpreter="latex")
ylabel("slope $w_2$",Interpreter="latex")
legend("$\sigma$","$2 \sigma$", "$\mu$")
xlim([-2.25,2.25])
ylim([-2.25,2.25])


f=@(x) -1*x.^2+0.315;
x=[-4,-1.034,0.53,1.75,4,-2,2,3,-0.25];
Y=f(x)'+randn(9,1)*sig_n;
X=[ones(1,9);x];

fig2=figure(2);
fig2.Position = fig.Position;
%f=@(x) -1*x.^2+0.315;
plot(X1,f(X1),'k--')
hold on
plot(x,Y,'bx')
grid on
xlabel("$\mathcal{X}$")
ylabel("f")
legend("f","y")

A = sig_n^-2*(X*X')+eye(2);
o = sig_n^-2*inv(A)*X*Y;
r = diag(inv(A));

figure(3)
[x,y]=ellipse(X1,sqrt(r),o);
plot([x,flip(x,2)],[y(1,:),flip(y(2,:),2)],'k--')
hold on
[x,y]=ellipse(X1,2*sqrt(r),o);
plot([x,flip(x,2)],[y(1,:),flip(y(2,:),2)],'k')
plot(o(1),o(2),'rx')
grid on
xlabel("offset $w_1$",Interpreter="latex")
ylabel("slope $w_2$",Interpreter="latex")
legend("$\sigma$","$2 \sigma$", "$\mu$")
xlim([-2.25,2.25])
ylim([-2.25,2.25])

figure(4)
x_s = linspace(-5,5,1000);
X_s = [ones(1,1000);x_s];
mu = X_s'*o;
sig = sqrt(diag(X_s'*inv(A)*X_s));
plot(x_s,mu+2*sig,x_s,mu-2*sig)
hold on 
plot(X1,f(X1),'k--')
grid on
xlabel("$\mathcal{X}$",Interpreter="latex")
ylabel("f",Interpreter="latex")
legend("$\mu(x_*) + 2\sigma(x_*)$","$\mu(x_*) - 2\sigma(x_*)$","f")



function [x,y]=ellipse(X,r,o)
    y(1,:) = r(2)*sqrt(1-(X-o(1)).^2/(r(1)^2))+o(2);
    y(2,:) = -r(2)*sqrt(1-(X-o(1)).^2/(r(1)^2))+o(2);
    I = y(1,:) == real(y(1,:));
    y = y(:,I);
    x = X(I);
end