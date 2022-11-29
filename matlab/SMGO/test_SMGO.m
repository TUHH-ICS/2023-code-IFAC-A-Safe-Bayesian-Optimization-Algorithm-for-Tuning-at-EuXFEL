clear all
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'defaultTextInterpreter','latex')
set(groot, 'defaultAxesFontsize',12)
set(groot, 'defaultTextFontsize',12)
set(groot, 'defaultLegendFontsize',12)

% obj = @(x) 2*sin(1*pi*x(:,1))-cos(0.21*x(:,1))+x(:,1).^2+1.4*x(:,2).^2;
% cond = [-4, 4;
%     -2.5, 1];
% x0 = [0.5,0];
% xs = build_testP(cond,100);
% figure(4)
% plot3(xs(:,1),xs(:,2),obj(xs))
% 
obj = @(x) 2*sin(1*pi*x(:,1))-cos(0.21*x(:,1))+x(:,1).^2;
cond = [-4, 4];
x0 = -2;
xs = build_testP(cond,10000);
fig=figure(4);
fig.Units='centimeters';
fig.Position(3:end)= [12.5,10.5];
plot(xs,obj(xs))


opts.L_f = 8;
opts.L_g = 8;
opts.start_points = 0;
opts.maxIt = 75;
opts.plot = 1;
opts.delta = 0.2;
opts.B = 10;
opts.safety.threshold = 5;
g = {@(x) opts.safety.threshold-obj(x)};

[X,Y,C,L_f,L_g] = SMGO(cond,obj,g,x0,opts);
[yopt,I] = min(Y);
X(I,:)
[ybad,I2]=max(Y);
% x = cond(1):0.1:cond(2);
% figure(2)
% plot(x,obj(x))
% hold on
[y_d,I_d]=min(obj(xs));
if size(cond,1) == 1
    figure(4)
    hold on
    plot(X(I,:),yopt,'gx',X(I2,:),ybad,'rx',xs(I_d),y_d,'bx')
    yline(opts.safety.threshold,'--','threshold',LineWidth=2)
    legend("","SMGO optimum","SMGO highest value","true optimum","")
    hold off
end


function Xs=build_testP(cond, samples)
    D = size(cond,1);
    xs=zeros(samples,D);
    for i=1:D
        xs(:,i)=linspace(cond(i,1),cond(i,2),samples);
    end
    n = size(xs,1);
    if D==2
        Xs = zeros(n^D,D);
        for i=1:n
            if mod(i,2)==0
                Xs((i-1)*n+1:i*n,2)=flip(xs(:,2));
            else
                Xs((i-1)*n+1:i*n,2)=xs(:,2);
            end
            Xs((i-1)*n+1:i*n,1)=ones(n,1)*xs(i,1); 
        end
    else 
        Xs=xs;
    end
end