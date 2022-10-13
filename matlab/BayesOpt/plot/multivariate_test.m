mvnpdf([0,1],[0,1],[1,0.99;0.99,1])
[X,Y] = meshgrid(-3:0.1:3);
x=[reshape(X,[],1),reshape(Y,[],1)];
%x=build_testP([-5,5;-5,5],200);
y=mvnpdf(x,[0,0],[1,0.6;0.6,1]);
Z=reshape(y,size(X));
surface(X,Y,Z)
colorbar
hold on
arg = find(X(1,:)==-1);
[zmax1,I1] = max(Z(:,arg));
[zmax2,I2] = max(Z(arg,:));
plot3(X(1:I1,arg),Y(1:I1,arg),ones(I1,1),'r--',X(I1,1:arg),Y(I1,1:arg),ones(arg,1),'r--','LineWidth',2)
title('probability density function of a bivariate Gaussian distribution','Interpreter','latex','FontSize',13)
xlabel('$$f(x_1)$$','Interpreter','latex')
ylabel('$$f(x_2)$$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex');

%plot3(x(:,1),x(:,2),y)

figure(2)
plot(Y(:,arg),Z(:,arg))

figure(3)
[X,Y] = meshgrid(-1:0.01:1);
X=gpuArray(X);
Y=gpuArray(Y);
x=[reshape(X,[],1),reshape(Y,[],1)];
%x=build_testP([-5,5;-5,5],200);
y=mvnpdf(x,[0,0],[1,0.6;0.6,1]);
Z=reshape(y,size(X));
arg = find(X(1,:)==-1);
[zmax1,I1] = max(Z(:,arg));
plot(gather(X(Z==zmax1)),gather(Y(Z==zmax1)),'.')

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