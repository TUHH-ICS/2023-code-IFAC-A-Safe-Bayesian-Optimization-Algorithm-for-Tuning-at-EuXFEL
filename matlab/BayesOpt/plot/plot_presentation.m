clear all
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'defaultTextInterpreter','latex')
set(groot, 'defaultAxesFontsize',15)
set(groot, 'defaultTextFontsize',20)
set(groot, 'defaultLegendFontSize',15)
% Kp_max = 3e1;
% Kp_min = 0.2;
% Ki_max = 6e1;
% Ki_min = 0;
% cond=[Kp_min, Kp_max];
%       %Ki_min, Ki_max];
% 
% cov_matern = @(varargin)covMaternard(3,varargin);
% scale = 1/5;
% inf_ = {@infGaussLik};
% mean_ = {@meanConst};
% lik_ = {@likGauss};
% cov_ = {@(varargin)covMaternard(3,varargin{:})};
% hyp.lik = log(1e-5);
% hyp.mean = 40;
% hyp.cov = log([(cond(:,2)-cond(:,1)).*scale;15]);
% acq = {@EI};
% xs = build_testP(cond,2e2);
% algo_data.start = 1;
% algo_data.l = 1;
% xt = [15;10;5;16];
% yt = [50; 35; 55;60];
% plot_post(hyp,inf_,mean_,cov_,lik_,xt,yt,xs,algo_data)
%%
xt = [0;2;4;5];
yt = [-1;1.5;-1.5;0];

xt = -100;
yt = 0;

cond = [0 5];

%f=@(x) (0.3375*x.^3-2.7125*x.^2+5.325*x-1);
f=@(x) (1343776/42771*x.^9-333476768/384939*x.^8+3918154000/384939*x.^7-25920084800/384939*x.^6+106835597842/384939*x.^5-285283315814/384939*x.^4+494455981150/384939*x.^3-59683351490/42771*x.^2+12306893312/14257*x-3306812579/14257).*(x>=1.5&x<=3.5)+0*x;
xs = linspace(0,5,400);
L=max(abs(diff(f(xs))/(xs(2)-xs(1))));
inf_ = {@infGaussLik};
mean_ = {@meanConst};
lik_ = {@likGauss};
cov_ = {@(varargin)covMaternard(3,varargin{:})};
hyp.mean = [10];
hyp.lik = log(2);
hyp.cov = log([(cond(:,2)-cond(:,1))/5; 8]);
acq = {@EI};
algo_data.start = 1;
algo_data.l = 1;
algo_data.f = f;
opts.plot = 1;
opts.acqFunc.xi = 0;
opts.acqFunc.beta = 2;
opts.termCondAcq = 0.05;
opts.safeOpt = 1;
opts.newSafeOpt = 1;
opts.safeOpts.threshold = 10;
opts.safeOpts.thresholdOffset = 100;
opts.safeOpts.searchCond = 4;
opts.safeOpts.thresholdPer = 0.2;
opts.safeOpts.thresholdOrder = 1;
opts.safeOpts.Optimize = 0;
opts.trainGP.train = 0;
opts.trainGP.acqVal = 5;
opts.samples = 4e2;
opts.dir_timeData = "/home/jannis/master/master_thesis/matlab/own_lab/test";
bayesOptima(hyp,inf_,mean_,cov_,lik_,acq,f,cond,opts,0.75)
%plot_post(hyp,inf_,mean_,cov_,lik_,xt,yt,build_testP(cond,100),algo_data)


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

function plot_post(hyp,inf_,mean_,cov_,lik_,xt,yt,xs,algo_data)
    start = algo_data.start;
    l = algo_data.l;
    D = length(l);
    if size(xt,2) ~= D
        x_vec = repmat(algo_data.x_vec,[size(xs,1),1]);
        x_vec(:,l) = xs;
        xs = x_vec;
    end
    [mu,var,~,~] = gp(hyp,inf_,mean_,cov_,lik_,xt,yt,xs);
    se = 2*sqrt(var);
    if D == 2
        fig=figure(2);
        plot3(xs(:,l(1)),xs(:,l(2)),se+mu,xs(:,l(1)),xs(:,l(2)),-se+mu)
        %fill3([Xs(:,1);flip(Xs(:,1))],[Xs(:,2);flip(Xs(:,2))],[se+mu;flip(-se+mu,1)],[9 9 9]/10)
        hold on
        %plot(xs,mu)
        plot3(xt(start:end-1,l(1)),xt(start:end-1,l(2)),yt(start:end-1),'k*', 'MarkerSize',10)
        plot3(xt(end,l(1)),xt(end,l(2)),yt(end,:),'k*','Color','g', 'MarkerSize',20)
        xlabel Kp
        ylabel Ki
%         if any(yt > 150)
%             zlim([0,150])
%         end
        hold off
    else
        T = 1;
        fig=figure(2);
        p1=fill([xs(:,l);flip(xs(:,l),1)],[se+mu;flip(-se+mu,1)],[9 9 9]/10);
        hold on
        p2=plot(xs(:,l),mu, 'Color','r');
        p3=plot(xt(:,l),yt(:),'k*', 'MarkerSize',10);
        %p4=plot(xt(end,l),yt(end),'k*','Color','r', 'MarkerSize',15);
        S = xs(se+mu < T);
        p4=plot(S(2:end-1),ones(length(S)-2,1)*T+5,'gs','MarkerFaceColor','g');
        y_vert = T:0.1:T+2;
        plot(S(1)*ones(size(y_vert)),y_vert,'k:',S(end)*ones(size(y_vert)),y_vert,'k:')
        p5=plot(S([2,end-1]),ones(2,1)*T+7,'bs','MarkerFaceColor','b');
        [u_max,I] = min(se+mu);
        M = xs(mu-se <= u_max);%-se(I)+mu(I) & se+mu<T);
        p6=plot(M(2:end-1),ones(size(M(2:end-1,1)))*T+9,'rs','MarkerFaceColor','r');
        yline(T,'--','threshold','LineWidth',2)
        legend([p1,p2,p3,p4,p5,p6],"$$2\sigma$$ range","mean","test point","$$\mathcal{S}$$","$$\mathcal{G}$$","$$\mathcal{M}$$", "Interpreter",'latex', 'location', 'northoutside','NumColumns',8)
        xlabel('X','interpreter','latex')
        ylabel('f','interpreter','latex')
        ylim([-3,3])
        
        hold off
    end
end