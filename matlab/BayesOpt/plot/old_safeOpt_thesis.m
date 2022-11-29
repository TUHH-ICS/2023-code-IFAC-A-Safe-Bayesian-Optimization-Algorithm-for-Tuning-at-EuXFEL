clear all
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'defaultTextInterpreter','latex')
set(groot, 'defaultAxesFontsize',15)
set(groot, 'defaultTextFontsize',20)
set(groot, 'defaultLegendFontSize',14)

cond = [0 5];

f=@(x) 0.3375*x.^3-2.7125*x.^2+5.325*x-1.0;
xs = linspace(0,5,400);
L=max(abs(diff(f(xs))/(xs(2)-xs(1))));
inf_ = {@infGaussLik};
mean_ = {@meanConst};
lik_ = {@likGauss};
cov_ = {@covSEard};
hyp.mean = [0.5];
hyp.lik = log(1e-5);
hyp.cov = log([(cond(:,2)-cond(:,1))/7; 0.75]);
acq = {@EI};
algo_data.start = 1;
algo_data.l = 1;
algo_data.f = f;
opts.plot = 1;
opts.acqFunc.xi = 0;
opts.acqFunc.beta = 2;
opts.termCondAcq = 0.01;
opts.safeOpt = 1;
opts.safeOpts.threshold = 1.3;
opts.safeOpts.thresholdOffset = 5;
opts.safeOpts.thresholdPer = 0.2;
opts.safeOpts.thresholdOrder = 1;
opts.safeOpts.Optimize = 0;
opts.trainGP.train = 0;
opts.samples = 4e2;
opts.dir_timeData = "/home/jannis/master/master_thesis/matlab/own_lab/test";

xs = build_testP(cond,400);
xt = [xs(201)];
yt = f(xt);
S = xt;
%%
[S,xn] = plot_post(hyp,inf_,mean_,cov_,lik_,xt,yt,xs,algo_data,opts,S);
xt = [xt;xn];
yt = f(xt);

function [S,xn] = plot_post(hyp,inf_,mean_,cov_,lik_,xt,yt,xs,algo_data,opts,S)
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
        hold on
        plot3(xt(start:end-1,l(1)),xt(start:end-1,l(2)),yt(start:end-1),'k*', 'MarkerSize',10)
        plot3(xt(end,l(1)),xt(end,l(2)),yt(end,:),'k*','Color','g', 'MarkerSize',20)
        set(gca,'TickLabelInterpreter','latex');
        xlabel("subdimension 1",'Interpreter','latex',"FontSize",15)
        ylabel("subdimension 2",'Interpreter','latex',"FontSize",15)
        hold off
    else
        fig=figure(2);
        if opts.safeOpt
            se = opts.acqFunc.beta*sqrt(var);
%             s = sprintf("$$%d%ssigma$$ confidence",opts.acqFunc.beta,"\");
            s = "$$\beta^{1/2}\sigma$$ confidence";
        else
            s = "$$2\sigma$$ confidence";
        end
        p1 = fill([xs(:,l);flip(xs(:,l),1)],[se+mu;flip(-se+mu,1)],[9 9 9]/10);
        hold on
        %p2 = plot(xs(:,l),mu, 'Color','r');
        %title('Posterior','Interpreter','latex')
        p3 = plot(xt(1:end,l),yt(1:end),'k*', 'MarkerSize',10);
        [yopt,I]=min(yt);

        algo_data.f=@(x) 0.3375*x.^3-2.7125*x.^2+5.325*x-1.0;
        f = algo_data.f(xs);
        p4 = plot(xs,f,'k-.');
        T = opts.safeOpts.threshold;
        Y = min(abs(xs-S'),[],2)*5.2911+mu+se;
        S = union(S,xs(Y <= T),"rows");
        [G,Ia,Ib] = setxor(xs,S);
        d_xx = abs(S-G');
        log_ar = mu(Y <= T)-se(Y <= T)+5.2911*d_xx <= T;
        u_min =  min(mu(Y <= T)+se(Y <= T));
        IM = mu(Y <= T)-se(Y <= T) <= u_min;
        IG = any(log_ar,2);
        M = S(IM);
        G = S(IG);
        se_S = se(Y <= T);
        if max(se_S(IG))>=max(se_S(IM))
            [~,Ixn] = max(se_S(IG));
            xn = G(Ixn);
        else
            [~,Ixn] = max(se_S(IM));
            xn = M(Ixn);
        end
        
        p8 = plot(xs,Y,'b--');
        p5 = plot(S(4:end-2),ones(length(S)-5,1)*T+0.5,'gs','MarkerFaceColor','g');
        if length(G)>14
            p7 = plot(G(4:end-2),ones(length(G)-5,1)*T+1,'bs','MarkerFaceColor','b');
        else
            p7 = plot(G+0.03,ones(length(G),1)*T+1,'bs','MarkerFaceColor','b');
        end
        if length(M)>3
            p6=plot(M(3:end-2),ones(size(M(3:end-2,1)))*T+1.5,'rs','MarkerFaceColor','r');
        else
            p6=plot(M,ones(size(M))*T+1.5,'rs','MarkerFaceColor','r');
        end
        p2 = plot(xs(:,l),mu, 'Color','r');
        y_vert = T:0.1:T+2;
        plot(min(S)*ones(size(y_vert)),y_vert,'k:',max(S)*ones(size(y_vert)),y_vert,'k:');

        if opts.safeOpt
            yline(opts.safeOpts.threshold,'--','threshold','LineWidth',2,'FontSize',13);
        end
        
        legend([p1,p2,p3,p4,p5,p6,p7,p8],s,"mean","test point","target function",'$$\mathcal{S}$$','$$\mathcal{M}$$','$$\mathcal{G}$$',"L", "Interpreter",'latex', 'location', 'northoutside','NumColumns',4,"FontSize",14)
        
        xlabel('$$x_*$$','interpreter','latex',"FontSize",15)
        ylabel('$f(x)$','interpreter','latex',"FontSize",15)
        set(gca,'TickLabelInterpreter','latex');
        hold off
        ylim([-1.75,3.25]);
        %exportgraphics(fig,"plots/"+sprintf('post_Safe_Nr_%d.pdf',length(xt)),"ContentType","vector")
        %pause(2)
    end
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
