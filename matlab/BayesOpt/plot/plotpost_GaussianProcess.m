Kp_max = 5;
Kp_min = 0;
Ki_max = 6e1;
Ki_min = 0;
cond=[Kp_min, Kp_max];
      %Ki_min, Ki_max];

cov_matern = @(varargin)covMaternard(3,varargin);
scale = 1/5;
inf_ = {@infGaussLik};
mean_ = {@meanConst};
lik_ = {@likGauss};
cov_ = {@(varargin)covMaternard(3,varargin{:})};
hyp.lik = log(1e-5);
hyp.mean = 0;
hyp.cov = log([(cond(:,2)-cond(:,1)).*scale;1]);
acq = {@EI};
xs = build_testP(cond,1e2);
algo_data.start = 1;
algo_data.l = 1;
 xt = [2.5;1;0;3];
 yt = [1; -0.2; 0.1;-0.5];
%xt = 2.5;
%yt = 1;
f=@(x) 0.3375*x.^3-2.7125*x.^2+5.325*x-1.0;
algo_data.f = f;
plot_post(hyp,inf_,mean_,cov_,lik_,xt,yt,xs,algo_data)



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
%     se = 2*ones(size(var));
%     mu = zeros(size(mu));
%     mu
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
        f=algo_data.f(xs);
        T = 60;
        fig=figure(2);
        p1=fill([xs(:,l);flip(xs(:,l),1)],[se+mu;flip(-se+mu,1)],[9 9 9]/10);
        hold on
        p2=plot(xs(:,l),mu, 'Color','r','MarkerSize',2.5);
        %plot(xs,f,'k--')
        title('Prior','Interpreter','latex')
        xlabel('$$\mathcal{X}$$','interpreter','latex')
        ylabel('f(x)','interpreter','latex')
        c = mu + se/2.*mvnrnd(mu,feval(cov_{:},hyp.cov,xs),3)';
        p=plot(xs,c(:,1));
        plot(xs,c(:,2:3));
        plot(xt,yt,'k*', 'MarkerSize', 8)
        legend("$$2\sigma$$ range","mean","$$g_1$$","$$g_2$$","$$g_3$$","test point", "Interpreter",'latex', 'location', 'northeast','NumColumns',2)
        %% 
        %"f1","f2","f3",
        axes = gca;
%         axes.XTick=xs;
%         xl = cell(1,length(xs));
%         for i = 1:length(xs)
%             xl{i}=sprintf('$$x_{%d*}$$',i);
%         end
%     
%         axes.XTickLabel = xl;
        axes.TickLabelInterpreter = 'latex';
        ylim([-3 3])
        hold off
        figure(1)
        post = yt;
        opts.maxProb = 0;
        opts.xi = 0;
        acq=EI(xs,hyp,inf_,mean_,cov_,lik_,xt,post,yt,opts,[]);
        plot(xs,-acq);
        hold on 
        [min_nacq,I] = min(acq);
        p1=plot(xs(I),-min_nacq, 'r*')
        legend(p1,'$$\mbox{max}(\alpha(x))$$', 'interpreter','latex')
        xlabel('$$\mathcal{X}$$','interpreter','latex')
        ylabel('$$\alpha(x)$$','interpreter','latex')
        axes = gca;
        axes.TickLabelInterpreter = 'latex';
    end
end