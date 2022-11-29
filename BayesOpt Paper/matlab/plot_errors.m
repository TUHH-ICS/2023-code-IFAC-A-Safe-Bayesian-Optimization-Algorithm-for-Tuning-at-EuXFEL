clear all
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'defaultTextInterpreter','latex')
set(groot, 'defaultAxesFontsize',10)
set(groot, 'defaultTextFontsize',10)
set(groot, 'defaultLegendFontSize',8)

t = -5:0.05:5;
f=@(x)sin(3*pi*x).*exp(-0.3*(x).^2)+0.75*(x).^2;
% plot(t,f(t))
% hold on
% xt=-3 + (3+3)*rand(10,1);
xt=(-5:1/3:5)';
% plot(x,f(x),'+')
f2=@(x)0.75*(x).^2;
% hold off
% plot(t,f2(t))
yt=f(xt);

cond=[-3,3];

inf_ = {@infGaussLik};
mean_ = {@meanConst};
lik_ = {@likGauss};
cov_ = {@(varargin)covMaternard(3,varargin{:})};
hyp.mean = [5];
hyp.lik = log(0.00001);
hyp.cov = log([(cond(:,2)-cond(:,1))/7; 5]);

opts.plot = 1;
opts.termCondAcq = .005;
opts.samples = 1e3;
opts.minFunc.mode=2;
opts.minFunc.MaxFunEvals=1000;
opts.moSaOpt = 0;
opts.safeOpt = 0;

algo_data.start = 1;
algo_data.f = f;

hyp2 = gpTrain(hyp,inf_,mean_,cov_,lik_,xt,yt,opts,[]);
plot_post(hyp,inf_,mean_,cov_,lik_,xt,yt,build_testP(cond,opts.samples),algo_data,opts,hyp2)


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

function plot_post(hyp,inf_,mean_,cov_,lik_,xt,yt,xs,algo_data,opts,hyp2)
    start = algo_data.start;
    l = 1
    D = length(l);
    if size(xt,2) ~= D
        x_vec = repmat(algo_data.x_vec,[size(xs,1),1]);
        x_vec(:,l) = xs;
        xs = x_vec;
    end
    [mu,var,~,~,~,post] = gp(hyp,inf_,mean_,cov_,lik_,xt,yt,xs);
    v = post.L'\(repmat(post.sW,1,size(xs,1)).*feval(cov_{:},hyp.cov,xt,xs));
    K = feval(cov_{:},hyp.cov,xs)-v'*v;
    est_L = 100;
    max_dR=zeros(est_L,1);
    for i = 1:est_L
        R=mvnrnd(mu,K);
        dR=diff(R)/(xs(2)-xs(1));
        max_dR(i) = max(abs(dR));
    end
    max_dR = sort(max_dR,'ascend');
    L = mean(max_dR(1:1));
    [mu,var,~,~,] = gp(hyp2,inf_,mean_,cov_,lik_,xt,yt,xs);
    se = 2*sqrt(var);
%     mu = zeros(size(xs));
%     se = 2*ones(size(xs));
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
        fig = figure(1);
        fig.Units = 'points';
        fig.Position(3:4)=[250,150];
        if opts.safeOpt
            se = opts.acqFunc.beta*sqrt(var);
            s = sprintf("$$%2$smu(%2$smathbf{x}_*) %2$spm %1$d%2$ssigma(%2$smathbf{x}_*)$$",opts.acqFunc.beta,"\");
        else
            s = "$$\mu(\mathbf{x}_*) \pm 2\sigma(\mathbf{x}_*)$$";
        end
        p1 = fill([xs(:,l);flip(xs(:,l),1)],[se+mu;flip(-se+mu,1)],[9 9 9]/10);
        hold on
        p2 = plot(xs(:,l),mu, 'Color','r');
        %title('Posterior','Interpreter','latex')
%         p3 = plot(xt(187,l),yt(187),'r*');
%         [yopt,I]=min(yt);
% 
%         T = 1;
%         algo_data.f=@(x) 0.3375*x.^3-2.7125*x.^2+5.325*x-1.0;
        f = algo_data.f(xs);
        p4 = plot(xs,f,'k-.');
        xk = permute(xt,[3,2,1]);
        x_norm=squeeze(vecnorm(xs-xk,2,2))';
        f_u = max(yt-L*x_norm,[],1)';
        p5=plot(xs,f_u,':','Color',[0.8500 0.3250 0.0980]);
        val2=[];
        xc=[];
%         ranii = randperm(19); 
%         xt=xt(ranii);
%         yt=yt(ranii);
        for i = 1:length(xt)
%             ranii(i)
            val =[];
            for j = 1:length(xt)
                if j==i
                    continue;
                end
                if xt(i)>= xt(j)
                    val(end+1)=abs(xt(j)-xt(i))/2+(yt(i)-yt(j))/(2*L);
                else
                    val(end+1)=abs(xt(j)-xt(i))/2+(yt(j)-yt(i))/(2*L);
                end
%                 val(end+1)=abs(yt(i)+yt(j))/L;

            end
            [~,I]=min(abs(val));
            xc(end+1) = val(I)-xt(i);
            val2(end+1) = yt(i) - (val(I))*L;
        end
        plot(xc,val2,'r+')
        min(val2)

%           
%         S = xs(se+mu < T);
%         [u_max,I] = min(se+mu);
%         y_vert = T:0.1:T+2;
%         M = xs(mu-se <= u_max & se + mu < T);%-se(I)+mu(I) & se+mu<T);
%         p5 = plot(S(3:end-2),ones(length(S)-4,1)*T+0.5,'gs','MarkerFaceColor','g');
%         if length(M)>10
%             p6=plot(M(3:end-2),ones(size(M(3:end-2,1)))*T+1.5,'rs','MarkerFaceColor','r');
%         else
%             p6=plot(M(2:end-1),ones(size(M(2:end-1,1)))*T+1.5,'rs','MarkerFaceColor','r');
%         end
%         p7 = plot(S([3,end-2]),ones(2,1)*T+1,'bs','MarkerFaceColor','b');
% 
%         plot(S(1)*ones(size(y_vert)),y_vert,'k:',S(end)*ones(size(y_vert)),y_vert,'k:');
        if opts.safeOpt
            yline(opts.safeOpts.threshold,'--','threshold','LineWidth',2);
        end
         %yline(yopt,'b:',"$f_n^*$",'Interpreter','latex',"FontSize",15,"LineWidth",1.5)
       if start < size(xt,1)
            p8 = plot(xt(start:end,l),yt(start:end),'b*');
%             p3 = plot(xt(187,l),yt(187),'r*');
            legend([p1,p2,p8,p4,p5],s,"$\mu(\mathbf{x}_*)$","$y_i$","$f(\mathbf{x})$",'$\tau(\mathbf{x})$','$$\mathcal{S}$$','$$\mathcal{M}$$','$$\mathcal{G}$$', "Interpreter",'latex', 'location', 'north','NumColumns',2,"FontSize",8) %"target function",'$$\mathcal{S}$$','$$\mathcal{M}$$','$$\mathcal{G}$$',
%         else
%            legend([p1,p2],s,"$\mu(\mathbf{x}_*)$","$y_i$","$f(\mathbf{x})$",'$$\mathcal{S}$$','$$\mathcal{M}$$','$$\mathcal{G}$$', "Interpreter",'latex', 'location', 'northoutside','NumColumns',3,"FontSize",8) %"target function",'$$\mathcal{S}$$','$$\mathcal{M}$$','$$\mathcal{G}$$',
       end
       
        
        xlabel('$$\mathbf{x}_* \in \mathcal{X}$$','interpreter','latex',"FontSize",10)
        ylabel("$f_* \vert \mathbf{y},X,X_*$","FontSize",10)
        xticks(-3:1:3)
        xlim([-3,3])
        set(gca,'TickLabelInterpreter','latex');
        grid on
%         ylim([9,17]);
        %exportgraphics(fig,"plots/"+sprintf('post_Safe_Nr_%d.pdf',length(xt)),"ContentType","vector")
        %pause(2)
    end
end