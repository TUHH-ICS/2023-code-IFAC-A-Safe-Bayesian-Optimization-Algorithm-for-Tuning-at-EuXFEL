function [x,y,c,L_f,L_g] = SMGO(cond_old,obj,g,x0,opts)

    oldOpts.maxIt = 250;
    oldOpts.L_f = 0;
    oldOpts.delta = 0.1;
    oldOpts.beta = 0.1;
    oldOpts.alpha = 0.05;
    oldOpts.mu = 1.4;
    oldOpts.B = 5;
    oldOpts.L_g = [];
    oldOpts.start_points = 0;
    oldOpts.safety.threshold = 10;
    oldOpts.samples = 100;
    oldOpts.plot = 0;
    oldOpts.coordTransform_sort = 0;
    oldOpts.coordTransform_interval = [-1 1];
    
    opts = getopts(oldOpts,opts);
    maxIt = opts.maxIt;
    S = length(g);
    D = size(cond_old,1);
    if isempty(opts.L_g)
        opts.L_g = zeros(1,S);
    end
    cond = repmat([-1 1],D,1);
    fTransf = @(x) forwardCoordTransf(cond_old,x,opts.coordTransform_sort,opts.coordTransform_interval);
    bTransf = @(x) backwardCoordTransf(cond_old,x,opts.coordTransform_sort,opts.coordTransform_interval);
    safeOpts = opts.safety;
    
    x0 = fTransf(x0);
  
    E = x0;
    E = buildSamplePoints(E,x0,opts.B,cond);
    L_f = opts.L_f;
    L_g = opts.L_g;
    x = zeros(maxIt,D);
    y = zeros(maxIt,1);
    c = zeros(maxIt,S);
    
    if opts.start_points
        for i = 1:size(x0,1)
            x(i,:) = x0(i,:);
            y(i) = obj(bTransf(x(i,:)));
            for k = 1:S
                c(i,k) = feval(g{k},bTransf(x(i,:)));
            end
        end
        start = size(x0,1);
        x_new = x0(end,:);
    else
        x_new = x0;
        start = 1;
    end
    
    if size(cond_old,1) <= 2
        xs = build_testP(cond, opts.samples);
    else
        warning("plots are disabled")
        opts.plot = 0;
    end
    
    for i = start:maxIt
        disp("Iteration: " + num2str(i))
        x(i,:) = x_new;
        y(i) = obj(bTransf(x(i,:)));
        disp("y = " + num2str(y(i)))
            
        [yopt,Iopt] = min(y(1:i));
        E = buildSamplePoints(E,x(1:i,:),opts.B,cond);
        xopt = bTransf(x(Iopt,:));
        disp("yopt = " + num2str(yopt) + newline +"xopt = " + num2str(xopt))
        for k = 1:S
            c(i,k) = feval(g{k},bTransf(x(i,:)));
        end
    
        if i > 1
            [L_f,L_g]=updateLipschitz(x(1:i,:),y(1:i),c(1:i,:),L_f,L_g);
        end

        disp("Lipschitz Objective = " + num2str(L_f))
        disp("Lipschitz Constraint = " + num2str(L_g))
    
        [f_app,g_app,f_unc,g_unc,f_u,f_l] = buildBounds(y(1:i),c(1:i,:),L_f,L_g,x(1:i,:),E,opts.mu);
        
        logic_arr = logical(prod(g_app >= 0,2));

        const_logic = logical(prod(logic_arr,2));
        
        [val,I]=min(f_app(const_logic)-opts.beta*f_unc(const_logic));
        E_l = E(const_logic,:);
        x_new = E_l(I,:);
        f_est = f_l(const_logic);
        f_est = f_est(I);
        
        if opts.plot
            plotObjective(x(1:i,:),y(1:i),xs,L_f,opts.mu,safeOpts,opts.alpha,yopt,cond_old)
            plotFirstConstraint(x(1:i,:),c(1:i,1),xs,L_g,opts.mu,cond_old)
            figure(5)
            if size(E,2) == 2, plot(E(:,1),E(:,2),'k*'); else, plot(E,ones(length(E),1),'k*'); end
        end

        if ~checkOptimizationPoint(f_est,yopt,opts.alpha,L_f)
            [f_app,g_app,f_unc,g_unc,f_u,f_l] = buildBounds(y(1:i),c(1:i,:),L_f,L_g,x(1:i,:),E,opts.mu);
            logic_arr = g_app >= 0;%, f_u <= safeOpts.threshold
            const_logic = logical(prod(logic_arr,2));

            disp("Exploration mode")
            w_f = zeros(size(f_unc));
            w_f(const_logic) = f_unc(const_logic);
            w_pi = sum(g_unc ./ L_g,2);
            
            w_g = prod(double(logic_arr) + ones(size(logic_arr)),2);
    
            [val2,I2] = max((1-opts.delta)*w_f+opts.delta*w_pi.*w_g);
            x_new = E(I2,:);
        else
            disp("Exploitation mode")
        end
        disp("new parameter xnew = " + num2str(bTransf(x_new)))
        fprintf("\n")
        pause(2)
    end
    x = bTransf(x);
end

function plotObjective(xk,yk,xs,L_f,mu,safeOpts,alpha,yopt,cond)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'defaultTextInterpreter','latex')
set(groot, 'defaultAxesFontsize',12)
set(groot, 'defaultTextFontsize',12)
set(groot, 'defaultLegendFontsize',12)
    [f_app,~,f_u,f_l] = buildObjectiveBounds(yk,L_f,xk,xs,mu);
    xk = backTransform(cond,xk);
    xs = backTransform(cond,xs);
    if size(xs,2) == 1
        fig=figure(1);
        fig.Units='centimeters';
        fig.Position(3:end)= [10.5,9];
        obj = @(x) 2*sin(1*pi*x(:,1))-cos(0.21*x(:,1))+x(:,1).^2;
        xss = linspace(-4,4,1000);
        p1 = plot(xss',obj(xss'),'k:','LineWidth',0.3);
        hold on
        p=plot(xs,f_app,'k-',xs,f_u,'b-',xs,f_l,'r-');
        p4=plot(xk(end,:),yk(end),'r*');
        yline(safeOpts.threshold,'--',"safety threshold","LineWidth",2)
        exploi_thresh = yopt-alpha*L_f;
        yline(exploi_thresh,'--',"exploitation threshold", "LineWidth",2)
        if size(xk,1) > 1
            p5=plot(xk(1:end-1,:),yk(1:end-1),'k*');
            legend([p',p4,p5,p1],"$$1/2(\overline{f}+\underline{f})$$","$$\overline{f}$$","$$\underline{f}$$","last point","test point","target function",'Interpreter','Latex',"NumColumns",2,'Location','north')
        else
            legend([p',p4,p1],"$$1/2(\overline{f}+\underline{f})$$","$$\overline{f}$$","$$\underline{f}$$","last point","target function",'Interpreter','Latex',"NumColumns",2);
        end
        ylim([min(min(f_l),exploi_thresh)*1.2,max(max(f_u),safeOpts.threshold)*1.5])
        hold off
        xlabel("$x_*$")
        ylabel("$J$ [fs]")
    else
        figure(1)
        plot3(xs(:,1),xs(:,2),f_u,xs(:,1),xs(:,2),f_l)
        hold on
        plot3(xk(end,1),xk(end,2),yk(end),'r*')
        if size(xk,1) > 1
            plot3(xk(1:end-1,1),xk(1:end-1,2),yk(1:end-1),'k*')
            legend("upper bound f","lower bound f","last point","test point",'Interpreter','Latex','NumColumns',2,'Location','best')
        else
            legend("upper bound f","lower bound f","last point",'Location','best')
        end
        hold off
    end
end

function plotFirstConstraint(xk,ck,xs,L_g,mu,cond)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'defaultTextInterpreter','latex')
set(groot, 'defaultAxesFontsize',12)
set(groot, 'defaultTextFontsize',12)
set(groot, 'defaultLegendFontsize',12)
    [g_app,~,g_u,g_l] = buildObjectiveBounds(ck,L_g,xk,xs,mu);
    fig=figure(2);
    fig.Units='centimeters';
    fig.Position(3:end)= [10.5,9];
    if size(xs,2) == 1
        xs = backTransform(cond,xs);
        p=plot(xs,g_u,'b',xs,g_l,'r',xs,g_app,'k');
        hold on
        yline(0,'--','LineWidth',2)
        legend("$$\overline{g}$$","$$\underline{g}$$","$$\tilde{g}$$",'NumColumns',3,'Location','north')
        ylim([min(g_l)*1.3,max(g_u)*1.5])
        xlabel("$x_*$")
        ylabel("g")
        hold off
    else
        plot3(xs(:,1),xs(:,2),g_app)
        legend("$$\tilde{g}$$",'Location','best')
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

function [x_transf]=forwardTransform(cond,x)
    mu_ = (cond(:,2)'+cond(:,1)')/2;
    slope_ = (cond(:,2)'-cond(:,1)')/2;
    x_transf = (x-mu_)./slope_;
end

function [x] = backTransform(cond,x_transf)
    mu_ = (cond(:,2)'+cond(:,1)')/2;
    slope_ = (cond(:,2)'-cond(:,1)')/2;
    x = (slope_.*x_transf)+mu_;
end