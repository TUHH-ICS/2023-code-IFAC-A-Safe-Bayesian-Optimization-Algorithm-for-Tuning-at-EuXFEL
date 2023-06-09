%---------------------------------------------------------------------------------------------
% For Paper,
% "A Safe Bayesian Optimization Algorithm for Tuning the Optical Synchronization System at European XFEL"
% by Jannis O. Lübsen, Maximilian Schütte, Sebastian Schulz, Annika Eichler
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Jannis Lübsen
%--------------------------------------------------------------------------------------------
% This is the main file of the Bayesian Optimization implementation
% struct hyp includes the hyperparameters
% inf_ the inference method
% mean_ the mean function
% cov_ covariance function
% lik_ likelihood distribution
% acq_func handle to acquisition function
% obj_func handle to the evaluation function
% cond are lower and upper boundaries for every parameter [lower, upper]
% opts is a struct with several options
% varargin can be used to define an x0 or other options as LineBO

function [xopt,yopt,xt,yt]=bayesOptima(hyp,inf_,mean_,cov_,lik_,acq_func,obj_func,cond,opts,varargin)

%%% standard options %%%%
oldOpts.plot = 0;  % defines if the predictive distributions, acq function etc shall be plotted
oldOpts.maxIt = 40; % max iterations
oldOpts.termCondAcq = 5e-4; % terminates BO if max(acqFunc(x)) < termCondAcq
oldOpts.safeOpt = 0; % use SafeOpt
oldOpts.moSaOpt = 0; % use MoSaOpt
oldOpts.maxProb = 0; % maximization problem
oldOpts.samples = 1e2; % samples used for the plot or if not in safeOpts defined also for safe options

trainGP.acqVal = 0.5; % max(acqFunc(x)) < acqVal the GP model hyperparameter are fitted (not if (modified) safe options are used)
trainGP.It = 40; % every 40 iterations the hyp params are fitted (not if (modified) safe options are used)
trainGP.train = false; % defined whether the hyperparameters should be fitted
trainGP.coolDown = 2; % minimum number of iterations after the hyperparameters can be fitted again
oldOpts.trainGP = trainGP;
skip_dim = 0; % skip dimension (debug: pause execution and set value to 1)


if isempty(opts) || ~isstruct(opts)
    opts = oldOpts;
else
    opts=getopts(oldOpts,opts);
end

fun_in = obj_func;
coolDown = opts.trainGP.coolDown;
gp_cd = 0;
checkGP_train = false;
maxProb = opts.maxProb;
opts.acqFunc.maxProb = maxProb;

if maxProb
    obj_eval = @max;
else
    obj_eval = @min;
end

D = size(cond,1);
xt = zeros(opts.maxIt,D);
yt = zeros(opts.maxIt,1);
x_new = zeros(1,D);
algo = 'naive';


if length(varargin) < 1 || isempty(varargin{1})
    randx = rand(1,D);
    xt(1,:) = bsxfun(@plus,bsxfun(@times,randx,(cond(:,2)'-cond(:,1)')),cond(:,1)');
    yt(1,:)= fun_in(xt(1,:));
end
if length(varargin) == 1
    print_vector("X0 = ", varargin{1})
    xt(1,:)=varargin{1};
    yt(1,:) = fun_in(xt(1,:));
end

if length(varargin) == 2
    AlgoStruct = varargin{2};
    if ~isstruct(AlgoStruct)
        error("varargin{2} need to be struct")
    end
    algo = AlgoStruct.algo;
end

switch algo
    case 'lineBO'
        x_vec = varargin{1};
        l = AlgoStruct.subspace;
        xt = AlgoStruct.post.x;
        yt = AlgoStruct.post.y;
        n = size(xt,1);
        counter = n;
        xt = [xt;zeros(opts.maxIt,D)];
        yt = [yt;zeros(opts.maxIt,1)];
        cond_acq = [cond(l,1), cond(l,2)];
        D = length(l);
        algo_data.l = l;
        algo_data.x_vec = x_vec;
        algo_data.start = counter;
        algo_data.name = 'lineBO';
        algo_data.plot = opts.plot;
        algo_data.cond = cond;
        algo_data.y0 = AlgoStruct.post.yopt;
        algo_data.x0 = AlgoStruct.post.xopt;
    otherwise
        algo_data.plot = opts.plot;
        counter = 1;
        algo_data.l = (1:D)';
        cond_acq = cond;
        algo_data.start = 1;
        algo_data.name = 'naive';
end
algo_data.observed = 0;

fprintf("y_start = %.5f", obj_eval(yt(1:counter)));
if D > 2 && opts.plot
    error("Plot only for (subspace) dim <= 2 supported")
end

if (opts.safeOpt || opts.moSaOpt) && D <= 2
    if ~isfield(opts,'safeOpts')
        error("At least the threshold of the safe options must be defined: opts.safeOpts.threshold=...")
    end
    if ~isfield(opts.safeOpts,'samples')
        opts.safeOpts.samples = opts.samples;
    end
    algo_data.observed = false;
    opts.safeOpts.xs_safe = build_testP(cond_acq, opts.safeOpts.samples);
elseif (opts.safeOpt || opts.moSaOpt) && D > 2
    error("safeOpt only for (subspace) dim <= 2 supported")
else
    xs_safe = [];
end

if opts.moSaOpt
    opts.trainGP.train = 1;
end
opts.safeOpts.useSafeOpt = opts.safeOpt;
opts.safeOpts.useMoSaOpt = opts.moSaOpt;


if D <= 2
    xs = build_testP(cond_acq, opts.samples);
else
    xs = [];
end

opts.algo_data = algo_data;

for i=1:opts.maxIt
    algo_data.it = i;
    fprintf("\niteration: %d\n", i);
    start = tic;
    try
        [x_new, min_nacq, algo_data] = find_new_parameters(hyp,inf_,mean_,cov_,lik_,xt(1:counter,:),yt(1:counter,:),acq_func,cond_acq,opts,algo_data,xs);
    catch ME
        if strcmp(ME.identifier,'MATLAB:posdef') || strcmp(ME.identifier, 'MATLAB:unassignedOutputs')
            fprintf("\nCovariance Matrix lost full rank!\nAbort BO and return optimal values!\n\n")
            print_vector("xnew = ", x_new)
            break;
        else
            rethrow(ME)
        end
    end
    if (abs(min_nacq) < opts.trainGP.acqVal || mod(i,opts.trainGP.It)==0)  && opts.trainGP.train && (algo_data.observed || ~opts.safeOpt)
        checkGP_train = true;
    end

    if checkGP_train && gp_cd == 0
        disp(newline+"cov hyp: "+num2str(hyp.cov')+"    lik hyp: "+num2str(hyp.lik))
        [hyp,~,~]=fitGP(hyp,inf_,mean_,cov_,lik_,xt(1:counter,:),yt(1:counter,:),opts,algo_data);
        gp_cd = coolDown;
        disp(newline+"cov hyp: "+num2str(hyp.cov')+"    lik hyp: "+num2str(hyp.lik))
        try
            [x_new, min_nacq, algo_data] = find_new_parameters(hyp,inf_,mean_,cov_,lik_,xt(1:counter,:),yt(1:counter,:),acq_func,cond_acq,opts,algo_data,xs);
        catch ME
            if strcmp(ME.identifier,'MATLAB:posdef') || strcmp(ME.identifier, 'MATLAB:unassignedOutputs')
                fprintf("\nCovariance Matrix lost full rank!\nAbort BO and return optimal values!\n\n")
                print_vector("xnew = ", x_new)
                break;
            else
                rethrow(ME)
            end
        end
    end

    print_vector("new parameter: ", x_new)

    if abs(min_nacq) < opts.termCondAcq || skip_dim
        break;
    end

    xt(counter+1,:) = x_new;
    yt(counter+1,:)=fun_in(xt(counter+1,:));
    counter = counter + 1;
    gp_cd = gp_cd - 1;
    if gp_cd < 0
        gp_cd = 0;
    end

    stop=toc(start);
    disp("time per interation: "+num2str(stop))
    if isfield(opts,'dir_timeData') && ~isempty(opts.dir_timeData)
        time_str = load(opts.dir_timeData+'/time_data.mat',"time");
        time = time_str.time;
        time(end+1)=stop+time(end);
        save(opts.dir_timeData+'/time_data.mat',"time")
    end
end


% crop array at true length
xt = xt(1:counter,:);
yt = yt(1:counter);

% find optimal parameters
[~,~,mu,~]=gp(hyp,inf_, mean_, cov_, lik_, xt(algo_data.start+1:end,:), yt(algo_data.start+1:end),xt(algo_data.start+1:end,:));
beta = 0.3;
[yopt,I] = obj_eval((1-beta)*yt(algo_data.start+1:end)+beta*(mu));
if yopt < algo_data.y0
    xopt = xt(algo_data.start+I,:);
else
    xopt=algo_data.x0;
    yopt = algo_data.y0;
end
print_vector("x0 = ", algo_data.x0)
print_vector("xopt = ", xopt)
print_vector("y0 = ",algo_data.y0)
print_vector("yopt = ", yopt)
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

function [x_new, min_nacq, algo_data]=find_new_parameters(hyp,inf_,mean_,cov_,lik_,xt,yt,acq_func,cond,opts,algo_data,xs)
l = algo_data.l;
if isscalar(l)
    [x_new, min_nacq, nacq, algo_data]=find_maximum(hyp,inf_,mean_,cov_,lik_,xt,yt,acq_func,cond,opts,algo_data,xs);
    if algo_data.plot
        plot_post(hyp,inf_,mean_,cov_,lik_,xt,yt,xs,algo_data,opts)
        if ~opts.safeOpt
            figure(3);
            plot(xs,-nacq,'Color','b')
            hold on
            p=plot(x_new(l),-min_nacq,'k*','Color','r', 'MarkerSize',15);
            legend(p, "$$\mbox{max}(\alpha)$$", "interpreter", 'latex', 'location', 'best',"FontSize",12)
            xlabel("$$x_*$$","Interpreter","latex","FontSize",14)
            ylabel("$$\alpha(x)$$", "Interpreter","latex","FontSize",14)
            set(gca,'TickLabelInterpreter','latex');
            hold off
        end
    end
else
    [x_new, min_nacq,algo_data]=find_maximum(hyp,inf_,mean_,cov_,lik_,xt,yt,acq_func,cond,opts,algo_data);
    if algo_data.plot
        plot_post(hyp,inf_,mean_,cov_,lik_,xt,yt,xs,algo_data,opts)
    end
end
fprintf("max acqusition value: %f\n", -min_nacq);
end

function plot_post(hyp,inf_,mean_,cov_,lik_,xt,yt,xs,algo_data,opts)
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
    figure(2);
    plot3(xs(:,l(1)),xs(:,l(2)),se+mu,xs(:,l(1)),xs(:,l(2)),-se+mu)
    hold on
    plot3(xt(start:end-1,l(1)),xt(start:end-1,l(2)),yt(start:end-1),'k*', 'MarkerSize',10)
    plot3(xt(end,l(1)),xt(end,l(2)),yt(end,:),'k*','Color','g', 'MarkerSize',20)
    set(gca,'TickLabelInterpreter','latex');
    xlabel("subdimension 1",'Interpreter','latex',"FontSize",15)
    ylabel("subdimension 2",'Interpreter','latex',"FontSize",15)
    hold off
else
    figure(2);
    if opts.safeOpt
        se = opts.acqFunc.beta*sqrt(var);
        s = sprintf("$$%d%ssigma$$ confidence",opts.acqFunc.beta,"\");
    else
        s = "$$2\sigma$$ confidence";
    end
    p1 = fill([xs(:,l);flip(xs(:,l),1)],[se+mu;flip(-se+mu,1)],[9 9 9]/10);
    hold on
    p2 = plot(xs(:,l),mu, 'Color','r');
    p3 = plot(xt(end,l),yt(end),'k*','Color','r', 'MarkerSize',15);
    if opts.safeOpt||opts.moSaOpt
        yline(opts.safeOpts.threshold,'--','threshold','LineWidth',2);
    end
    if start < size(xt,1)
        p8 = plot(xt(start:end-1,l),yt(start:end-1),'k*', 'MarkerSize',10);
        legend([p1,p2,p3,p8],s,"mean","new point","test point", "Interpreter",'latex', 'location', 'northoutside','NumColumns',4,"FontSize",14) %"target function",'$$\mathcal{S}$$','$$\mathcal{M}$$','$$\mathcal{G}$$',
    else
        legend([p1,p2,p3],s,"mean","new point", "Interpreter",'latex', 'location', 'northoutside','NumColumns',4,"FontSize",14) %"target function",'$$\mathcal{S}$$','$$\mathcal{M}$$','$$\mathcal{G}$$',
    end

    xlabel('$$x_*$$','interpreter','latex',"FontSize",15)
    ylabel('$f(x)$','interpreter','latex',"FontSize",15)
    set(gca,'TickLabelInterpreter','latex');
    hold off
end
end

function print_vector(msg,vector)
fprintf("%s[",msg)
if length(vector) > 1
    fprintf("%.5f, ",vector(1:end-1))
end
fprintf("%.5f]\n", vector(end))
end


