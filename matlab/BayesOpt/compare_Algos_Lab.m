clear all

addresses={'LAB.SYNC/LASER.LOCK/26A1.L3.MENHIR/ADV_CTRL_MANAGER.0.PID.1.P_PARAM.WR',...
    'LAB.SYNC/LASER.LOCK/26A1.L3.MENHIR/ADV_CTRL_MANAGER.0.PID.1.I_PARAM.WR',...
    'LAB.SYNC/LINK.LOCK/26A.AMC7.CONTROLLER/LSU.1.SS_PID.P_PARAM.WR',...
    'LAB.SYNC/LINK.LOCK/26A.AMC7.CONTROLLER/LSU.1.SS_PID.I_PARAM.WR',...
    'LAB.SYNC/LASER.LOCK/26A2.L2.ORIGAMI10/ADV_CTRL_MANAGER.0.PID.2.P_PARAM.WR',...
    'LAB.SYNC/LASER.LOCK/26A2.L2.ORIGAMI10/ADV_CTRL_MANAGER.0.PID.2.I_PARAM.WR',...
    'LAB.SYNC/LINK.LOCK/26A.AMC7.CONTROLLER/LSU.0.CTRL_IN.STD_DEV.RD'};

lock_status={'LAB.SYNC/LASER.LOCK/26A2.L2.ORIGAMI10/LOCK_STATUS.VALUE.RD',...
    'LAB.SYNC/LASER.LOCK/26A1.L3.MENHIR/LOCK_STATUS.VALUE.RD',...
    'LAB.SYNC/LINK.LOCK/26A.AMC7.CONTROLLER/LSU.1.LOCK_STATUS.VALUE.RD'};

jitter_addr={'LAB.SYNC/LASER.LOCK/26A1.L3.MENHIR/CURRENT_INPUT_JITTER.RD',
    'LAB.SYNC/LINK.LOCK/26A.AMC7.CONTROLLER/LSU.1.TIMING_JITTER_FS.RD',
    'LAB.SYNC/LASER.LOCK/26A2.L2.ORIGAMI10/CURRENT_INPUT_JITTER.RD'};

cond_t=[
    -10 -0.31;
    -0.031 0;
    0 1.2;
    0.00001 0.21;
    -0.12 -0.0145;
    -0.00026 -0.000008
    ];

cond = repmat([-1,1],size(cond_t,1),1);
scale = [1/7;1/7];
scale = repmat(scale,[3,1]);

inf_ = {@infGaussLik};
mean_ = {@meanConst};
lik_ = {@likGauss};
cov_ = {@(varargin)covMaternard(3,varargin{:})};
hyp.lik = log(0.1);
hyp.mean = 11;
hyp.cov = log([(cond(:,2)-cond(:,1)).*scale;7]);
acq = {@EI};

x0 = [-4,-0.002,0.2,0.002,-0.04,-0.0001];
%x0 = [-4,-0.002,0.4,0.002,-0.04,-0.0001];
opts.plot=1;
opts.minFunc.mode=2;
opts.acqFunc.maxProb = 0;
opts.acqFunc.xi = 0.0;
opts.acqFunc.beta = 2;
opts.trainGP.acqVal = 10;%0.055;%0.5 %1D       %%% 0.05 D=1 with EI; 0.5 D = 1
opts.termCondAcq = 0.05;%0.05;%0.25;%0.5 %1D    %%% 0.05 D=1 with EI; 0.25 or 0.5 D=1; 0.2 D=2 with EI sf = 5
opts.maxIt = 50;
opts.trainGP.It = 10000;
opts.trainGP.train = 1;
opts.safeOpt = 1;
opts.newSafeOpt = 1;
opts.safeOpts.threshold = 11;
opts.safeOpts.thresholdOffset = 3;
opts.safeOpts.searchCond = 1;
opts.safeOpts.thresholdPer = 0.3;
opts.safeOpts.thresholdOrder = 1;
opts_lBO.maxIt = 100;
opts_lBO.sharedGP = true;
opts_lBO.subspaceDim = 1;
%opts_lBO.dim_combinations = [3;4;5;6];
%opts_lBO.oracle = 'descent';
opts_lBO.alpha = 0.02;
opts.dir_timeData="/home/luebsen/master/master_thesis/matlab/own_lab/data";
meas_mean = 1;
dir_ = "/home/luebsen/master/master_thesis/matlab/own_lab/data";
stepsize = 1 * (cond_t(:,2)-cond_t(:,1));
stepsize(stepsize > 1) = 1;

fun = @(params) LAB_readWrite(params, addresses, 4, 11, lock_status, jitter_addr,stepsize,1,2656,cond_t,dir_, meas_mean);

X0=[-3.39324394886329,-0.00935872887441044,0.799833496083104,0.0374160340995269,-0.0716150568183287,-8.23174052013711e-06;
    -4.48273891323462,-0.0250813654476746,0.513903591575263,0.101229812596080,-0.0720999393290359,-0.000111444113856845;
    -5.48680442245460,-0.0190768071445483,0.699583659297209,0.0528867676579610,-0.0609761364898579,-0.000104493097132908;
    -5.22527071337090,-0.00544433573264612,1.17919607966634,0.153354943888273,-0.0574760562304686,-0.000112814528013830;
    -4,-0.002,0.2,0.002,-0.04,-0.0001];

counter = 0;
save(dir_+'/counter_unlock.mat','counter');
time = 0;
save(opts.dir_timeData+'/time_data.mat','time');
jitter_data = [];
save(dir_+'/jitter_data.mat','jitter_data');

data = cell(5,5);
for i = 1:size(X0,1)
    x0 = forwardCoordTransf(cond_t,X0(i,:));
    [xopt,X,Y,DIM]=lineBO(hyp,inf_,mean_,cov_,lik_,acq, fun,cond,opts,opts_lBO,x0);
    data{i,1}=X;
    data{i,2}=Y;
    data{i,3}=DIM;
    jitter_data_str = load(dir_+'/jitter_data.mat');
    data{i,4}=jitter_data_str.jitter_data;
    count_str=load(dir_+'/counter_unlock.mat');
    counter = count_str.counter;
    data{i,5} = counter;
    time_str = load(opts.dir_timeData+'/time_data.mat');
    data{i,6} = time_str.time;
    counter = 0;
    save(dir_+'/counter_unlock.mat','counter');
    time = 0;
    save(opts.dir_timeData+'/time_data.mat','time');
    jitter_data = [];
    save(dir_+'/jitter_data.mat','jitter_data');

end
%save("Lab_data_dim1_random.mat",'data')

%%
% x0 = [];
% while size(x0,1) < 4
%     randx = rand(1,size(cond,1));
%     temp = bsxfun(@plus,bsxfun(@times,randx,(cond(:,2)'-cond(:,1)')),cond(:,1)');
%     if fun(temp) <= 40
%         x0(end+1,:) = temp;
%     end
% end
%%
%fun(x(I(1),:))
fun(X0(end,:))
%fun([-2.26859, -0.00157, 0.01212, 0.00213, -0.10082, -0.00000])
%xopt = [-4.00000, -0.00094, 0.20000, 0.00200, -0.06879, -0.00001];
%setParams(xopt,addresses,5,11)

% function varargout= setParams(params,addr,duration,freq, lock_addr, jitter_addr,stepsize,t)
%     D = length(params);
%     params_old = zeros(1,D);
%     for i = 1:D
%         data_str = doocsread(addr{i});
%         params_old(i)=round(data_str.data,6);
%     end
%     counter=1:D;
%     counter=counter(abs(params_old-params) > 10e-6);
%     %disp(counter)
%     for i=counter
%        writeData(addr{i},params(i),params_old(i),stepsize,t)
%     end
%     pause(0.2)
%     %pause(2)
%     [data,jitters] = readData(addr{end},duration,freq,lock_addr,jitter_addr);
%     varargout{1} = data;
%     if nargout == 2
%        varargout{2} = jitters;
%     end
% end
% 
% function writeData(addr,param,param_old,stepsize,t)
%     num_step=fix((param-param_old)/stepsize);
%     for i = 1:abs(num_step)
%         if num_step > 0
%             m=doocswrite(addr,param_old+i*stepsize);
%         else
%             m=doocswrite(addr,param_old-i*stepsize);
%         end
%         pause(t)
%     end
%     m=doocswrite(addr,param);
% end
% 
% function [data,jitters]=readData(addr,duration,freq,lock_addr,jitter_addr)
%     times = 1/freq;
%     time = duration;
%     iter = ceil(time/times);
%     jitter = zeros(iter,1);
%     jitters = zeros(iter,length(jitter_addr));
%         for i = 1:iter
%              if ~checkLockStatus(lock_addr)
%                  st=tic;
%                 while ~checkLockStatus(lock_addr)
%                     pause(3)
%                     if toc(st) > 31
%                         s = input("Can't lock system, set value manually by writing 'user' or continue by pressing 'enter': ",'s');
%                         if strcmp(s,'user')
%                             data = input("Set value: ");
%                             jitters=100*ones(1,length(jitter_addr)+1);
%                         else
%                             st = tic;
%                             continue;
%                         end
%                         return
%                     end
%                 end
%             end
%            data_str = doocsread(addr);
%            jitter(i) = data_str.data*10e3;
%            for l = 1:length(jitter_addr)
%                temp_str = doocsread(jitter_addr{l});
%                jitters(i,l) = temp_str.data;
%            end
%            pause(times)
%         end
%     data = mean(jitter);
%     std_dev = std(jitter);
%     jitters = [mean(jitters,1),std_dev];
%     disp(data)
%     disp(jitters)
%     disp(std_dev)
% end
% 
% function b=checkLockStatus(lock_addr)
%     c = 0;
%     for i =1:length(lock_addr)
%         c_str = doocsread(lock_addr{i});
%         c = c + c_str.data;
%     end
%     if c == 0
%         b = true;
%     else
%         b = false;
%     end
% end  