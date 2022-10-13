clear all

addresses={'XFEL.SYNC/LASER.LOCK.XLO/XTIN.MLO1/ADV_CTRL_MANAGER.0.PID.2.P_PARAM.WR',...
     'XFEL.SYNC/LASER.LOCK.XLO/XTIN.MLO1/ADV_CTRL_MANAGER.0.PID.2.I_PARAM.WR',...
     'XFEL.SYNC/LINK.LOCK/XTIN.AMC8.CONTROLLER/LSU.2.SS_PID.P_PARAM.WR',...
     'XFEL.SYNC/LINK.LOCK/XTIN.AMC8.CONTROLLER/LSU.2.SS_PID.I_PARAM.WR',...
     'XFEL.SYNC/LASER.LOCK.XLO/XHEXP1.SLO1/ADV_CTRL_MANAGER.0.PID.2.P_PARAM.WR',...
     'XFEL.SYNC/LASER.LOCK.XLO/XHEXP1.SLO1/ADV_CTRL_MANAGER.0.PID.2.I_PARAM.WR',...
     'XFEL.SYNC/LINK.LOCK/XHEXP1.AMC6.CONTROLLER/LSU.0.SS_PID.P_PARAM.WR',...
     'XFEL.SYNC/LINK.LOCK/XHEXP1.AMC6.CONTROLLER/LSU.0.SS_PID.I_PARAM.WR',...
     'XFEL.SDIAG/BAM/1932S.TL/EXPERT_STATISTICS.MACRO_PULSE.interMpStdev.1'
    };

lock_status={'XFEL.SYNC/LASER.LOCK.XLO/XTIN.MLO1/LOCK_STATUS.VALUE.RD',...
    'XFEL.SYNC/LINK.LOCK/XTIN.AMC8.CONTROLLER/LSU.2.LOCK_STATUS.VALUE.RD',...
    'XFEL.SYNC/LASER.LOCK.XLO/XHEXP1.SLO1/LOCK_STATUS.VALUE.RD',...
    'XFEL.SYNC/LINK.LOCK/XHEXP1.AMC6.CONTROLLER/LSU.0.LOCK_STATUS.VALUE.RD'};

jitter_addr={'XFEL.SYNC/LASER.LOCK.XLO/XTIN.MLO1/ADV_CTRL_MANAGER.0.PID_INPUT_JITTER.2.RD',
    'XFEL.SYNC/LINK.LOCK/XTIN.AMC8.CONTROLLER/LSU.2.TIMING_JITTER_FS.RD',
    'XFEL.SYNC/LASER.LOCK.XLO/XHEXP1.SLO1/ADV_CTRL_MANAGER.0.PID_INPUT_JITTER.2.RD',
    'XFEL.SYNC/LINK.LOCK/XHEXP1.AMC6.CONTROLLER/LSU.0.TIMING_JITTER_FS.RD'};

kphi = 1;
add = addresses;

% for i = 1:length(add)
%     c=doocsread(add{i});
%     d=c.data;
%     disp(d)
% end
%%

cond_t=[
    -1.7 -0.37;
    -0.005 -0.0002;
    0 0.8;
    0.0001 0.06;
    0.05 10;
    0.0005 0.05;
    -5 0;
    -0.16 -0.0001
    ];

cond = repmat([-1,1],size(cond_t,1),1);
scale = [1/10;1/10;1/10;1/10;1/10;1/10;1/10;1/10];


inf_ = {@infGaussLik};
mean_ = {@meanConst};
lik_ = {@likGauss};
cov_ = {@(varargin)covMaternard(3,varargin{:})};
hyp.lik = log(0.5);             %%%% must be set %%%%%
hyp.mean = 25;                  %%%% must be set %%%%%
hyp.cov = log([(cond(:,2)-cond(:,1)).*scale;10]); %%%% must be set %%%%%
acq = {@EI};


x0 = forwardCoordTransf(cond_t,[-1,-0.002,0,0.015,4.7,0.007,-0.1,-0.008]);

opts.plot=1;
opts.minFunc.mode=2;
opts.acqFunc.maxProb = 0;
opts.acqFunc.xi = 0.0;
opts.acqFunc.beta = 2;
opts.trainGP.acqVal = 10;
opts.termCondAcq = 0.05;
opts.maxIt = 35;
opts.trainGP.It = 10000;
opts.trainGP.train = 1;
opts.safeOpt = 1;
opts.newSafeOpt = 1;
opts.safeOpts.threshold = 25;       %%%% set %%%%
opts.safeOpts.thresholdOffset = 7;  %%%% set %%%%
opts.safeOpts.searchCond = 4;       %%%% set %%%%
opts.safeOpts.thresholdPer = 0.3;
opts.safeOpts.thresholdOrder = 1;
opts_lBO.maxIt = 100;
opts_lBO.sharedGP = true;
opts_lBO.subspaceDim = 1;
opts_lBO.dim_combinations = [2];
%opts_lBO.oracle = 'descent';
opts.dir_timeData="/home/jannis/master";
opts_lBO.alpha = 0.02;
use_mean = 1;

counter = 0;
save(opts.dir_timeData+'/counter_unlock.mat','counter');
time = 0;
save(opts.dir_timeData+'/time_data.mat','time');
jitter_data = [];
save(opts.dir_timeData+'/jitter_data.mat','jitter_data');

stepsize = 0.1 * (cond_t(:,2)-cond_t(:,1));
% stepsize(stepsize > 0.5) = 0.5;
fun = @(params) BAM_readWrite(params, addresses, 5, 11, lock_status, jitter_addr,stepsize,1,kphi,cond_t,opts.dir_timeData,use_mean); %%% write new setParams

start=tic;
[xopt,X,Y,DIM]=lineBO(hyp,inf_,mean_,cov_,lik_,acq, fun,cond,opts,opts_lBO,data{1,1}{16,1},data{1,2}{16,1});
toc(start);
counter_str = load(opts.dir_timeData+'/counter_unlock.mat');
jitters_str = load(opts.dir_timeData+'/jitter_data.mat');
time_str = load(opts.dir_timeData+'/time_data.mat');

data{1} = X;
data{2} = Y;
data{3} = DIM;
data{4} = jitters_str.jitter_data;
data{5} = counter_str.counter;
data{6} = time_str.time;
save(opts.dir_timeData+'/BAM_data.mat','data');
