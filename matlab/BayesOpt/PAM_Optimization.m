clear all

addresses={'XFEL.SYNC/LASER.LOCK.XLO/XTIN.MLO1/ADV_CTRL_MANAGER.0.PID.2.P_PARAM.WR',...
     'XFEL.SYNC/LASER.LOCK.XLO/XTIN.MLO1/ADV_CTRL_MANAGER.0.PID.2.I_PARAM.WR',...
     'XFEL.SYNC/LINK.LOCK/XTIN.AMC8.CONTROLLER/LSU.2.SS_PID.P_PARAM.WR',...
     'XFEL.SYNC/LINK.LOCK/XTIN.AMC8.CONTROLLER/LSU.2.SS_PID.I_PARAM.WR',...
     'XFEL.SYNC/LASER.LOCK.XLO/XHEXP1.SLO1/ADV_CTRL_MANAGER.0.PID.2.P_PARAM.WR',...
     'XFEL.SYNC/LASER.LOCK.XLO/XHEXP1.SLO1/ADV_CTRL_MANAGER.0.PID.2.I_PARAM.WR',...
     'XFEL.SYNC/LINK.LOCK/XHEXP1.AMC5.CONTROLLER/LSU.0.SS_PID.P_PARAM.WR',...
     'XFEL.SYNC/LINK.LOCK/XHEXP1.AMC5.CONTROLLER/LSU.0.SS_PID.I_PARAM.WR',...
     'XFEL.SYNC/LASER.LOCK.EXP/XHEXP1.SASE1.PPL.OSC/ADV_CTRL_MANAGER.0.PID.2.P_PARAM.WR',...
     'XFEL.SYNC/LASER.LOCK.EXP/XHEXP1.SASE1.PPL.OSC/ADV_CTRL_MANAGER.0.PID.2.I_PARAM.WR'
    };

lock_status={'XFEL.SYNC/LASER.LOCK.XLO/XTIN.MLO1/LOCK_STATUS.VALUE.RD',...
    'XFEL.SYNC/LINK.LOCK/XTIN.AMC8.CONTROLLER/LSU.2.LOCK_STATUS.VALUE.RD',...
    'XFEL.SYNC/LASER.LOCK.XLO/XHEXP1.SLO1/LOCK_STATUS.VALUE.RD',...
    'XFEL.SYNC/LINK.LOCK/XHEXP1.AMC5.CONTROLLER/LSU.0.LOCK_STATUS.VALUE.RD',...
    'XFEL.SYNC/LASER.LOCK.EXP/XHEXP1.SASE1.PPL.OSC/LOCK_STATUS.VALUE.RD'};

jitter_addr={'XFEL.SYNC/LASER.LOCK.XLO/XTIN.MLO1/ADV_CTRL_MANAGER.0.PID_INPUT_JITTER.2.RD',...
    'XFEL.SYNC/LINK.LOCK/XTIN.AMC8.CONTROLLER/LSU.2.TIMING_JITTER_FS.RD',...
    'XFEL.SYNC/LASER.LOCK.XLO/XHEXP1.SLO1/ADV_CTRL_MANAGER.0.PID_INPUT_JITTER.2.RD',...
    'XFEL.SYNC/LINK.LOCK/XHEXP1.AMC5.CONTROLLER/LSU.0.TIMING_JITTER_FS.RD',...
    'XFEL.SYNC/LASER.LOCK.EXP/XHEXP1.SASE1.PPL.OSC/CURRENT_INPUT_JITTER.RD'};

add = addresses;

for i = 1:length(add)
    c=doocsread(add{i});
    d=c.data;
    disp(d)
end
%%

cond_t=[
    -1.7 -0.37;
    -0.005 -0.0002;
    0 0.8;
    0.0001 0.06;
    0.05 10;
    0.0005 0.05;
    -0.11 0;
    -0.025 -0.00001;
    -0.25 -0.01
    -0.0006 -0.00001 
    ];

cond = repmat([-1,1],size(cond_t,1),1);
scale = [1/10;1/10;1/10;1/10;1/10;1/10;1/10;1/10;1/10;1/10];

inf_ = {@infGaussLik};
mean_ = {@meanConst};
lik_ = {@likGauss};
cov_ = {@(varargin)covMaternard(3,varargin{:})};
hyp.lik = log(0.25);             %%%% must be set %%%%%
hyp.mean = 18;                  %%%% must be set %%%%%
hyp.cov = log([(cond(:,2)-cond(:,1)).*scale;9]); %%%% must be set %%%%%
acq = {@EI};

x0 = forwardCoordTransf(cond_t,[-1,-0.002,0.1,0.015,4.5,0.001,-0.003,-0.001,-0.15,-0.0003]);

opts.plot=1;
opts.minFunc.mode=2;
opts.acqFunc.maxProb = 0;
opts.acqFunc.xi = 0.0;
opts.acqFunc.beta = 2;
opts.trainGP.acqVal = 10;
opts.termCondAcq = 0.05;
opts.maxIt = 50;
opts.trainGP.It = 10000;
opts.trainGP.train = 1;
opts.safeOpt = 1;
opts.newSafeOpt = 1;
opts.safeOpts.threshold = 18;       %%%% set %%%%
opts.safeOpts.thresholdOffset = 5;  %%%% set %%%%
opts.safeOpts.searchCond = 3;       %%%% set %%%%
opts.safeOpts.thresholdPer = 0.2;
opts.safeOpts.thresholdOrder = 1;
opts_lBO.maxIt = 50;
opts_lBO.sharedGP = true;
opts_lBO.subspaceDim = 1;
%opts_lBO.dim_combinations = [3;4;5;6];
opts_lBO.oracle = 'descent';
opts_lBO.alpha = 0.02;

counter = 0;
save('/home/luebsen/master/master_thesis/matlab/own_lab/data/counter_unlock.mat','counter');
time = 0;
save('/home/luebsen/master/master_thesis/matlab/own_lab/data/time_data.mat','time');
jitter_data = [];
save('/home/luebsen/master/master_thesis/matlab/own_lab/data/jitter_data.mat','jitter_data');

stepsize = 0.1 * (cond_t(:,2)-cond_t(:,1));
%stepsize(stepsize > 0.5) = 0.5;

fun = @(params) readWritePAM(params, addresses, 10, 11, lock_status, jitter_addr, stepsize, 0.2, cond_t);

start=tic;
[xopt,X,Y,DIM]=lineBO(hyp,inf_,mean_,cov_,lik_,acq, fun,cond,opts,opts_lBO,x0);
toc(start);
counter_str = load('/home/luebsen/master/master_thesis/matlab/own_lab/data/counter_unlock.mat');
jitters_str = load('/home/luebsen/master/master_thesis/matlab/own_lab/data/jitter_data.mat');
time_str = load('/home/luebsen/master/master_thesis/matlab/own_lab/data/time_data.mat');

data{1} = X;
data{2} = Y;
data{3} = DIM;
data{4} = jitters_str.jitter_data;
data{5} = counter_str.counter;
data{6} = time_str.time;

