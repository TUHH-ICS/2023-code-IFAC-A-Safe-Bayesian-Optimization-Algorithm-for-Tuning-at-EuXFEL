%% Example script to optimize test plant at the desy laboratory

clear all
addresses={...'LAB.SYNC/LASER.LOCK/26A1.L3.MENHIR/ADV_CTRL_MANAGER.0.PID.0.P_PARAM.WR',...
    ...'LAB.SYNC/LASER.LOCK/26A1.L3.MENHIR/ADV_CTRL_MANAGER.0.PID.0.I_PARAM.WR',...
    'LAB.SYNC/LINK.LOCK/26A.AMC7.CONTROLLER/LSU.1.SS_PID.P_PARAM.WR',...
    'LAB.SYNC/LINK.LOCK/26A.AMC7.CONTROLLER/LSU.1.SS_PID.I_PARAM.WR',...
    'LAB.SYNC/LASER.LOCK/26A1.L1.ORIGAMI15/ADV_CTRL_MANAGER.0.PID.2.P_PARAM.WR',...
    'LAB.SYNC/LASER.LOCK/26A1.L1.ORIGAMI15/ADV_CTRL_MANAGER.0.PID.2.I_PARAM.WR',...
...     'LAB.SYNC/LINK.LOCK/26A.AMC7.CONTROLLER/LSU.0.CTRL_IN.STD_DEV.RD'
    'LAB.SYNC/LASER.LOCK/26A1.L1.ORIGAMI15/DCS_7.SPEC',...
    'LAB.SYNC/LASER.LOCK/26A1.L1.ORIGAMI15/DCS_8.SPEC'
    };

% out-of-loop oxc: kphi = 4000
% out-of-loop oxc:
% 'LAB.SYNC/LINK.LOCK/26A.AMC7.CONTROLLER/LSU.0.CTRL_IN.STD_DEV.RD' 2656
% harmonic complete jitter :
% 'LAB.SYNC/LASER.LOCK/26A1.L3.MENHIR/DC8.STD_DEV.RD' kphi = 9000;
lock_status={'LAB.SYNC/LASER.LOCK/26A1.L1.ORIGAMI15/LOCK_STATUS.VALUE.RD',...
    'LAB.SYNC/LASER.LOCK/26A1.L3.MENHIR/LOCK_STATUS.VALUE.RD',...
    'LAB.SYNC/LINK.LOCK/26A.AMC7.CONTROLLER/LSU.1.LOCK_STATUS.VALUE.RD'};

jitter_addr={'LAB.SYNC/LASER.LOCK/26A1.L3.MENHIR/CURRENT_INPUT_JITTER.RD',...
    'LAB.SYNC/LINK.LOCK/26A.AMC7.CONTROLLER/LSU.1.TIMING_JITTER_FS.RD',...
    'LAB.SYNC/LASER.LOCK/26A1.L1.ORIGAMI15/CURRENT_INPUT_JITTER.RD'};

%%% MENHIR included %%%
% cond_t=[
%     -19.7 -0.31;
%     -0.031 0;
%     0 1.2;
%     0.00001 0.21;
%     -0.12 -0.0145;
%     -0.00026 -0.000008
%     ];

%%% MENHIR not included %%%
cond_t=[
    0 1.2;
    0.00001 0.21;
    -0.12 -0.0145;
    -0.00026 -0.000008
    ];

cond = repmat([-1,1],size(cond_t,1),1);

scale = [1/5;1/8;1/7;1/5];

x0 = [0.051,0.03,-0.01,-0.00003];
% x0 = [-4,-0.002,0.2,0.002,-0.04,-0.0001];

inf_ = {@infGaussLik};
mean_ = {@meanConst};
lik_ = {@likGauss};
cov_ = {@(varargin)covMaternard(3,varargin{:})};
hyp.lik = log(0.1);
hyp.mean = 20;
hyp.cov = log([(cond(:,2)-cond(:,1)).*scale;7]);
acq = {@EI};

x0 = forwardCoordTransf(cond_t,x0);

opts.plot=1;
opts.minFunc.mode=2;
opts.maxProb = 0;
opts.acqFunc.xi = 0.0;
opts.acqFunc.beta = 2;
opts.trainGP.acqVal = 10;%0.055;%0.5 %1D       %%% 0.05 D=1 with EI; 0.5 D = 1
opts.termCondAcq = 0.05;%0.05;%0.25;%0.5 %1D    %%% 0.05 D=1 with EI; 0.25 or 0.5 D=1; 0.2 D=2 with EI sf = 5
opts.maxIt = 50;
opts.trainGP.It = 10000;
opts.trainGP.train = 1;
opts.safeOpt = 1;
opts.moSaOpt = 1;
opts.safeOpts.threshold = 20;
opts.safeOpts.thresholdOffset = 8;
opts.safeOpts.searchCond = 2;
opts.safeOpts.thresholdPer = 0.1;
opts.safeOpts.thresholdOrder = 1;
opts_lBO.maxIt = 100;
opts_lBO.sharedGP = true;
opts_lBO.subspaceDim = 1;
%opts_lBO.dim_combinations = [3;4;5;6];
opts_lBO.oracle = 'random';
opts_lBO.alpha = 0.02;
opts.dir_timeData="/home/luebsen/bayesopt/bayesianoptimization/data/";
meas_mean = 1;
dir_ = "/home/luebsen/bayesopt/bayesianoptimization/data/";
kphi = 4000;

counter = 0;
save(dir_+"counter_unlock.mat",'counter');
time = 0;
save(opts.dir_timeData+"time_data.mat",'time');
jitter_data = [];
save(dir_+'jitter_data.mat','jitter_data');

stepsize = 0.2 * (cond_t(:,2)-cond_t(:,1));
% stepsize(stepsize > 1) = 1;
fun = @(params) LAB_readWrite(params, addresses, 5, 10.5, lock_status, jitter_addr,stepsize,1,kphi,cond_t,dir_, meas_mean);

tic
[xopt,X,Y,DIM]=lineBO(hyp,inf_,mean_,cov_,lik_,acq, fun,cond,opts,opts_lBO,x0);
toc
    