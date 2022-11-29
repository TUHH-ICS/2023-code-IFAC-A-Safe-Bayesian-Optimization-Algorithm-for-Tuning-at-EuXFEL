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
    -0.11 -0.009;
    -0.00021 -0.00005336
    ];

cond = repmat([-1,1],size(cond_t,1),1);

scale = [1/7;1/8;1/10;1/8];

x0 = [0.051,0.03,-0.01,-0.00003];

inf_ = {@infGaussLik};
mean_ = {@meanConst};
lik_ = {@likGauss};
cov_ = {@(varargin)covMaternard(3,varargin{:})};
hyp.lik = log(0.3);
hyp.mean = 25;
hyp.cov = log([(cond(:,2)-cond(:,1)).*scale;13]);
acq = {@EI};

x0 = forwardCoordTransf(cond_t,x0);

opts.plot=0;
opts.minFunc.mode=2;
opts.maxProb = 0;
opts.acqFunc.xi = 0.0;
opts.acqFunc.beta = 2;
opts.trainGP.acqVal = 10;%0.055;%0.5 %1D       %%% 0.05 D=1 with EI; 0.5 D = 1
opts.termCondAcq = -0.5;%0.05;%0.25;%0.5 %1D    %%% 0.05 D=1 with EI; 0.25 or 0.5 D=1; 0.2 D=2 with EI sf = 5
opts.maxIt = 80;
opts.trainGP.It = 10000;
opts.trainGP.train = 0;
opts.safeOpt = 1;
opts.moSaOpt = 1;
opts.safeOpts.threshold = 25;
opts.safeOpts.thresholdOffset = 5;
opts.safeOpts.searchCond = -0.5;
opts.safeOpts.thresholdPer = 0.1;
opts.safeOpts.thresholdOrder = 1;
opts_lBO.maxIt = 50;
opts_lBO.sharedGP = true;
opts_lBO.subspaceDim = 1;
opts_lBO.m = 4;
%opts_lBO.dim_combinations = [3;4;5;6];
opts_lBO.oracle = 'random';
opts_lBO.alpha = 0.04;
opts.dir_timeData="/home/jannis/Documents/BayesOpt/bayesianoptimization/data/";
meas_mean = 1;
dir_ = "/home/jannis/Documents/BayesOpt/bayesianoptimization/data/";
kphi = 4000;

counter = 0;
save(dir_+"counter_unlock.mat",'counter');
time = 0;
save(opts.dir_timeData+"time_data.mat",'time');
jitter_data = [];
save(dir_+'jitter_data.mat','jitter_data');

stepsize = 0.2 * (cond_t(:,2)-cond_t(:,1));
% stepsize(stepsize > 1) = 1;
fun = @(params) meas(params);%LAB_readWrite(params, addresses, 10, 10.5, lock_status, jitter_addr,stepsize,1,kphi,cond_t,dir_, meas_mean);

X0=[-0.915000000000000,-0.714367350826230,1,0.825396825396825;
    -0.732137498024058,-0.938220902510097,0.178283412139096,-0.397387870827216;
    0.325616123921949,-0.338342009593391,0.196972275668599,-0.763689603106579;
    0.824264948479246,-0.241976850441243,0.491092147403435,0.472534911193277;
    -0.2, -0.05, -0.5,0.45;
    -0.515000000000000,-0.214367350826230,0.4508530805687204,0.825396825396825;
    -0.713131313131313,-0.7,0.25,0.376767676767677;
    -0.67,-0.376767676767677,0.5478,0.976767676767677;
    -0.4, -0.43,0,0;
    -0.731754134342635,-0.574796933282314,0,-0.857094374425911];

% X0=[-3.39324394886329,-0.00935872887441044,0.799833496083104,0.0374160340995269,-0.0716150568183287,-8.23174052013711e-06;
%     -4.48273891323462,-0.0250813654476746,0.513903591575263,0.101229812596080,-0.0720999393290359,-0.000111444113856845;
%     -5.48680442245460,-0.0190768071445483,0.699583659297209,0.0528867676579610,-0.0609761364898579,-0.000104493097132908;
%     -5.22527071337090,-0.00544433573264612,1.17919607966634,0.153354943888273,-0.0574760562304686,-0.000112814528013830;
%     -4,-0.002,0.2,0.002,-0.04,-0.0001];

counter = 0;
save(dir_+'/counter_unlock.mat','counter');
time = 0;
save(opts.dir_timeData+'/time_data.mat','time');
jitter_data = [];
save(dir_+'/jitter_data.mat','jitter_data');

data = load(dir_+'data.mat');
data=data.data;
for i = 4
    x0 = X0(i,:);
    [xopt,X,Y,DIM]=lineBO(hyp,inf_,mean_,cov_,lik_,acq, fun,cond,opts,opts_lBO,x0);
    data{i,1}=X;
    data{i,2}=Y;
    data{i,3}=DIM;
    jitter_data_str = load(dir_+'jitter_data.mat');
    data{i,4}=jitter_data_str.jitter_data;
    count_str=load(dir_+'counter_unlock.mat');
    counter = count_str.counter;
    data{i,5} = counter;
    time_str = load(opts.dir_timeData+'time_data.mat');
    data{i,6} = time_str.time;
    data{i,7} = xopt;
    counter = 0;
    save(dir_+'counter_unlock.mat','counter');
    time = 0;
    save(opts.dir_timeData+'time_data.mat','time');
    jitter_data = [];
    save(dir_+'jitter_data.mat','jitter_data');
    save(dir_+'data.mat','data')
end
%save("Lab_data_dim1_random.mat",'data')

%%
function y = meas(params)
    pause(1);
    y=1;
end