clear all
addresses={
    'LAB.SYNC/LASER.LOCK/26A1.L3.MENHIR/ADV_CTRL_MANAGER.0.PID.1.P_PARAM.WR',...
    'LAB.SYNC/LASER.LOCK/26A1.L3.MENHIR/ADV_CTRL_MANAGER.0.PID.1.I_PARAM.WR',...
    'LAB.SYNC/LINK.LOCK/26A.AMC7.CONTROLLER/LSU.1.SS_PID.P_PARAM.WR',...
    'LAB.SYNC/LINK.LOCK/26A.AMC7.CONTROLLER/LSU.1.SS_PID.I_PARAM.WR',...
    'LAB.SYNC/LASER.LOCK/26A2.L2.ORIGAMI10/ADV_CTRL_MANAGER.0.PID.2.P_PARAM.WR',...
    'LAB.SYNC/LASER.LOCK/26A2.L2.ORIGAMI10/ADV_CTRL_MANAGER.0.PID.2.I_PARAM.WR',...
    'LAB.SYNC/LASER.LOCK/26A1.L3.MENHIR/DC8.STD_DEV.RD'
    };

% out-of-loop oxc:
% 'LAB.SYNC/LINK.LOCK/26A.AMC7.CONTROLLER/LSU.0.CTRL_IN.STD_DEV.RD' 2656
% harmonic complete jitter :
% 'LAB.SYNC/LASER.LOCK/26A1.L3.MENHIR/DC8.STD_DEV.RD' kphi = 9000;
lock_status={'LAB.SYNC/LASER.LOCK/26A2.L2.ORIGAMI10/LOCK_STATUS.VALUE.RD',...
    'LAB.SYNC/LASER.LOCK/26A1.L3.MENHIR/LOCK_STATUS.VALUE.RD',...
    'LAB.SYNC/LINK.LOCK/26A.AMC7.CONTROLLER/LSU.1.LOCK_STATUS.VALUE.RD'};

jitter_addr={'LAB.SYNC/LASER.LOCK/26A1.L3.MENHIR/CURRENT_INPUT_JITTER.RD',
    'LAB.SYNC/LINK.LOCK/26A.AMC7.CONTROLLER/LSU.1.TIMING_JITTER_FS.RD',
    'LAB.SYNC/LASER.LOCK/26A2.L2.ORIGAMI10/CURRENT_INPUT_JITTER.RD'};

cond=[
     -19.7 -0.31;
    -0.031 0;
    0 1.2;
    0.00001 0.21;
    -0.12 -0.0145;
    -0.00026 -0.000008
    ];

opts.L_f = 15;
opts.L_g = 15;
opts.start_points = 0;
opts.maxIt = 300;
opts.B = 2;
opts.alpha = 0.01;
opts.delta = 0.01;
opts.mu = 1.1;
opts.safety.threshold = 15;
g = {@(varargin) opts.safety.threshold-varargin{2}};
meas_mean = 0;
stepsize = 1 * (cond(:,2)-cond(:,1));
stepsize(stepsize > 1) = 1;

X0=[-3.39324394886329,-0.00935872887441044,0.799833496083104,0.0374160340995269,-0.0716150568183287,-8.23174052013711e-06;
    -4.48273891323462,-0.0250813654476746,0.513903591575263,0.101229812596080,-0.0720999393290359,-0.000111444113856845;
    -5.48680442245460,-0.0190768071445483,0.699583659297209,0.0528867676579610,-0.0609761364898579,-0.000104493097132908;
    -5.22527071337090,-0.00544433573264612,1.17919607966634,0.153354943888273,-0.0574760562304686,-0.000112814528013830;
    -4,-0.002,0.2,0.002,-0.04,-0.0001];

dir_ = "/home/luebsen/master/master_thesis/SMGO/Data";
obj = @(params) LAB_readWrite_SMGO(params, addresses, 60, 11, lock_status, jitter_addr,stepsize,1,9000,dir_, meas_mean);

counter = 0;
save(dir_+'/counter_unlock.mat','counter');
jitter_data = [];
save(dir_+'/jitter_data.mat','jitter_data');

%load(dir_+"/SMGO_menhir_fix_LAB.mat");
data=cell(1,7);
for i = 1:1
    %x0 = X0(i,:);
    x0 = [-11.2780800000000,-0.000310000000000000,0.200000000000000,0.00200000000000000,-0.0400000000000000,-0.000100000000000000];
    [X,Y,C,L_f,L_g] = SMGO(cond,obj,g,x0,opts);
    [yopt,I] = min(Y);
    xopt=X(I,:);
    [ybad,I2]=max(Y);
    
    data{i,1} = X;
    data{i,2} = Y;
    data{i,3} = C;
    data{i,4} = opts;
    data{i,5} = L_f;
    data{i,6} = L_g;
    load(dir_+"/counter_unlock.mat",'counter');
    data{i,7} = counter;
    load(dir_+"/jitter_data.mat",'jitter_data');
    data{i,8} = jitter_data;
    save(dir_+"/SMGO_DC_scheme_LAB.mat","data")
end
%obj(X0(end,:))
