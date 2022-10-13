addresses={'LAB.SYNC/LASER.LOCK/26A1.L3.MENHIR/ADV_CTRL_MANAGER.0.PID.1.P_PARAM.WR',...
    'LAB.SYNC/LASER.LOCK/26A1.L3.MENHIR/ADV_CTRL_MANAGER.0.PID.1.I_PARAM.WR',...
    'LAB.SYNC/LINK.LOCK/26A.AMC7.CONTROLLER/LSU.1.SS_PID.P_PARAM.WR',...
    'LAB.SYNC/LINK.LOCK/26A.AMC7.CONTROLLER/LSU.1.SS_PID.I_PARAM.WR',...
    'LAB.SYNC/LASER.LOCK/26A2.L2.ORIGAMI10/ADV_CTRL_MANAGER.0.PID.2.P_PARAM.WR',...
    'LAB.SYNC/LASER.LOCK/26A2.L2.ORIGAMI10/ADV_CTRL_MANAGER.0.PID.2.I_PARAM.WR',...
    'LAB.SYNC/LASER.LOCK/26A1.L3.MENHIR/DC8.STD_DEV.RD'};

lock_status={'LAB.SYNC/LASER.LOCK/26A2.L2.ORIGAMI10/LOCK_STATUS.VALUE.RD',...
    'LAB.SYNC/LASER.LOCK/26A1.L3.MENHIR/LOCK_STATUS.VALUE.RD',...
    'LAB.SYNC/LINK.LOCK/26A.AMC7.CONTROLLER/LSU.1.LOCK_STATUS.VALUE.RD'};

jitter_addr={'LAB.SYNC/LASER.LOCK/26A1.L3.MENHIR/CURRENT_INPUT_JITTER.RD',
    'LAB.SYNC/LINK.LOCK/26A.AMC7.CONTROLLER/LSU.1.TIMING_JITTER_FS.RD',
    'LAB.SYNC/LASER.LOCK/26A2.L2.ORIGAMI10/CURRENT_INPUT_JITTER.RD'};

cond=[
    -20 -0.31;
    -0.031 0;
    0 1.2;
    0.00001 0.21;
    -0.12 -0.0145;
    -0.00026 -0.000008
    ];
stepsize = 1 * (cond(:,2)-cond(:,1));
stepsize(stepsize > 0.5) = 0.5;
fun = @(params) setParams_NelderMead(params, addresses, 60, 11, lock_status, jitter_addr,stepsize,1,cond,9000);

x0=[-4,-0.002,0.2,0.002,-0.04,-0.0001];

options = optimset('FunValCheck','on');
options.Display ='iter';
options.PlotFcns = {@optimplotfval,@optimplotfunccount};
options.MaxFunEvals = 250;
options.TolFun = 0.1;
options.TolX = 0.1;

X0=[-3.39324394886329,-0.00935872887441044,0.799833496083104,0.0374160340995269,-0.0716150568183287,-8.23174052013711e-06;
    -4.48273891323462,-0.0250813654476746,0.513903591575263,0.101229812596080,-0.0720999393290359,-0.000111444113856845;
    -5.48680442245460,-0.0190768071445483,0.699583659297209,0.0528867676579610,-0.0609761364898579,-0.000104493097132908;
    -5.22527071337090,-0.00544433573264612,1.17919607966634,0.153354943888273,-0.0574760562304686,-0.000112814528013830;
    -4,-0.002,0.2,0.002,-0.04,-0.0001];

data = cell(5,5);
counter = 0;
save('/home/luebsen/master/master_thesis/matlab/own_lab/test/counter_unlock.mat','counter');
for i = 5:size(X0,1)
    x0 = sinTransform(X0(i,:),cond);
    start=tic;
    [x]=fminsearch(fun,x0,options);
    t = toc(start);
    val_str=load('/home/luebsen/master/master_thesis/matlab/own_lab/test/val.mat');
    val = val_str.val;
    data{i,1} = val(:,1);
    data{i,3} = t;
    data{i,2} = val(:,2);
    data{i,4} = sinBackTrans(x,cond);
    count_str=load('/home/luebsen/master/master_thesis/matlab/own_lab/test/counter_unlock.mat');
    counter = count_str.counter;
    data{i,5} = counter;
    counter = 0;
    save('/home/luebsen/master/master_thesis/matlab/own_lab/test/counter_unlock.mat','counter');
    val = [];
    save('/home/luebsen/master/master_thesis/matlab/own_lab/test/val.mat','val');
end
%% Sine transformation to restrict nelder-mead to the 'cond' intervall

function y=sinTransform(x,cond)
    D=size(cond,1);
    y = zeros(size(x));
    for i = 1:D
        if x(i) <= cond(i,1)
            y(i) = -pi/2;
        elseif x(i) >= cond(i,2)
            y(i) = pi/2;
        else
            y(i) = 2*(x(i)-cond(i,1))/(cond(i,2)-cond(i,1))-1;
            y(i) = 2*pi+asin(max(-1,min(1,y(i))));
        end
    end
end

function x=sinBackTrans(y,cond)
    x=(cond(:,2)'-cond(:,1)').*(sin(y)+1)/2 + cond(:,1)';
end
