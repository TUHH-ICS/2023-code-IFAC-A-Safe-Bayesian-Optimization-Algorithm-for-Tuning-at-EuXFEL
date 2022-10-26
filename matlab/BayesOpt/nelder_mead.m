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
cond=[
    0 1.2;
    0.00001 0.21;
    -0.11 -0.009;
    -0.00021 -0.00005336
    ];

stepsize = 1 * (cond(:,2)-cond(:,1));
% stepsize(stepsize > 0.5) = 0.5;
fun = @(params) setParams_NelderMead(params, addresses, 1, 11, lock_status, jitter_addr,stepsize,1,cond,4000);


options = optimset('FunValCheck','on');
options.Display ='iter';
options.PlotFcns = {@optimplotfval,@optimplotfunccount};
options.MaxFunEvals = 200;
options.TolFun = 0.01;
options.TolX = 0.001;

    X0=[0.0509999999999999,0.0300000000000000,-0.00900000000000000,-6.70349206349207e-05;
        0.160717501185565,0.00649649634095237,-0.0504966876869757,-0.000162803418043188;
        0.795369674353169,0.0694807807027419,-0.0495529000787358,-0.000191492169715307;
        1.09455896908755,0.0795986405879217,-0.0346998465561265,-9.46710657553426e-05;
        0.480000000000000,0.0997552500000000,-0.0847500000000000,-9.64360000000000e-05;
        0.291000000000000,0.0824975000000000,-0.0367319194312796,-6.70349206349207e-05;
        0.172121212121212,0.0315085000000000,-0.0468750000000000,-0.000102171555555556;
        0.198000000000000,0.0654462777777778,-0.0318361000000000,-5.51795555555555e-05;
        0.360000000000000,0.0598571500000000,-0.0595000000000000,-0.000131680000000000;
        0.160947519394419,0.0446541959900235,-0.0595000000000000,-0.000198807631405037];
% val = 0;
% save('/home/luebsen/bayesopt/bayesianoptimization/data/val.mat','val');
% data = cell(10,5);
% counter = 0;
% save('/home/luebsen/bayesopt/bayesianoptimization/data/counter_unlock.mat','counter');
for i = 1:10
%     x0 = sinTransform(X0(i,:),cond);
    x0 = X0(i,:);
    start=tic;
    [x]=fminsearch(fun,x0,options);
    t = toc(start);
    val_str=load('/home/luebsen/bayesopt/bayesianoptimization/data/val.mat');
    val = val_str.val;
    data{i,1} = val(:,1);
    data{i,3} = t;
    data{i,2} = val(:,2);
    data{i,4} = sinBackTrans(x,cond);
    count_str=load('/home/luebsen/bayesopt/bayesianoptimization/data/counter_unlock.mat');
    counter = count_str.counter;
    data{i,5} = counter;
    counter = 0;
    save('/home/luebsen/bayesopt/bayesianoptimization/data/counter_unlock.mat','counter');
    val = [];
    save('/home/luebsen/bayesopt/bayesianoptimization/data/val.mat','val');
    save('/home/luebsen/bayesopt/bayesianoptimization/data/data.mat','data');
end
%%
fun = @(params) setParams_NelderMead(params, addresses, 10, 11, lock_status, jitter_addr,stepsize,1,cond,4000);
for i = 1:10
    x=sinTransform(data{i,4},cond);
    data{i,6}=fun(x);
end
    %%
sinTransform(data{10,4},cond) 

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
%             y(i) = 2*(x(i)-cond(i,1))/(cond(i,2)-cond(i,1))-1;
%             y(i) = 2*pi+asin(max(-1,min(1,y(i))));
              y(i) = asin(2*(x(i)-cond(i,1))/(cond(i,2)-cond(i,1))-1);
        end
    end
end

function x=sinBackTrans(y,cond)
    x=(cond(:,2)'-cond(:,1)').*(sin(y)+1)/2 + cond(:,1)';
end
