clear all; close all;
%% Configuration
N = 3;  % Number of repetitions of G in system chain
N_L = 2;
%% Build model from parameters
lbsync = load('lbsync.mat');  % Load parameter database
tunit = 'seconds'; tunitexp = 0;  % Overwrite time scale (to evaluate effect on optimization problem)
scaling = 15;  % Exponent of the model output scaling, 0 = s, 12 = ps, 15 = fs, etc.
sys = build_laser_model(lbsync.sim.laser.origami, scaling, tunit);

% Plant
ctrl_gain = 1;
G = balreal(ss(series(sys.G_pzt, sys.G_l) / ctrl_gain));
G.u = 'u';
G.y = 'phi';
% Reference noise coloring filter
Fr = sys.Fr;
%Fr.P{1}(1) = -1e1 * 2*pi * 10^(-tunitexp);  % Change integral behaviour to frequency region of interest
Fr = balreal(ss(Fr));  % Alt: ss, balreal, prescale
Fr.D = zeros(size(Fr.D));  % Make proper

% Plant output disturbance coloring filter
Fd = sys.Fd;
%Fd.P{1}(1) = -1e1 * 2*pi * 10^(-tunitexp);  % Change integral behaviour to frequency region of interest
Fd = balreal(ss(Fd));  % Alt: ss, balreal, prescale
Fd.D = zeros(size(Fd.D));  % Make proper

Glaser = connect(G,Fd,sumblk('y = phi + d'),{'w','u'},{'y'});
Glaser2 = sys.G;

% Link model
if N_L > 0
    
    sys_link = lbsync.sim.link.short;
    sys_link.Fd.P{1}(end) = -1e-1 * 2*pi;
    G_pz = zpk(sys_link.G_pz);
    G_pz.P{1} = G_pz.P{1}(1:2);
    sys_link.G_pz = ss(G_pz);
    clear G_pz;
    
    sys_link = build_link_model(sys_link, scaling, tunit);
    Glink = sys_link.Gpade;
    %Glink = repmat({Glink}, 1, N_L);
end

% Connectivity
Fr.u = 'w(1)';
Fr.y = 'r';

if N_L < 1
    G = repmat({Glaser}, 1, N);
    sums = cell(1, N);
   
    for i = 1:N
        G{i}.u = {sprintf('w(%d)', i+1);sprintf('u(%d)', i)};
        G{i}.y = sprintf('y(%d)', i);
        if i == 1
            sums{i} = sumblk('e(1) = r - y(1)');
        else
            sums{i} = sumblk(sprintf('e(%1$d) = y(%2$d) - y(%1$d)', i, i-1));
        end
    end
    sums{end+1} = sumblk(sprintf('z = r - y(%d)', N));
else
    G = cell(1, N + N_L);
    G{1} = Glaser;
    for i=2:2:N+N_L
        G{i} = Glink;
        G{i+1} = Glaser;
    end
    for i = 1:N+N_L
        if mod(i,2) ~= 0
            G{i}.u = {sprintf('w(%d)', i+1) ;sprintf('u(%d)', i)};
            G{i}.y = sprintf('y(%d)', i);
            if i == 1
                sums{i} = sumblk('e(1) = r - y(1)');
            else
                sums{i} = sumblk(sprintf('e(%1$d) = y(%2$d) - y(%1$d)', i, i-1));
            end
        else
            G{i}.u = {sprintf('y(%d)',i-1); sprintf('w(%d)', i+1); sprintf('u(%d)', i)};
            G{i}.y = {sprintf('l(%d)', i/2);sprintf('y(%d)', i)};
            sums{i} = sumblk(sprintf('e(%d) = y(%d) - l(%d)', i, i-1, i/2));
        end
    end
    sums{end+1} = sumblk(sprintf('z = r - y(%d)', N+N_L));
end
Gg = connect(G{:}, Fr, sums{:}, {'u','w'},{'e','z'});
%%
Kp_max = 3e1;
Kp_min = 0.2;
Ki_max = 6e1;
Ki_min = 0;

cond=[Kp_min, Kp_max;
      Ki_min, Ki_max;
      %1.2, 3.8;
      0, 0.000105*350;
      0, 3;
      Kp_min, Kp_max;
      Ki_min, Ki_max;
      %1.2, 3.8;
      0, 0.000105*350;
      0, 3;
      Kp_min, Kp_max;
      Ki_min, Ki_max];

x0 = [23.71978, 13.51946, 0.01, 2.5, 9.57467, 26.47217, 0.02, 3, 3.19737, 12.91408];
X0 = [20.6963   21.5537    0.0271    1.1841   20.5658   42.2428    0.0163    0.0587   10.0596   25.4586
   11.4464   12.9611    0.0290    2.8479    9.9614   40.2759    0.0161    2.5005   23.1119   10.0352
   23.7409   19.1115    0.0196    0.2699    3.5288    8.1776    0.0249    1.4855    5.8534   29.7003
   17.3249   50.7107    0.0271    1.7580    7.5527   39.9850    0.0031    1.8779   19.8961   43.7851
   25.3767   12.5643    0.0203    1.8897    1.1533   36.8828    0.0133    0.1486   14.7892   11.5506
    3.8679   12.3297    0.0054    0.5672    1.4710   38.1119    0.0104    1.6158   20.9159   29.9470
   12.6275   12.3585    0.0348    0.2462    3.3501    8.5225    0.0061    1.8629   17.2966    3.1247
    5.4926   23.9154    0.0049    0.0927   28.1864   18.0784    0.0109    0.9988   14.1186   38.8919
    0.9518   50.5324    0.0205    2.5623   10.5668   26.7616    0.0020    0.5313   19.9517   19.8497
   24.5825    6.0133    0.0065    1.0789    1.8898   31.3131    0.0123    0.5270    6.4266   54.3092];
% X0 = x0;
opts.L_f = 15;
opts.L_g = 15;
opts.start_points = 0;
opts.maxIt = 750;
opts.B = 2;
opts.alpha = 0.2;
opts.delta = 0.2;
opts.mu = 1.1;
opts.safety.threshold = 50;
g = {@(varargin) opts.safety.threshold-varargin{2}};
obj = @(params) connect_PI(params, Gg, [1/sys.k_phi 1/sys_link.k_phi]);
data = cell(size(X0,1),4);
for i = 1:size(X0)
    x0 = X0(i,:);
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
    save("new_algorithm/Data/SMGO_test_different_inits2.mat","data")
end


%%

function [y] = connect_PI(pi_params, Gg, scale)
    N = length(pi_params)/2;
    C = cell(1,N);
    len_scale = length(scale);
    for i=1:N
        C{i} = pid(pi_params(1,2*i-1)*scale(len_scale-mod(i,len_scale)),pi_params(1,2*i)*scale(len_scale-mod(i,len_scale)));
        C{i}.y = sprintf('u(%d)', i);
        C{i}.u = sprintf('e(%d)', i);
    end
    Gcl = connect(Gg,C{:}, 'w','z');
    y = norm(Gcl,2);
    if y == inf
        y = 1e4;
    end
end
