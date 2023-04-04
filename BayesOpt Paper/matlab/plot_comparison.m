clear all
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'defaultTextInterpreter','latex')
set(groot, 'defaultAxesFontsize',8)
set(groot, 'defaultTextFontsize',8)
set(groot, 'defaultLegendFontSize',8)

data1 = "LineBO_Lab_1sec/data";

data2="LineBO_descent_Lab_1sec/data";

data3="Nelder_Mead_Lab_1sec/data";

data4 = "LineBO_Lab/data";

data5="LineBO_descent_Lab/data";

data6="Nelder_Mead_Lab/data";

% data3="Nelder_Mead_Lab_1sec/data";

% data4="data_dim2_2_2.mat";

% data5="data_dim2_notSafe2_2.mat";

% data6="data_dim2_optimize2_2.mat";
Lab = 1;

x=[];
if Lab
    I_d = 2;
    if ~(data1=="")&&exist("data1",'var')
        load(data1)
        data_dim1_1 = data(:,2);
        par_dim = data(:,1);
        %unlocks(1) = sum(cell2mat(data(:,end)));
        for i = 1:size(data_dim1_1,1)
            par_temp = par_dim{i};
            data_temp = data_dim1_1{i};
            c=cellfun("isempty",data_temp);
            par_temp = par_temp(~c);
            data_temp = data_temp(~c);
            [y(1,i),I]=min(cell2mat(data_temp(end)));
            par_temp = cell2mat(par_temp(end));
            x(end+1,:) = par_temp(I,:);
        end
    end

    if ~(data2=="")&&exist("data2",'var')
        load(data2)
        data_dim1_2 = data(:,2);
        par_dim = data(:,1);
        for i = 1:size(data_dim1_2,1)
            par_temp = par_dim{i};
            data_temp = data_dim1_2{i};
            c=cellfun("isempty",data_temp);
            data_temp = data_temp(~c);
            par_temp = par_temp(~c);
            [y(2,i),I]=min(cell2mat(data_temp(end)));
            par_temp = cell2mat(par_temp(end));
            x(end+1,:) = par_temp(I,:);
        end
    end
    if ~(data3=="")&&exist("data3",'var')
        load(data3)
        data_dim1_3 = data(1:end,1);
    end
    if exist("data4",'var')&&~(data4=="")
        load(data4)
        data_dim2_1 = data(:,2);
        par_dim = data(:,1);
        for i = 1:size(data_dim2_1,1)
            par_temp = par_dim{i};
            data_temp = data_dim2_1{i};
            c=cellfun("isempty",data_temp);
            par_temp = par_temp(~c);
            data_temp = data_temp(~c);
            [y(1,i),I]=min(cell2mat(data_temp(end)));
            par_temp = cell2mat(par_temp(end));
            x(end+1,:) = par_temp(I,:);
        end
    end
    if exist("data5",'var')&&~(data5=="")
        load(data5)
        data_dim2_2 = data(:,2);
        par_dim = data(:,1);
        for i = 1:size(data_dim2_2,1)
            par_temp = par_dim{i};
            data_temp = data_dim2_2{i};
            c=cellfun("isempty",data_temp);
            data_temp = data_temp(~c);
            par_temp = par_temp(~c);
            [y(2,i),I]=min(cell2mat(data_temp(end)));
            par_temp = cell2mat(par_temp(end));
            x(end+1,:) = par_temp(I,:);
        end
    end
    if exist("data6",'var')&&~(data6=="")
        load(data6)
        data_dim2_3 = data(1:end,1);
    end
else
    I_d = 9;
    if exist("data1",'var')&&~(data1=="")
        load(data1)
        data_dim1_1 = data(:,9);
        data_x1 = data(:,8);
        yopts1=cell2mat(data(:,11));
        x(1,:) = cell2mat(data(:,2));
    end
    if exist("data2",'var')&&~(data2=="")
        load(data2)
        data_dim1_2 = data(:,9);
        data_x2 = data(:,8);
        yopts2=cell2mat(data(:,11));
        x(2,:) = cell2mat(data(:,2));
    end
    if exist("data3",'var')&&~(data3=="")
        load(data3)
        data_dim1_3 = data(:,9);
        x(3,:) = cell2mat(data(:,2));
    end
    if exist("data4",'var')&&~(data4=="")
        load(data4)
        data_dim2_1 = data(:,9);
        x(4,:) = cell2mat(data(:,2));
    end
    if exist("data5",'var')&&~(data5=="")
        load(data5)
        data_dim2_2 = data(:,9);
        x(5,:) = cell2mat(data(:,2));
    end
    if exist("data6",'var')&&~(data6=="")
        load(data6)
        data_dim2_3 = data(:,9);
        x(6,:) = cell2mat(data(:,2));
    end
end
%%

f_alpha = 0.2;
fig1=figure(1);
fig1.Units = 'centimeters';
fig1.Position(3:4)=[8.5,4];
ax=gca;
ax.Box='on';
ax.FontSize=10;

hold on
y = data_dim1_1;
[Y,std_Y]=getvals(y);
X=1:length(Y);
[yopt,xopt]=min(Y);
p1=plot(X,Y,'-','Color',[0 0.4470 0.7410],LineWidth=1.5);
p11 = plot(xopt,yopt,'*','Color',[0 0.4470 0.7410],'MarkerSize',10);
p12=fill([X,flip(X,2)],[Y+std_Y,flip(Y-std_Y,2)],[0 0.4470 0.7410],'FaceAlpha',2*f_alpha,'EdgeColor','none');


y = data_dim1_2;
[Y,std_Y]=getvals(y);
X=1:length(Y);
p2=plot(X,Y,'-','Color',[0.8500 0.3250 0.0980],LineWidth=1.5);
[yopt,xopt]=min(Y);
p21 = plot(xopt,yopt,'*','Color',[0.8500 0.3250 0.0980],'MarkerSize',10);
p22 = fill([X,flip(X,2)],[Y+std_Y,flip(Y-std_Y,2)],[0.8500 0.3250 0.0980],'EdgeColor','none');
p22.FaceAlpha = 2*f_alpha;

y = data_dim1_3;
[Y,std_Y]=getvals(y);
X=1:length(Y);
p3=plot(X,Y,'-','Color','k'	,LineWidth=1.5);
[yopt,xopt]=min(Y);
p31 = plot(xopt,yopt,'*','Color','k'	,'MarkerSize',10);
p32 = fill([X,flip(X,2)],[Y+std_Y,flip(Y-std_Y,2)],[0 0 0],'EdgeColor','none');
p32.FaceAlpha = f_alpha;
grid on
xlim([0 130])
ax = gca;
set(gca,'TickLabelInterpreter','latex');
xlabel("Iteration $n$","Interpreter","latex")
ylabel("$$J_{opt}(n)$$ [fs]",'Interpreter','latex')


fig=figure(2);
fig.Units = 'centimeters';
fig.Position(3:4)=[8.5,4];
ax=gca;
ax.Box='on';
ax.FontSize=10;
% fig.Position(1)=fig1.Position(1)+fig1.Position(3)+0.8
hold on
y = data_dim2_1;
[Y,std_Y]=getvals(y);
X=1:length(Y);
p4=plot(X,Y,'-','Color',[0 0.4470 0.7410]	,LineWidth=1.5);
[yopt,xopt]=min(Y);
p41 = plot(xopt,yopt,'*','Color',[0 0.4470 0.7410]	,'MarkerSize',10);
p42 = fill([X,flip(X,2)],[Y+std_Y,flip(Y-std_Y,2)],[0 0.4470 0.7410],'EdgeColor','none');
p42.FaceAlpha = 2*f_alpha;

y = data_dim2_2;
[Y,std_Y]=getvals(y);
X=1:length(Y);
p5=plot(X,Y,'-','Color',[0.8500 0.3250 0.0980]	,LineWidth=1.5);
[yopt,xopt]=min(Y);
p51 = plot(xopt,yopt,'*','Color',[0.8500 0.3250 0.0980]	,'MarkerSize',10);
p52=fill([X,flip(X,2)],[Y+std_Y,flip(Y-std_Y,2)],[0.8500 0.3250 0.0980],'FaceAlpha',2*f_alpha,'EdgeColor','none');
y = data_dim2_3;
[Y,std_Y]=getvals(y);
X=1:length(Y);
p6=plot(X,Y,'-','Color','k'	,LineWidth=1.5);
[yopt,xopt]=min(Y);
p61 = plot(xopt,yopt,'*','Color','k'	,'MarkerSize',10);
p62 = fill([X,flip(X,2)],[Y+std_Y,flip(Y-std_Y,2)],'k','EdgeColor','none');
p62.FaceAlpha = f_alpha;
ax=gca;
ax.Box='on';
grid on
xlim([0 130])
% 
if Lab
    l1=legend([p1,p2,p3],'\texttt{LineBO} random','\texttt{LineBO} descent', 'Nelder-Mead','Interpreter','latex','NumColumns',2,'Location','Northeast','fontsize',8)
    l2=legend([p4,p5,p6],'\texttt{LineBO} random','\texttt{LineBO} descent', 'Nelder-Mead','Interpreter','latex','NumColumns',2,'Location','Northeast','fontsize',8)
    ax = gca;
    set(gca,'TickLabelInterpreter','latex');
    xlabel("Iteration $n$","Interpreter","latex")
    ylabel("$$J_{opt}(n)$$ [fs]",'Interpreter','latex')
    l1.Position(1)=l1.Position(1)+0.02;
    l2.Position(1)=l2.Position(1)+0.02;
%     ax.Box='on';
else
    ylim([10,30])
    hold off
    grid on
    xlim([0,1000])
    ylim([8,35])
    ax = gca;
    set(gca,'TickLabelInterpreter','latex');
    xlabel("iteration $n$","Interpreter","latex")
    ylabel("$$J_{opt}(n)$$ [fs]",'Interpreter','latex')
    l1=legend([p1,p2],'SafeOpt','MoSaOpt','Interpreter','latex','NumColumns',1,'Location','Northeast','fontsize',8);
    title(l1,'LineBO +','Interpreter','latex')
    a=axes('Position',get(ax,'position'),'visible','off');
    ax.Box = 'on';
    l2 = legend(a,[p4,p6],'SafeOpt','MoSaOpt','Interpreter','latex','NumColumns',1,'fontsize',8);
    title(l2,'PlaneBO + ','Interpreter','latex')
    l2.Units='centimeters';
    l1.Units='centimeters';
    l2.Position(1:2)=[l1.Position(1)-3,l1.Position(2)];
end
%%
y1opt = mean(yopts1);
y2opt = mean(yopts2);
x= categorical({'$\kappa = 0.9$','$\kappa = 0$'})
bar([1;2;3],[12.18;y1opt;y2opt]);
hold on
errorbar([2,3],[y1opt,y2opt],[min(yopts1)-y1opt,min(yopts2)-y2opt],[max(yopts1)-y1opt,max(yopts2)-y2opt],'s')
ax = gca;
ax.XTickLabel={'noise free','$\kappa = 0.9$','$\kappa = 0$'};


%% Configuration
N = 3;  % Number of repetitions of G in system chain
N_L = 2;
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

Kp_max = 3e1;
Kp_min = 0.2;
Ki_max = 6e1;
Ki_min = 0;

cond_t=[Kp_min, Kp_max;
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

scale=[1/sys.k_phi 1/sys_link.k_phi];
%%
for i=1:10
    temp = data_x1{i};
    temp=temp(~cellfun('isempty',temp));
    temp = temp{end};
    len = length(temp);
    y=zeros(len,1);
    for j = 1 : len
        y(j) = connect_PI(temp(j,:),Gg,scale,cond_t);
    end
    data{i,12} = y;
end
%%
fig =figure(5);
hold on
f_alpha = 0.2;
y = data_dim1_1;
[Y,std_Y]=getvals(y)
X=1:length(Y);
[yopt,xopt]=min(Y);
p1=plot(X,Y,'-','Color',[0 0.4470 0.7410],LineWidth=1.2)
p11 = plot(xopt,yopt,'*','Color',[0 0.4470 0.7410],'MarkerSize',10)
fill([X,flip(X,2)],[Y+std_Y,flip(Y-std_Y,2)],[0 0.4470 0.7410],'FaceAlpha',2*f_alpha,'EdgeColor','none');


y = data_dim1_2(:,1);
[Y,std_Y]=getvals(y);
% Y_true = cellfun(@(x)connect_PI(x,Gg,scale,cond_t),X_d);
X=1:length(Y);
p2=plot(X,Y,'-','Color',[0.8500 0.3250 0.0980],LineWidth=1)
[yopt,xopt]=min(Y);
p21 = plot(xopt,yopt,'*','Color',[0.8500 0.3250 0.0980],'MarkerSize',10)
p22 = fill([X,flip(X,2)],[Y+std_Y,flip(Y-std_Y,2)],[0.8500 0.3250 0.0980],'EdgeColor','none');
p22.FaceAlpha = 2*f_alpha;
hold off
%%


function varargout = getvals(y,varargin)
    yt = find_cells(y);
    len=max(cellfun('length',yt));
    yt1=zeros(length(yt),len);
    if ~isempty(varargin)
        x=varargin{1};
        xt=find_cells(x);
        D=size(xt{1},2);
        xt1 = cell(length(yt),len);
    end
    for i = 1:length(yt)
        temp = ones(len,1)*min(yt{i});
        temp(1:length(yt{i})) = yt{i};
        if ~isempty(varargin)
            tempx = xt{i};
        end
        for j=1:len
            [yt1(i,j),I]=min(temp(1:j));
            if ~isempty(varargin)
                xt1{i,j} = tempx(I,:); 
            end
        end
    end
    varargout{1} = mean(yt1,1);
    varargout{2} = std(yt1,0);
    if ~isempty(varargin)
        varargout{3} = xt1;
    end
end

function yt=find_cells(y)
    yt=cell(1,length(y));
    for i = 1:length(y)
        temp = y{i};
        if iscell(temp)
            temp = temp(~cellfun('isempty',temp));
            temp = temp{end};
            yt{i}=temp;
        else
            yt{i} = temp;
        end
    end
end

function [y] = connect_PI(pi_params, Gg, scale,cond)
    pi_params=backwardCoordTransf(cond,pi_params);
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
