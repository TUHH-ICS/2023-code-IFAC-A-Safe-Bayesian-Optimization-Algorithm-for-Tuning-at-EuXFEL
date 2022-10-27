clear all
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'defaultTextInterpreter','latex')
set(groot, 'defaultAxesFontsize',10)
set(groot, 'defaultTextFontsize',10)
set(groot, 'defaultLegendFontSize',10)

x = zeros(6,10);
load("data_dim1_2.mat")
data_dim1_1 = data;
x(1,:) = cell2mat(data_dim1_1(:,2));
load("data_dim1_notSafe2.mat")
data_dim1_2 = data;
x(2,:) = cell2mat(data_dim1_2(:,2));
load("data_dim1_optimize2_2.mat")
data_dim1_3 = data;
x(3,:) = cell2mat(data_dim1_3(:,2));
load("data_dim2_2_2.mat")
data_dim2_1 = data;
x(4,:) = cell2mat(data_dim2_1(:,2));
load("data_dim2_notSafe2_2.mat")
data_dim2_2 = data;
x(5,:) = cell2mat(data_dim2_2(:,2));
load("data_dim2_optimize2_2.mat")
data_dim2_3 = data;
x(6,:) = cell2mat(data_dim2_3(:,2));

y = x;
x=1:10;
figure(3)
plot(x,y)
legend("D = 1; safe","D = 1; naive", "D = 1; safe \& optimized","D = 2; safe","D = 2; naive", "D = 2; safe \& optimized",'interpreter','latex','NumColumns',2,'Location','northoutside')
ylabel("$$\mbox{J}_{opt}$$",'Interpreter',"latex")
xlabel("start point",'Interpreter','latex')
xlim([x(1), x(end)])
set(gca,'TickLabelInterpreter','latex');

fig=figure(2)
%fig.Position = [500,500, 1500,500];
yt=mean(y,2)';
std_y = std(y,[],2);
std_y=std_y([1,2,3;4,5,6]);
y=[yt(1:3);yt(4:end)];
%x=categorical({'D = 1','D = 2',});
hold on
b=bar(y,'grouped');

xtips=[b(1).XEndPoints;b(2).XEndPoints;b(3).XEndPoints]';
ytips=[b(1).YEndPoints,b(2).YEndPoints,b(3).YEndPoints];
labels=cell(length(ytips),1);
for i = 1:length(yt)
    labels{i}=sprintf("$$%.2f \\pm %.3f$$",ytips(i),std_y(i));
end
%text(xtips(:),ytips,labels,'HorizontalAlignment','center','VerticalAlignment','bottom','Interpreter','latex','FontSize',11)
% ngroups = 3;
% nbars = 2;
% x = nan(ngroups, nbars);
% for i = 1:ngroups
%     x(i,:) = b(i).XEndPoints;
% end
% Plot the errorbars
errorbar(xtips,y,2*std_y,'k','linestyle','none');
yt2=round(yt,2);
yt2=num2cell(yt2);
yt2=cellfun(@num2str,yt2,'UniformOutput',false);
text(reshape(xtips',[1,6]),ones(1,6)*13.5,yt2,'HorizontalAlignment','center','FontSize',12)
ax = gca;
ax.Box = 'on';
ax.XTick = [1, 2];
ax.XTickLabel = {'$$\dim(\mathcal{L}) = 1$$', '$$\dim(\mathcal{L}) = 2$$'};
grid on
hold off
set(gca,'TickLabelInterpreter','latex');
%hold off
ylim([0 17])
ylabel("$$J_{opt}$$ [fs]",'Interpreter','latex')
legend('SafeOpt','Naive','MoSaOpt','interpreter','latex','numColumns',3,'Location','north')

%%
y=[];
y(1,:) = cell2mat(data_dim1_1(:,3));
y(2,:) = cell2mat(data_dim1_2(:,3));
y(end+1,:) = cell2mat(data_dim1_3(:,3));
y(end+1,:) = cell2mat(data_dim2_1(:,3));
y(5,:) = cell2mat(data_dim2_2(:,3));
y(end+1,:) = cell2mat(data_dim2_3(:,3));
l=size(y,1);

yt=max(y,[],2)';
y=[yt(1:l/2);yt(l/2+1:end)];
%x=categorical({'D = 1','D = 2',});
fig=figure(1);
b1=bar(y, 'grouped');
%b1(2).FaceColor = [0.93,0.69,0.13];
hold on
y1=yline(50,'--','threshold', 'LineWidth',2,'Fontsize',12);
yl.LabelHorizontalAlignment = 'left';
ax = gca;
ax.XTick = [1, 2];
ax.XTickLabel = {'$$\dim(\mathcal{L}) = 1$$', '$$\dim(\mathcal{L}) = 2$$'};
hold off
ylim([30 60])
ylabel("$$J_{max}$$ [fs]",'Interpreter','latex')
legend("SafeOpt","Naive","MoSaOpt",'interpreter','latex', 'Location','north','NumColumns',1)
grid on 
% y=[];
% y(1,:) = cell2mat(data_dim1_2(:,3));
% y(2,:) = cell2mat(data_dim2_2(:,3));
% y=max(y,[],2);
% subplot(1,2,2)
% b2=bar([1;2],y, 'grouped')
% b2.FaceColor=[0.85,0.33,0.10];
% hold on
% y1=yline(50,'--','threshold', 'LineWidth',2,'FontSize',12);
% yl.LabelHorizontalAlignment = 'left';
% ax = gca;
% %ax.XTick = [1, 2];
% ax.XTickLabel = {'$$\dim(\mathcal{L}) = 1$$', '$$\dim(\mathcal{L}) = 2$$'};
% hold off
% ylabel("$$J_{max}$$",'Interpreter','latex')
% legend("Naive",'interpreter','latex', 'Location','northwest','NumColumns',2,'Fontsize',15)
%%
y=[];
load("data_dim1_optimize2_2.mat")
data_dim1_1 = data;
load("Lab_data_dim1_descent2.mat")
data_dim1_2 = data;
load("data_nelder_mead.mat")
data_dim1_3 = data;
t = data(:,1);
for i = 1:length(t)
    temp=t{i};
    t{i}=temp(1:400);
end
data_dim1_3(:,1)=t;

% load("SMGO_test_different_inits2.mat")
% data_dim1_3 = data;
% load("SMGO_test_different_inits_LAB.mat")
% data_dim2_1 = data;

load("Lab_data_NelderMead2.mat")
data_dim2_1 = data;
load("test_many_dim_lineBO.mat")
data_dim2_2 = data(1,:);
load("test_many_dim_nelder_mead.mat")
data_dim2_3 = data;
 

fig = figure(3)
fig.Units='centimeters';
fig.Position(3:end)= [12.5,10.5];
%fig.Position(3:end)=[1050,650];
%fig.Position = [0,0, 1000,500];
hold on
y = data_dim1_1(:,9);
Y=getvals(y);
X=1:length(Y);
[yopt,xopt]=min(Y);
p1=plot(X,Y,'--','Color',[0 0.4470 0.7410],LineWidth=1.5)
p11 = plot(xopt,yopt,'*','Color',[0 0.4470 0.7410],'MarkerSize',10)
fill([X,flip(X,2)],[Y+std_Y,flip(Y-std_Y,2)],[0 0.4470 0.7410],'FaceAlpha',2*f_alpha,'EdgeColor','none');

y = data_dim1_2(:,2);
Y=getvals(y);
X=1:length(Y);
p2=plot(X,Y,'-','Color',[0 0.4470 0.7410],LineWidth=1.5)
[yopt,xopt]=min(Y);
p21 = plot(xopt,yopt,'*','Color',[0 0.4470 0.7410],'MarkerSize',10)
fill([X,flip(X,2)],[Y+std_Y,flip(Y-std_Y,2)],[0 0.4470 0.7410],'FaceAlpha',2*f_alpha,'EdgeColor','none');

y = data_dim1_3(:,1);
Y=getvals(y);
X=1:length(Y);
p3=plot(X,Y,'--','Color',[0.8500 0.3250 0.0980],LineWidth=1.5)
[yopt,xopt]=min(Y);
p31 = plot(xopt,yopt,'*','Color',[0.8500 0.3250 0.0980],'MarkerSize',10)
fill([X,flip(X,2)],[Y+std_Y,flip(Y-std_Y,2)],[0.8500 0.3250 0.0980],'FaceAlpha',2*f_alpha,'EdgeColor','none');


y = data_dim2_1(:,1);
Y=getvals(y);
X=1:length(Y);
p4=plot(X,Y,'Color',[0.8500 0.3250 0.0980],LineWidth=1.5)
[yopt,xopt]=min(Y);
p41 = plot(xopt,yopt,'*','Color',[0.8500 0.3250 0.0980],'MarkerSize',10)
fill([X,flip(X,2)],[Y+std_Y,flip(Y-std_Y,2)],[0.8500 0.3250 0.0980],'FaceAlpha',2*f_alpha,'EdgeColor','none');

% y = data_dim2_2(:,9);
% Y=getvals(y);
% X=1:length(Y);
% p5=plot(X,Y,'Color',[0 0.4470 0.7410],LineWidth=1.5)
% [yopt,xopt]=min(Y);
% p51 = plot(xopt,yopt,'*','Color',[0 0.4470 0.7410],'MarkerSize',10)
% 
% y = data_dim2_3(:,1);
% Y=getvals(y);
% X=1:length(Y);
% p6=plot(X,Y,'Color',[0.8500 0.3250 0.0980],LineWidth=1.5)
% [yopt,xopt]=min(Y);
% p61 = plot(xopt,yopt,'*','Color',[0.8500 0.3250 0.0980],'MarkerSize',10)
hold off
grid on
% xlim([0,2500])
% ylim([11.5,28])
ax = gca;
set(gca,'TickLabelInterpreter','latex');
%legend("D = 1; safe","D = 1; naive", "D = 1; safe \& optimized","D = 2; safe","D = 2; naive", "D = 2; safe \& optimized",'interpreter','latex','NumColumns',2)
xlabel("iteration $n$","Interpreter","latex")
ylabel("$$J_{opt}(n)$$ [fs]",'Interpreter','latex')
l1=legend([p1,p2],'Simulation','Laboratory','Interpreter','latex','NumColumns',3);
title(l1,'LineBO + MoSaOpt','Interpreter','latex')
a=axes('Position',get(ax,'position'),'visible','off');
ax.Box = 'on';
l2 = legend(a,[p3,p4],'Simulation','Laboratory','Interpreter','latex','NumColumns',3);
title(l2,'Nelder-Mead','Interpreter','latex')
l2.Units='centimeters';
l1.Units='centimeters';
l2.Position(1:2)=[l1.Position(1),7];
%%
d01=[2,7,8,9];
y=zeros(2,4);
y(1,:)=[data_dim2_3{d01,2}];
y(2,:)=[data_dim1_3{d01,2}];
fig=figure(3)
fig.Position(3:4)=[900,250];
b=bar(y,'grouped');
xtips=[b(1).XEndPoints;b(2).XEndPoints;b(3).XEndPoints;b(4).XEndPoints];
xtips=reshape(xtips,[1,8]);
y=reshape(y',1,8);
yt2=round(y,2);
yt2=num2cell(yt2);
yt2=cellfun(@num2str,yt2,'UniformOutput',false);
text(xtips,y+0.1,yt2,'HorizontalAlignment','center','FontSize',14,'VerticalAlignment','bottom')
for i=1:4
    b(i).EdgeColor='k';
    b(i).FaceColor=col(i);
end
ylabel("$J_{opt}$ [fs]")
ax=gca(figure(3));
ax.XTickLabel = {'$$\dim(\mathcal{L}) = 2$$', '$$\dim(\mathcal{L}) = 1$$'};
grid on

%%
% load("/home/jannis/master/master_thesis/matlab/own_lab/20220509_BayOpt_PAM_data/data_measurement/XFEL_PAM_measurement_corrected.mat") %PAM
%load("/home/jannis/master/master_thesis/matlab/own_lab/BAM_data/BAM_data.mat")
load("Lab_data_dim1_descent2.mat")
% x=data_dim2_3(:,2);
% x=cell2mat(x);
y=[];
x=[];
x01=[];
k=[1,2,3,5];
for i = 1:4
    id=~cellfun(@isempty,data{k(i),2});
    Y=data{k(i),2};
    t=Y(id);
    t=t{end};
    [y(i),I]=min(t);
    X=data{k(i),1};
    t=X(id);
    t=t{end};
    x(i,:)=t(I,:);
    x01(i,:) = t(1,:);
end
% d01=1;
% Kp_max = 3e1;
% Kp_min = 0.2;
% Ki_max = 6e1;
% Ki_min = 0;

% cond_t=[Kp_min, Kp_max;
%       Ki_min, Ki_max;
% %       1.2, 3.8;
%       0, 0.000105*350;
%       0, 3;
%       Kp_min, Kp_max;
%       Ki_min, Ki_max;
% %       1.2, 3.8;
%       0, 0.000105*350;
%       0, 3;
%       Kp_min, Kp_max;
%       Ki_min, Ki_max];

% for i=1:size(x,1)
%     x(i,:)=forwardCoordTransf(cond_t,x(i,:));
% end

% X0 = [20.6963   21.5537    0.0271    1.1841   20.5658   42.2428    0.0163    0.0587   10.0596   25.4586
%    11.4464   12.9611    0.0290    2.8479    9.9614   40.2759    0.0161    2.5005   23.1119   10.0352
%    23.7409   19.1115    0.0196    0.2699    3.5288    8.1776    0.0249    1.4855    5.8534   29.7003
%    17.3249   50.7107    0.0271    1.7580    7.5527   39.9850    0.0031    1.8779   19.8961   43.7851
%    25.3767   12.5643    0.0203    1.8897    1.1533   36.8828    0.0133    0.1486   14.7892   11.5506
%     3.8679   12.3297    0.0054    0.5672    1.4710   38.1119    0.0104    1.6158   20.9159   29.9470
%    12.6275   12.3585    0.0348    0.2462    3.3501    8.5225    0.0061    1.8629   17.2966    3.1247
%     5.4926   23.9154    0.0049    0.0927   28.1864   18.0784    0.0109    0.9988   14.1186   38.8919
%     0.9518   50.5324    0.0205    2.5623   10.5668   26.7616    0.0020    0.5313   19.9517   19.8497
%    24.5825    6.0133    0.0065    1.0789    1.8898   31.3131    0.0123    0.5270    6.4266   54.3092];

for i=1:size(x01,1)
    x01(i,:)=forwardCoordTransf(cond_t,x01(i,:),1,[0,1]);
end

% [x01,d01] = find_dispersedParams(x0);
% for i=1:size(x,1)
%     x(i,:)=backwardCoordTransf(cond_t,x(i,:));
% end
for i=1:size(x,1)
    x(i,:)=forwardCoordTransf(cond_t,x(i,:),1,[0,1]);
end
% [~,id]=min(y{end});
% id = d01;
%x01=x(1,:);
% x = x(id,:);
fig=figure(4);
fig.Units='centimeters';
fig.Position(3:4) = [12,12];
col = ['r','b','m','g'];
h1=subplot(2,2,1);
pos1 = h1.Position;
title('MENHIR')
hold on
hp = plot(x(1,1),x(1,2),'rx',x01(1,1),x01(1,2),'bo');
hp(1).MarkerSize = 8;
hp(2).MarkerSize = 8;
for i = 1:size(x,1)
    h(1:2,i)=plot(x(i,1),x(i,2),sprintf("%sx",col(i)),x01(i,1),x01(i,2),sprintf("%so",col(i)));
    h(1,i).MarkerSize = 8;
    h(2,i).MarkerSize = 8;
end
% h1.XScale='log';
% h1.YScale='log';
% h1.XMinorTick="off";
% h1.YMinorTick="off";
% h1.XMinorGrid="off";
% h1.YMinorGrid="off";
h1.PlotBoxAspectRatioMode='manual';
h1.PlotBoxAspectRatio=[1,1,1];
xlabel("P")
ylabel("I")
xlim([0,1]);
ylim([0,1])
grid on
hold off

x02 = x01;%find_dispersedParams(x0(:,3:4));
h2=subplot(2,2,2);
pos2 = h2.Position;
title('LSU 1.0')
hold on
for i = 1:size(x,1)
    h(3:4,i)=plot(x(i,3),x(i,4),sprintf("%sx",col(i)),x02(i,3),x02(i,4),sprintf("%so",col(i)));
    h(3,i).MarkerSize = 8;
    h(4,i).MarkerSize = 8;
end
% h2.XScale='log';
% h2.YScale='log';
% h2.XMinorTick="off";
% h2.YMinorTick="off";
% h2.XMinorGrid="off";
% h2.YMinorGrid="off";
h2.PlotBoxAspectRatioMode='manual';
xlabel("P")
ylabel("I")
xlim([0,1]);
ylim([0,1])
grid on
hold off

x03 = x01;%find_dispersedParams(x0(:,5:6));
h3=subplot(2,2,3);
pos3 = h3.Position;
title('SLO')
hold on
for i = 1:size(x,1)
    h(5:6,i)=loglog(x(i,5),x(i,6),sprintf("%sx",col(i)),x03(i,5),x03(i,6),sprintf("%so",col(i)));
    h(5,i).MarkerSize = 8;
    h(6,i).MarkerSize = 8;
end
% h3.XScale='log';
% h3.YScale='log';
% h3.XMinorTick="off";
% h3.YMinorTick="off";
% h3.XMinorGrid="off";
% h3.YMinorGrid="off";
h3.PlotBoxAspectRatioMode='manual';
xlabel("P")
ylabel("I")
xlim([0,1]);
ylim([0,1])
grid on
hold off

% x04 = x01;%find_dispersedParams(x0(:,7:8));
% h4=subplot(2,2,4);
% pos4 = h4.Position;
% title('LSU 2.2')
% hold on
% for i = 1:size(x,1)
%     h(7:8,i)=plot(x(i,7),x(i,8),sprintf("%sx",col(i)),x04(i,7),x04(i,8),sprintf("%so",col(i)));
%     h(7,i).MarkerSize = 8;
%     h(8,i).MarkerSize = 8;
% end
% h4.XScale='log';
% h4.YScale='log';
% h4.XMinorTick="off";
% h4.YMinorTick="off";
% h4.XMinorGrid="off";
% h4.YMinorGrid="off";
% h4.PlotBoxAspectRatioMode='manual';
% xlabel("P")
% ylabel("I")
% xlim([0.0001,1]);
% ylim([0.0001,1])
% grid on
% hold off
% 
% x05 = x01;%find_dispersedParams(x0(:,9:10));
% h5=subplot(2,2,5);
% pos5 = h5.Position;
% title('PPL')
% hold on
% for i = 1:size(x,1)
%     h(9:10,i)=plot(x(i,9),x(i,10),sprintf("%sx",col(i)),x05(i,9),x05(i,10),sprintf("%so",col(i)));
%     h(9,i).MarkerSize = 8;
%     h(10,i).MarkerSize = 8;
% end
% h5.XScale='log';
% h5.YScale='log';
% h5.XMinorTick="off";
% h5.YMinorTick="off";
% h5.XMinorGrid="off";
% h5.YMinorGrid="off";
% h5.PlotBoxAspectRatioMode='manual';
% xlabel("P")
% ylabel("I")
% xlim([0.0001,1]);
% ylim([0.0001,1])
% grid on
% hold off

lg=legend(hp,"optimum","start point",'numColumns',1,'location','northoutside','Fontsize',11 );
lg.Position(2)=0.95;


h1.Box = 'on';
h2.Box = 'on';
h3.Box = 'on';
% h4.Box = 'on';
% h5.Box = 'on';
% h1.Units='centimeters';
% h2.Units='centimeters';
% h3.Units='centimeters';
% h4.Units='centimeters';
% h5.Units='centimeters';
% h1.Position(2)=h1.Position(2)-0.01;
% h2.Position(2)=h2.Position(2)-0.015;
% h3.Position(2)=h3.Position(2)-0.02;
% h4.Position(2)=h4.Position(2)-0.025;
% h5.Position(2)=h5.Position(2)-0.03;
    

%%
load("data_dim1_coordinateOracle.mat")
data_dim1_4 = data;
load("data_dim1_descentOracle.mat")
data_dim1_5 = data;

fig = figure(1)
hold on
y = data_dim1_2(:,9);
Y=getvals(y);
X=1:length(Y);
plot(X,Y)
y = data_dim1_4(:,9);
Y=getvals(y);
X=1:length(Y);
plot(X,Y)
y = data_dim1_5(:,9);
Y=getvals(y);
X=1:length(Y);
plot(X,Y)
hold off
legend("random","coordinate","descent",'interpreter','latex')
xlabel("iteration n","Interpreter","latex")
ylabel("mean(jitter) [fs]",'Interpreter','latex')


function Y = getvals(y)
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
    len=max(cellfun('length',yt));
    for i = 1:length(yt)
        temp = ones(len,1)*min(yt{i});
        temp(1:length(yt{i})) = yt{i};
        for j=1:len
            yt1(i,j)=min(temp(1:j));
        end
    end
    Y = mean(yt1,1);
end

function varargout = find_dispersedParams(x0)
    t=x0;
    t=reshape(t,[10,1,2,5]);
    t1=permute(t,[2,1,3,4]);
    t2 = bsxfun(@minus,t1,t);
    t3=reshape(t2,100,2,5);
    t3=squeeze(vecnorm(t3,2,2));
    t3=reshape(t3,10,10,5);
    t3=sum(t3,1);
    t3=sum(t3,3);

    %t1=squeeze(vecnorm(t1-t,2,2));
%     combis = nchoosek(1:size(t1,1),4);
%     for i = 1:length(combis)
%         temp(i) = sum(t1(combis(i,1),combis(i,2:end))+t1(combis(i,2),combis(i,[1,3,4]))+t1(combis(i,3),combis(i,[1,2,4]))+t1(combis(i,4),combis(i,1:end-1)));
%     end
%     [~,combId] = max(temp);
%     vals = combis(combId,:);
    %t1=sum(t1,2);
    [~,vals] = sort(t3,'descend');
    vals = sort(vals(1:4));
    varargout{1}=x0(vals,:);
    varargout{2}=vals;
end