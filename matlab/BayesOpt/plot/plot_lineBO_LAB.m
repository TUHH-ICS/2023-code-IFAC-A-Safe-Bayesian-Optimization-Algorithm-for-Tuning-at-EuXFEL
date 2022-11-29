clear all
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'defaultTextInterpreter','latex')
set(groot, 'defaultAxesFontsize',12)
set(groot, 'defaultTextFontsize',12)

load("/home/jannis/master/master_thesis/matlab/own_lab/20220509_BayOpt_PAM_data/data_measurement/XFEL_PAM_measurement_corrected.mat")
data_nr = 1;
xt = data(data_nr,1);
xt = xt{1,1};
yt = data(data_nr,2);
yt = yt{1,1};
% t = data(data_nr,4);
% t=t{1,1};
% jitters = cell2mat(t(1));
% jitters = jitters(:,1:end-1);
% t = t(:,end);
% t=cell2mat(t);
% t=round(t);
% t=t+7;
% t(607:end)=t(607:end)-320;
jitters = data{data_nr,4};
c=cellfun("isempty",yt);
yt = yt(~c);
xt = xt(~c);
xt=xt{end};
yt=yt{end};

xp=1:size(yt,1);
yp = zeros(size(yt));
jit = zeros(size(jitters));
for i = 1:size(yt,1)
    [yp(i),I] = min(yt(1:i));
    jit(i,:) = jitters(I,:);
end

val = 40;
[min_y,I] = min(yp);
temp=[];
%temp = find(yp <= (1.0+per)*globOpt);
%msg1_1 = sprintf('$$\\left(n_{opt} / J_{opt}\\right)$$');
msg2_1 = sprintf("$$\\left(\\mbox{n}_{opt} / \\mbox{J}_{opt}\\right) \\widehat{=} \\left(%d / %.3f\\right)$$",I,min_y);
msg1 = {msg2_1, '$$ \big\downarrow$$'};
if ~isempty(temp)
    I2 = temp(1);
    min_y2 = yp(I2);
    msg2_1 = sprintf("$$\\left(\\mbox{n}_{%3$.1f%4$s} / \\mbox{J}_{%3$.1f%4$s}\\right) \\widehat{=} \\left(%1$d / %2$.3f\\right)$$",I2,min_y2,per*100,'\%');
    msg2 = {msg2_1, '$$ \big\downarrow$$'};
    %msg2 = sprintf("$$\\longleftarrow \\left(\\mbox{n}_{%3$.1f%4$s} / \\mbox{J}_{%3$.1f%4$s}\\right) \\widehat{=} \\left(%1$d / %2$.3f\\right)$$",I2,min_y2,per*100,'\%');
else
    msg2={};
end

% tick=[];
% for i=1:length(t)
%     if i == 1
%         xticlabel{i} = num2str(t(i));
%         tick(end+1) = i;
%         continue;
%     end
%     if mod(i,val) == 0
%         xticlabel{1+floor(i/val)} = num2str(t(i));
%         tick(end+1) = i;
% %     else
% %         xticlabel{i} = '';
%     end
% end
fig1 = figure(1);
fig1.Position = [0,0, 1800,800];
p=plot(xp,yp);

%plot([0,xp],ones(1,length(xp)+1)*globOpt)
%l=legend("$$J(n)$$","$$J_{glob,opt}$$",'Interpreter','latex','Fontsize', 12);

%plot(I,min_y,'r+', 'MarkerSize',8)
t1=text(I,min_y,msg1,'HorizontalAlignment','left', 'VerticalAlignment', 'bottom', 'Interpreter','latex');
if ~isempty(temp)
    t2=text(I2,min_y2,msg2,'HorizontalAlignment','center', 'VerticalAlignment', 'bottom', 'Interpreter','latex');
    t2.FontSize = 12;
end
t1.FontSize = 12;

grid on
xlim([0 xp(end)+20])
ylim([yp(end)-0.5 yp(1)+0.5])
xlabel("Iteration n")
ylabel("Jitter [fs]")
ax1 = gca;
ax1.XTick = tick;
ax1.FontSize = 12;
ax1.Layer = 'top';
ax2 = axes('Position',ax1.Position,'XAxisLocation','top','YAxisLocation','right','Layer','bottom');
ax2.YAxis.Visible = 'off';
ax2.XLim = ax1.XLim;

ax2.XTick = tick;
ax2.XTickLabel =xticlabel;
ax2.TickLength = ax1.TickLength/2;
ax2.FontSize=12;
ax2.YTick = [];
ax2.YLim = ax1.YLim;
ax2.XGrid = 'on';
ax2.GridLineStyle = ':';
ax2.Color = 'none';
%p.Box = 'off';
ax2.Box = 'off';
xlabel('time [s]');
%%

min_jit = (min(jit,[],1)'.*ones(size(xp)))';
fig2=figure(2);
fig2.Position=[1000,1564,560,620];
p=subplot(3,1,1);
plot(p,xp,jit(:,1),xp,min_jit(:,1),'k--')
legend("MENHIR $J(n)$","MENHIR $J_{min}$",Location="best")
p=subplot(3,1,2);
plot(p,xp,jit(:,2),xp,min_jit(:,2),'k--')
legend("LINK $J(n)$","LINK $J_{min}$",Location="best")
p=subplot(3,1,3);
plot(p,xp,jit(:,3),xp,min_jit(:,3),'k--')
legend("ORIGAMI $J(n)$","ORIGAMI $J_{min}$",Location="best")