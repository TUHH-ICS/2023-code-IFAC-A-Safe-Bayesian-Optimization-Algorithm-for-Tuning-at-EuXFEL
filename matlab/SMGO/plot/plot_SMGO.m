set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'defaultTextInterpreter','latex')
set(groot, 'defaultAxesFontsize',12)
set(groot, 'defaultTextFontsize',12)
per = 0/100;
sharedGP = true;
globOpt =3;% 12.18;
val = 1;
X = data{1};
Y = data{2};
y = Y{length(DIM)};
%y = Y;
xp = 1:length(y);
imp = zeros(1,length(y));
for i=1:length(xp)
    [yp(i),Iopt(i)]=min(y(1:i));
    if i > 1 && yp(i) ~= yp(i-1)
        imp(i)=1;
    end
end


it = 1:length(arg);
[min_y,I] = min(yp);
temp = find(yp <= (1.0+per)*globOpt);
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

fig = figure(1);
fig.Position = [0,0, 2000,1125];
p = subplot(2,1,1);
plot(p,xp,yp)
hold on 
plot(p,xp,yp+std_dev(Iopt)','k--')
plot(p,xp,yp-std_dev(Iopt)','k--')
plot(p,xp(logical(imp)),yp(logical(imp)),'r+')
hold off
%plot(p,[0,xp],ones(1,length(xp)+1)*globOpt)
l=legend("$$J(n)$$","$J(n) \pm \sigma$","","improvements",'Interpreter','latex','Fontsize', 12);

%plot(I,min_y,'r+', 'MarkerSize',8)
t1=text(I,min_y,msg1,'HorizontalAlignment','center', 'VerticalAlignment', 'bottom', 'Interpreter','latex');
if ~isempty(temp)
    t2=text(I2,min_y2,msg2,'HorizontalAlignment','center', 'VerticalAlignment', 'bottom', 'Interpreter','latex');
    t2.FontSize = 12;
end
t1.FontSize = 12;

grid on
xlim([0 xp(end)])
ylim([yp(end)-std_dev(Iopt(end))'-0.5 yp(1)+std_dev(Iopt(1))+0.5])
xlabel("Iteration n")
ylabel("Jitter [fs]")

ax1 = gca;
ax1.FontSize = 12;
ax1.Layer = 'top';
%ax1.XTick = 1:10:length(yp);
ax2 = axes('Position',p.Position,'XAxisLocation','top','YAxisLocation','right','Layer','bottom');
ax2.YAxis.Visible = 'off';
ax2.XLim = ax1.XLim;
ax2.XTick = lenDim(1:end);
for i=1:length(lenDim)
    if i == 1
        xticlabel{i} = num2str(i);
        continue;
    end
    if mod(i,val) == 0
        xticlabel{i} = num2str(i);
    else
        xticlabel{i} = '';
    end
end
ax2.XTickLabel =xticlabel;
ax2.TickLength = ax1.TickLength/2;
ax2.FontSize=8;
ax2.YTick = [];
ax2.YLim = ax1.YLim;
ax2.XGrid = 'on';
ax2.GridLineStyle = ':';
ax2.Color = 'none';
%p.Box = 'off';
ax2.Box = 'off';
xlabel('Subspace Iteration K');


logic = logical(DIM(:,subDIM+1));
effspace = DIM(logic,1:subDIM);
neffspace = DIM(~logic,1:subDIM);
it_p = it(logic);
it_m = it(~logic);
p=subplot(2,1,2);
plot(it_p,effspace,'o', 'Color', 'r')
grid on
hold on 
plot(it_m,neffspace,'o', 'Color', 'b')
xlabel("Subspace Iteration K")
ylabel("Subspace L")
ax=gca;
ax.FontSize = 12;
yticks(0:size(cond,1)+1)
xticks(1:length(DIM))
ylim([0.5 size(cond,1)+0.5])
xlim([0,length(DIM)+1])
if subDIM == 1
    legend("$$\mbox{Objective improved}$$", "$$\mbox{Objective not improved}$$", 'Interpreter','latex','Fontsize',12,'Location','east')
else
    legend("$$\mbox{Objective improved}$$","", "$$\mbox{Objective not improved}$$", 'Interpreter','latex','Fontsize',12,'Location','best')
end

fig2 = figure(2);
plot(1:length(Iopt),std_dev(Iopt))
ylabel("$\sigma(n)$")
xlabel("Iteration n")
grid on
%%
threshold = 30;
for i=1:length(y)
        yp2(i)=max(y(1:i));
end

figure(2)
plot(xp,yp2,xp,ones(length(xp),1)*threshold)
xlim([0 xp(end)])
ylim([yp2(1)-0.5 max(yp2(end),threshold)+0.5])
xlabel("Iteration n")
ylabel("Jitter [fs]")
l=legend("$$J_{max}(n)$$","$$J_{T}$$",'Interpreter','latex','Fontsize',12,'Location','best');
grid on
