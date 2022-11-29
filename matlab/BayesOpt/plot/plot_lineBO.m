per = 0/100;
sharedGP = true;
globOpt =3;% 12.18;
val = 2;

subDIM = size(DIM,2)-2;
arg =find(DIM(:,end));
DIM = DIM(arg,:);
dimS = zeros(size(DIM,1),1);
Y0 = zeros(size(DIM,1),2);
n = 1;
lenDim = DIM(:,3);
if ~sharedGP
    for i = 1:length(arg)
        if i < length(arg)
            lenDim(i+1) = lenDim(i) + size(X{i},1);
        end
        for j = 1:size(X{i},1)
            xp(n) = n;
            y = Y{i};
            yp(n) = min(y(1:j));
            if j == 1
                Y0(i,1) = yp(n);
            end
            n = n + 1;
        end
        dimS(i,1) = DIM(i,1);
        Y0(i,2) = yp(n-1);
    end
else
    y = Y{length(DIM)};
    %y = Y;
    xp = 1:length(y);
    for i=1:length(xp)
        yp(i)=min(y(1:i));
    end
    lenDim = DIM(:,end);
    lenDim(1) = 0;
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
plot(p,[0,xp],ones(1,length(xp)+1)*globOpt)
l=legend("$$J(n)$$","$$J_{glob,opt}$$",'Interpreter','latex','Fontsize', 12);

%plot(I,min_y,'r+', 'MarkerSize',8)
t1=text(I,min_y,msg1,'HorizontalAlignment','center', 'VerticalAlignment', 'bottom', 'Interpreter','latex');
if ~isempty(temp)
    t2=text(I2,min_y2,msg2,'HorizontalAlignment','center', 'VerticalAlignment', 'bottom', 'Interpreter','latex');
    t2.FontSize = 12;
end
t1.FontSize = 12;

grid on
xlim([0 xp(end)])
ylim([yp(end)-0.5 yp(1)+0.5])
xlabel("Iteration n")
ylabel("Jitter [fs]")
ax1 = gca;
ax1.FontSize = 12;
ax1.Layer = 'top';
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
p2=plot(it_p,effspace,'o', 'Color', 'r');
grid on
hold on 
p1=plot(it_m,neffspace,'o', 'Color', 'b');
xlabel("Subspace Iteration K")
ylabel("Subspace L")
ax=gca;
ax.FontSize = 12;
yticks(0:size(cond,1)+1)
xticks(1:length(DIM))
ylim([0.5 size(cond,1)+0.5])
xlim([0,length(DIM)+1])
if subDIM == 1
    legend("Objective improved", "Objective not improved", 'Interpreter','latex','Fontsize',12,'Location','east')
else
    legend([p1(1),p2(1)],"Objective improved", "Objective not improved", 'Interpreter','latex','Fontsize',12,'Location','best')
end

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
