clear all
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'defaultTextInterpreter','latex')
set(groot, 'defaultAxesFontsize',15)
set(groot, 'defaultTextFontsize',15)
set(groot, 'defaultLegendFontSize',14)
    
in_dir = "/home/jannis/master/master_thesis/matlab/own_lab/20220509_BayOpt_PAM_data/20220509_pam_data";
load("/home/jannis/master/master_thesis/matlab/own_lab/20220509_BayOpt_PAM_data/data_measurement/XFEL_PAM_measurement_corrected.mat")
fileInfo = dir(in_dir);
files = {fileInfo.name};
files = files(3:end);
id=strfind(files,'_');
arr = zeros(length(id),1);
for i = 1:length(id)
    temp_id = id{i};
    temp = files{i};
    arr(i) = str2double(temp(temp_id(end)+1:end-4));
end
[~,I] = sort(arr);
d_arr = diff (arr);
logic = d_arr == 40;
% count = 23:61;
% val = zeros(length(count),40);
% for i = count
%     val(i-count(1)+1,:)=rms(readPAM(in_dir+"/"+files{I(i)}),2)';
% end
% val_2m = zeros(floor(length(count)/2),1);
% std_dev_2 = zeros(size(val_2m));
% for i = 1:2:length(count)-1
%     val_2m(1+floor(i/2)) = mean(val(i:i+1,:),'all');
%     std_dev_2(1+floor(i/2)) = std(val(i:i+1,:),0,'all');
% end
% val_m = mean(val,2);
% std_dev = std(val,0,2);
% val_mm = mean(val_m);
% fig1 = figure(1)
% plot(val_m)
% hold on
% plot(ones(length(val_m),1)*val_mm)
% plot(2:2:length(count),val_2m)
% hold off
% 
% fig2 = figure(2)
% plot(std_dev)
% hold on
% plot(2:2:length(count),std_dev_2)
% hold off




% fig3 = figure(3);
% std_dev_meas = data{4};
% std_dev_meas = std_dev_meas(:,end);
% plot(std_dev_meas)
% title("Std Deviation of all Measurements")
Iopt=[1373623184+20;1373624723-20;1373625523+20];
Id_o = find(arr>Iopt(1) & arr<Iopt(2));
%Id_o = [Id_o;find(arr>Iopt(3))];

Istart = [1373624785+20;1373625523-20];
Id_s = find(arr>Istart(1) & arr<Istart(2));

k=1;
for i = (Id_s')
    val2_start(k,:)=std(readPAM(in_dir+"/"+files{I(i)}),1,2)';
    k=k+1;
end
k=1;
for i = (Id_o')
    val2_opt(k,:)=std(readPAM(in_dir+"/"+files{I(i)}),1,2)';
    k=k+1;
end
val2_start_m = mean(val2_start,2);
%%%% start before optimization %%%%%%
% count = 143:164;
% val = zeros(length(count),100);
% for i = count
%     val(i-count(1)+1,:)=rms(readPAM(in_dir+"/"+files{I(i)}),2)';
% end
% val2_start=val;
% val2_start_m = mean(val,2);
%%%%%                       %%%%%%%%%%
val2_opt_m = mean(val2_opt,2);
val2_start_m = val2_start_m(:);
x = linspace(0,35,1000);

val2_start = val2_start(:);
val2_opt = val2_opt(:);
%%
% val2_opt_m = ones(10,1);
% val2_start_m = ones(10,1);
% val2_opt = dist1(:);
% val2_start = dist2(:);

fig7=figure(7);
histogram(val2_start,20,'Normalization','pdf')
hold on
histogram(val2_opt,20,'Normalization','pdf')
pd = fitdist(val2_start,"Lognormal");
pd_opt = fitdist(val2_opt,"Lognormal");
pdf_val = pdf(pd,x);
pdf_opt = pdf(pd_opt,x);
plot(x,pdf_val,'b','LineWidth',1.5)
plot(x,pdf_opt,'r','LineWidth',1.5);
legend("$J_{0}$","$J_{opt}$","pdf $J_{0}$","pdf $J_{opt}$",'NumColumns',2)
xlim([0 35]);
hold off
xlabel("$J$ [fs]")
ylabel("$\phi(J)$")

% fig8=figure(8);
% histogram(val2_start_m,20,'Normalization','pdf')
% hold on
% histogram(val2_opt_m,20,'Normalization','pdf')
% pd = fitdist(val2_start_m,"Gamma");
% pd_opt = fitdist(val2_opt_m,"Gamma");
% pdf_val = pdf(pd,x);
% pdf_opt = pdf(pd_opt,x);
% plot(x,pdf_val,'b','LineWidth',1.5)
% plot(x,pdf_opt,'r','LineWidth',1.5);
% hold off

figure(9)
plot(x,pdf_val)
hold on
plot(x,pdf_opt);
legend("pdf before optimization","pdf after optimization")
hold off

gamma_opt = fitdist(val2_opt,"Lognormal");
gammaOpt_cdf = cdf(gamma_opt,x);
gamma_start = fitdist(val2_start(:),"Lognormal");
gammaInit_pdf = pdf(gamma_start,x);

int1_=trapz(x,gammaOpt_cdf.*gammaInit_pdf);

gamma_opt = fitdist(val2_opt,"Lognormal");
gammaOpt_cdf = cdf(gamma_opt,x);
gamma_start = fitdist(val2_start,"Lognormal");
gammaInit_pdf = pdf(gamma_start,x);

int2_=trapz(x,gammaOpt_cdf.*gammaInit_pdf);
int1_ = round(int1_,3);
int2_ = round(int2_,3);
figure(9)
bar([int1_,int2_])
text([1;2],[int1_;int2_],{num2str(int1_);num2str(int2_)},"HorizontalAlignment","center", "VerticalAlignment","bottom")
ylim([0,1.2])
legend()

%% compare rms and std dev
% data_str = data{5};
% val1 = zeros(length(data_str),100);
% val2 = zeros(size(val1));
% for i = 1:length(data_str)
%     val1(i,:)=rms(readPAM(in_dir+"/"+data_str{i}),2)';
%     val2(i,:)=std(readPAM(in_dir+"/"+data_str{i}),0,2)';
% end
% plot(mean(val1,2))
% hold on
% plot(mean(val2,2))
%%
%  y = data{2}{10};
%  i = 1;
%  j = 143;
%  j_old = j;
%  log_arr = zeros(1,length(files));
%  while i <= length(y)
%      upper = j +300;
%      while j <= upper
%          if mean(std(readPAM(in_dir+"/"+files{I(j)}),0,2)) == y(i)
%              j_old = j;
%              log_arr(j) = 1;
%              break;
%          end
%          j = j+1;
%      end
%      j = j_old;
%      i = i+1;
%  end
%  log_arr = logical(log_arr);
%  %%
%  c = 143:897;
%  y_all = zeros(size(c));
%  for i = c
%      y_all(i-c(1)+1) = mean(std(readPAM(in_dir+"/"+files{I(i)}),0,2));
%  end
%  fig7=figure(7);
%  fig7.Position = [0,0, 2400,1125];
%  title("all mean values during optimization")
%  plot(c,y_all)
%  hold on
%  plot(c(log_arr(143:897)),y_all(log_arr(143:897)),'r+')
%  legend("mean","mean used by algorithm")
%  hold off
     