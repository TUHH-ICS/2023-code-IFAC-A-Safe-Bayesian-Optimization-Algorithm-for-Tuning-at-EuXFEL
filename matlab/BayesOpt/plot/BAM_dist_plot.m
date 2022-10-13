set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'defaultTextInterpreter','latex')
set(groot, 'defaultAxesFontsize',15)
set(groot, 'defaultTextFontsize',15)
set(groot, 'defaultLegendFontSize',14)

load("/home/jannis/master/master_thesis/matlab/own_lab/BAM_data/jitter_opt.mat")
dist1 = jitter;
load("/home/jannis/master/master_thesis/matlab/own_lab/BAM_data/jitter_start.mat")
dist2 = jitter;

opt_dist = fitdist(dist1,'Lognormal');
start_dist = fitdist(dist2,'EpsilonSkewNormal');
x=0:0.01:35;
pdf_opt = pdf(opt_dist,x);
cdf_opt = cdf(opt_dist,x);
pdf_start = pdf(start_dist,x);

histogram(dist1,20,'Normalization','pdf')
hold on
histogram(dist2,20,'Normalization','pdf')
plot(x,pdf_opt,'b-','LineWidth',1.5)
plot(x,pdf_start,'r-','LineWidth',1.5)
legend("$J_{opt}$","$J_{0}$","pdf $J_{opt}$","pdf $J_{0}$",'NumColumns',2)
xlim([0 35]);
hold off
xlabel("$J$ [fs]")
ylabel("$\phi(J)$")

dist1_s=sort(dist1,'ascend');
dist2_s =sort(dist2,'ascend');

p=0;
for i = 1:length(dist2_s)
    p=p+sum(dist1_s < dist2_s(i))/(length(dist2_s)*length(dist1_s));
end
p=p*100;

p1=trapz(cdf_opt.*pdf_start);