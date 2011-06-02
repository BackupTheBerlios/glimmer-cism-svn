function [b,bint,stats,linmelt_out,resid_out] = may7_work(docean,dice)

figure(1);
clf;

fs = 15;
[m1,n1,k] = size(docean.bmelt(:,:,:));
%m1=m;
%n1=n;
m=m1-2;
n=n1-2;
t_stride=20;
times = 1:t_stride:k;
k = length(times);
lasttime = 806;
remove_edge=@(d) d(2:(m1-1),2:(n1-1),:);

num_plots = 6;

flatbmlt = reshape(remove_edge(docean.bmelt(:,:,times)),m*n*k,1);
flatspeed = reshape(sqrt(remove_edge(docean.su(:,:,times)) .^2 ...
                        +remove_edge(docean.sv(:,:,times)) .^ 2),m*n*k,1);
flattemp = reshape(remove_edge(docean.temp(:,:,times)),m*n*k,1);
flatsalt = reshape(remove_edge(docean.salt(:,:,times)),m*n*k,1);
flatgradx = reshape(remove_edge(docean.gradx(:,:,times)),m*n*k,1);
flatgrady = reshape(remove_edge(docean.grady(:,:,times)),m*n*k,1);
flatgrad = reshape(remove_edge(docean.grad(:,:,times)),m*n*k,1);
flatdraft = reshape(remove_edge(docean.draft(:,:,times)),m*n*k,1);
flatadv = reshape(remove_edge(-dice.flux_div(:,:,times) + ...
                               dice.vel_div(:,:,times)),m*n*k,1);
flatadv = reshape(remove_edge(-dice.flux_div(:,:,times)),m*n*k,1);
flatadv = reshape(remove_edge(-dice.y_adv(:,:,times)),m*n*k,1);                           

subplot(2,num_plots/2,1);
hold on
plot(flatspeed,flatbmlt,'b.');
set(gca,'FontSize',fs);
xlabel('plume speed (m/s)','FontSize',fs);
ylabel('melt rate (m/year)','FontSize',fs);
title('corr = 0.94','FontSize',fs);
subplot(2,num_plots/2,2);
hold on
xlabel('plume temperature (C)','FontSize',fs);
ylabel('melt rate (m/year)','FontSize',fs);
title('corr = 0.17','FontSize',fs);
plot(flattemp,flatbmlt,'b.');
set(gca,'FontSize',fs);
subplot(2,num_plots/2,3);
hold on
xlabel('plume salinity (psu)','FontSize',fs);
ylabel('melt rate (m/year)','FontSize',fs);
title('corr = -0.54','FontSize',fs);
plot(flatsalt,flatbmlt,'b.');
set(gca,'FontSize',fs);
subplot(2,num_plots/2,4);
hold on
xlabel('ice gradient (non-dimensional)','FontSize',fs);
ylabel('melt rate (m/year)','FontSize',fs);
title('corr = 0.46','FontSize',fs);
plot(flatgrad,flatbmlt,'b.');
set(gca,'FontSize',fs);
subplot(2,num_plots/2,5);
hold on
xlabel('v_0 H_y (m/year)','FontSize',fs);
xlim([-50 50]);
ylabel('melt rate (m/year)','FontSize',fs);
title('corr = 0.91','FontSize',fs);
plot(flatadv,flatbmlt,'b.');

set(gca,'FontSize',fs);
subplot(2,num_plots/2,6);
hold on
xlabel('ice draft (m)','FontSize',fs);
ylabel('melt rate (m/year)','FontSize',fs);
title('corr = 0.18','FontSize',fs);
plot(-flatdraft,flatbmlt,'b.');

set(gca,'FontSize',fs);
corrcoef([flatbmlt,flatspeed,flattemp,flatsalt,flatgrad,-flatdraft,flatadv])
%b
%bint
%stats
[b,bint,r,rint,stats] = regress(flatbmlt,[ones(m*n*k,1), ...
                                          flatgrad, ...
                                          flatgrad.^2, ...
                                          flatdraft, ...
                                          flatdraft.^2, ...
                                          flatadv]);

lookback = 10;
flatbmlt = reshape(remove_edge(mean(docean.bmelt(:,:,(lasttime-lookback):lasttime),3)),m*n,1);
flatspeed = reshape(remove_edge(mean(sqrt(docean.su(:,:,(lasttime-lookback):lasttime).^2 ...
                                         +docean.sv(:,:,(lasttime-lookback):lasttime).^2),3)),m*n,1);
flattemp = reshape(remove_edge(mean(docean.temp(:,:,(lasttime-lookback):lasttime),3)),m*n,1);
flatsalt = reshape(remove_edge(mean(docean.salt(:,:,(lasttime-lookback):lasttime),3)),m*n,1);
flatgrad = reshape(remove_edge(mean(docean.grad(:,:,(lasttime-lookback):lasttime),3)),m*n,1);
flatgradx = reshape(remove_edge(mean(docean.gradx(:,:,(lasttime-lookback):lasttime),3)),m*n,1);
flatgrady = reshape(remove_edge(mean(docean.grady(:,:,(lasttime-lookback):lasttime),3)),m*n,1);
flatdraft = reshape(remove_edge(mean(docean.draft(:,:,(lasttime-lookback):lasttime),3)),m*n,1);
flatadv = reshape(remove_edge(mean(-dice.flux_div(:,:,(lasttime-lookback):lasttime) + ...
                                    dice.vel_div(:,:,(lasttime-lookback):lasttime),3)),m*n,1);
flatadv = reshape(remove_edge(mean(1000*docean.grady(:,:,(lasttime-lookback):lasttime),3)),m*n,1);
%-dice.y_adv(:,:,(lasttime-lookback):lasttime),3)),m*n,1);
subplot(2,num_plots/2,1);
plot(flatspeed,flatbmlt,'r.');
set(gca,'FontSize',fs);
hold off
subplot(2,num_plots/2,2);
plot(flattemp,flatbmlt,'r.');
set(gca,'FontSize',fs);
hold off
subplot(2,num_plots/2,3);
plot(flatsalt,flatbmlt,'r.');
set(gca,'FontSize',fs);
hold off
subplot(2,num_plots/2,4);
plot(flatgrad,flatbmlt,'r.');
hold off
subplot(2,num_plots/2,5);
set(gca,'FontSize',fs);
plot(flatadv,flatbmlt,'r.');
plot(0:0.5:40,0:0.5:40,'k.');
hold off
subplot(2,num_plots/2,6);
plot(-flatdraft,flatbmlt,'r.');
set(gca,'FontSize',fs);
hold off

size(flatgradx)
size(flatdraft)

corrcoef([flatbmlt,flatspeed,flattemp,flatsalt,flatgrad,flatgradx,flatgrady,flatadv,-flatdraft])
[b,bint,r,rint,stats] = regress(flatbmlt,[ones(m*n,1), ...
                                          flatgrad, ...
                                          flatgrad.^2, ...
                                          flatgradx, ...
                                          flatgradx.^2,...
                                          flatgrady, ...
                                          flatgrady.^2,...
                                          flatdraft, ...
                                          flatdraft.^2, ...
                                          flatadv]);

linmelt = ones(m*n,1)*b(1) + b(2)*flatgrad + b(3)*flatgrad.^2 ...
                           + b(4)*flatgradx+ b(5)*flatgradx.^2 ...
                           + b(6)*flatgrady+ b(7)*flatgrady.^2 ...                           
                           + b(8)*flatdraft +b(9)*flatdraft.^2 ...
                           + b(10)*flatadv;
%b
%bint
%stats
linmelt_out = zeros(m1,n1);
linmelt_out(2:(m1-1),2:(n1-1)) = reshape(linmelt,m,n);

resid_out = zeros(m1,n1);
resid_out(2:(m1-1),2:(n1-1)) = reshape(r,m,n);
end