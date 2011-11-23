clear all

jobs = getenv('GC_JOBS');

j1name = 'oct30_perturb_usq.8_thk_perturb_k_7.0_amp_25.0';

f1_ice = strcat([jobs,'/',j1name,'/',j1name,'.out.nc']);
f1_ocean = strcat([jobs,'/',j1name,'/plume.',j1name,'.out.nc']);

f2_ice = strcat([jobs,'/',j1name,'/','no_homotopy.out.nc']);
f2_ocean = strcat([jobs,'/',j1name,'/plume.no_homotopy.nc']);

dice1 = nc_ice_read(f1_ice,1,20,421);
dice2 = nc_ice_read(f2_ice,20,1,21);

dplume1 = nc_plume_read(f1_ocean,1,20,421);
dplume2 = nc_plume_read(f2_ocean,20,1,21);

m_control = 0.789;
[m,in,out,acab,unsteady] = mass_balance(dice1,dplume1);
m1 = m/in
[m,in,out,acab,unsteady] = mass_balance(dice2,dplume2);
m2 = m/in

fig1 = figure(1);
set(fig1,'Position',[1 1 1200 800]);
fs = 16;
fs_label=18;
clf;

subset = @(M) M(:,1:120);

subplot(2,2,1);
hold on
set(gca,'FontSize',fs);
xlabel('across shelf distance (km)','FontSize',fs);
ylabel('along shelf distance (km)','FontSize',fs);
contourf(subset(dice1.Xgrid/1000),subset(dice1.Ygrid/1000), ...
         subset(dice1.bmlt(:,:,2)),30,'EdgeColor','None');
colorbar('FontSize',fs);
caxis([0 80]);

subplot(2,2,2);
hold on
set(gca,'FontSize',fs);
xlabel('across shelf distance (km)','FontSize',fs);
ylabel('along shelf distance (km)','FontSize',fs)
hold on
contourf(subset(dice1.Xgrid/1000),subset(dice1.Ygrid/1000),...
         subset(dice2.bmlt(:,:,1)),30,'EdgeColor','None');
colorbar('FontSize',fs);
caxis([0 80]);

subplot(2,2,3);
hold on
set(gca,'FontSize',fs);
xlabel('across shelf distance (km)','FontSize',fs);
ylabel('along shelf distance (km)','FontSize',fs)
hold on
contourf(subset(dice1.Xgrid/1000),subset(dice1.Ygrid/1000), ...
         subset(dice1.lsurf(:,:,2)),30,'EdgeColor','None');
colorbar('FontSize',fs);
caxis([-550 0]);

subplot(2,2,4);
hold on
set(gca,'FontSize',fs);
xlabel('across shelf distance (km)','FontSize',fs);
ylabel('along shelf distance (km)','FontSize',fs)
contourf(subset(dice2.Xgrid/1000),subset(dice2.Ygrid/1000), ...
         subset(dice2.lsurf(:,:,1)),30,'EdgeColor','None');
colorbar('FontSize',fs);
caxis([-550 0]);
