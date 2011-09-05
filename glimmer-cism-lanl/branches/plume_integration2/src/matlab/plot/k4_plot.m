f_ice6 = '/data/gladish/gc_output/igs_jobs/may29_highres_temp_-10.0_tau_25000.0_diff_10.0_k_4.0_amp_50.0_tempbot_-0.1_pmin_10.0/filtered.ice.nc';
f_plume6='/data/gladish/gc_output/igs_jobs/may29_highres_temp_-10.0_tau_25000.0_diff_10.0_k_4.0_amp_50.0_tempbot_-0.1_pmin_10.0/filtered.plume.nc';

%dice6 = nc_ice_read(f_ice6,1,51);
%dplume6 = nc_plume_read(f_plume6,1,51);

clf;

fs = 16;

subplot(1,2,1);
contourf(dice6.Xgrid/1000.0,dice6.Ygrid/1000.0,-dice6.lsurf(:,:,end),30,'EdgeColor','None');colorbar;
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);
title('Ice draft (m)','FontSize',fs);
caxis([0 600])
set(gca,'FontSize',fs);

subplot(1,2,2);
contourf(dice6.Xgrid/1000.0,dice6.Ygrid/1000.0,dice6.bmlt(:,:,end),30,'EdgeColor','None');colorbar;
xlabel('Across shelf distance (km)','FontSize',fs);
%ylabel('Along shelf distance (km)','FontSize',fs);
title('Basal melt rate (m/a)','FontSize',fs);
set(gca,'FontSize',fs);