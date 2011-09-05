fice_min = '/data/gladish/gc_output/paper_jobs/pmin/last_time.aug27_pmin.out.nc';
fplume_min = '/data/gladish/gc_output/paper_jobs/pmin/last_time.plume.aug27_pmin.out.nc';

fice_ustar = '/data/gladish/gc_output/paper_jobs/ustar/last_time.aug28_ustar.out.nc';
fplume_ustar = '/data/gladish/gc_output/paper_jobs/ustar/last_time.plume.aug28_ustar.out.nc';

fice_sgd = '/data/gladish/gc_output/paper_jobs/sgd/last_time.aug28_sgd.out.nc';
fplume_sgd = '/data/gladish/gc_output/paper_jobs/sgd/last_time.plume.aug28_sgd.out.nc';

fice_min = '/data/gladish/gc_output/1djobs_fromcluster/pmin/aug31_m20.out.nc';
fplume_min = '/data/gladish/gc_output/1djobs_fromcluster/pmin/plume.aug31_m20.out.nc';

fice_ustar = '/data/gladish/gc_output/1djobs_fromcluster/ustar/aug28_ustar.out.nc';
fplume_ustar = '/data/gladish/gc_output/1djobs_fromcluster/ustar/plume.aug28_ustar.out.nc';

fice_sgd = '/data/gladish/gc_output/1djobs_fromcluster/sgd/aug28_sgd.out.nc';
fplume_sgd = '/data/gladish/gc_output/1djobs_fromcluster/sgd/plume.aug28_sgd.out.nc';

dplume_min = nc_plume_read(fplume_min,1,1);
dice_min =   nc_ice_read(fice_min,1,1);
dplume_ustar = nc_plume_read(fplume_ustar,1,1);
dice_ustar =   nc_ice_read(fice_ustar,1,1);
dplume_sgd = nc_plume_read(fplume_sgd,1,1);
dice_sgd =   nc_ice_read(fice_sgd,1,1);

m10_fice_min = '/data/gladish/gc_output/paper_jobs/pmin/last_time.restart_t_m10.out.nc';
m10_fplume_min = '/data/gladish/gc_output/paper_jobs/pmin/last_time.plume.restart_t_m10.out.nc';

m10_fice_ustar = '/data/gladish/gc_output/paper_jobs/ustar/last_time.ustar_restart_m10.out.nc';
m10_fplume_ustar = '/data/gladish/gc_output/paper_jobs/ustar/last_time.plume.ustar_restart_m10.out.nc';

m10_fice_sgd = '/data/gladish/gc_output/paper_jobs/sgd/last_time.sgd_restart_m10.out.nc';
m10_fplume_sgd = '/data/gladish/gc_output/paper_jobs/sgd/last_time.plume.sgd_restart_m10.out.nc';

m10_fice_min = '/data/gladish/gc_output/1djobs_fromcluster/pmin/restart_t_m10_uniform_restart.out.nc';
m10_fplume_min = '/data/gladish/gc_output/1djobs_fromcluster/pmin/plume.restart_t_m10_uniform_restart.out.nc';

m10_fice_ustar = '/data/gladish/gc_output/1djobs_fromcluster/ustar/aug31_m10_ustar.out.nc';
m10_fplume_ustar = '/data/gladish/gc_output/1djobs_fromcluster/ustar/plume.aug31_m10_ustar.out.nc';

m10_fice_sgd = '/data/gladish/gc_output/1djobs_fromcluster/sgd/aug31_m10_sgd.out.nc';
m10_fplume_sgd = '/data/gladish/gc_output/1djobs_fromcluster/sgd/plume.aug31_m10_sgd.out.nc';


m10_dplume_min = nc_plume_read(m10_fplume_min,1,1);
m10_dice_min =   nc_ice_read(m10_fice_min,1,1);
m10_dplume_ustar = nc_plume_read(m10_fplume_ustar,1,1);
m10_dice_ustar =   nc_ice_read(m10_fice_ustar,1,1);
m10_dplume_sgd = nc_plume_read(m10_fplume_sgd,1,1);
m10_dice_sgd =   nc_ice_read(m10_fice_sgd,1,1);

lw = 1.5;
fs = 16;
fs2 = 16;

ymax = 25;
xmax = 62;
xmin = -2;

clf;

subplot(2,1,1);
hold on
plot(dice_min.ygrid/1000.0,mean(dplume_min.bmelt(:,:,end),1), ...
     'Linewidth',lw,'LineStyle','-','color','k');
plot(dice_ustar.ygrid/1000.0,mean(dplume_ustar.bmelt(:,:,end),1),...
     'Linewidth',lw,'LineStyle','-','color','b');
plot(dice_sgd.ygrid/1000.0,mean(dplume_sgd.bmelt(:,:,end),1),...
     'Linewidth',lw,'LineStyle','-','color','r');
plot(m10_dice_min.ygrid/1000.0,mean(m10_dplume_min.bmelt(:,:,end),1), ...
     'Linewidth',lw,'LineStyle','*','color','k');
plot(m10_dice_ustar.ygrid/1000.0,mean(m10_dplume_ustar.bmelt(:,:,end),1),...
     'Linewidth',lw,'LineStyle','*','color','b');
plot(m10_dice_sgd.ygrid/1000.0,mean(m10_dplume_sgd.bmelt(:,:,end),1),...
     'Linewidth',lw,'LineStyle','*','color','r');
xlim([xmin xmax]);
ylim([0 ymax]);
set(gca,'YColor','k');
set(gca,'YTick',0:5:ymax);
ylabel('basal melt rate (m/a)','FontSize',fs)
set(gca,'fontsize',fs2);
xlabel('along shelf distance (km)','FontSize',fs);

subplot(2,1,2);
hold on
plot( dice_min.ygrid/1000.0,mean(dice_min.lsurf(:,:,end),1),...
       'Linewidth',lw,'LineStyle','-','color','k');
plot( dice_ustar.ygrid/1000.0,mean(dice_ustar.lsurf(:,:,end),1),...
	   'Linewidth',lw,'LineStyle','-','color','b');
plot(dice_sgd.ygrid/1000.0,mean(dice_sgd.lsurf(:,:,end),1),...
     'Linewidth',lw,'LineStyle','-','color','r');
plot( m10_dice_min.ygrid/1000.0,mean(m10_dice_min.lsurf(:,:,end),1),...
       'Linewidth',lw,'LineStyle','*','color','k');
plot( m10_dice_ustar.ygrid/1000.0,mean(m10_dice_ustar.lsurf(:,:,end),1),...
	   'Linewidth',lw,'LineStyle','*','color','b');
plot(m10_dice_sgd.ygrid/1000.0,mean(m10_dice_sgd.lsurf(:,:,end),1),...
     'Linewidth',lw,'LineStyle','*','color','r');

xlim([xmin xmax]);
ylim([-600 0]);

set(gca,'YTick',-600:100:0);
set(gca,'YColor','k');
xlabel('along shelf distance (km)','FontSize',fs);
ylabel('ice draft (m)','FontSize',fs)
set(gca,'fontsize',fs2);


%legend([h1 h3 h5 h2 h4 h6], ...
%       'basal melt (min. thickness)',...
%       'basal melt (u_* offset)', ...
%       'basal melt (sub-glacial discharge)', ...
%       'ice draft (min thickness)', ...
%       'ice draft (u_* offset)',...       
%       'ice draft (sub-glacial discharge)',...
%       'Location','NorthEast');

