%focean_min = '/data/gladish/gc_output/paper_jobs/no_rot_jul4/no_rot_jul4.out.nc';
%fplume_min = '/data/gladish/gc_output/paper_jobs/no_rot_jul4/plume.no_rot_jul4.out.nc';
focean_min = '/data/gladish/gc_output/paper_jobs/pmin/last_time.aug27_pmin.out.nc';
fplume_min = '/data/gladish/gc_output/paper_jobs/pmin/last_time.plume.aug27_pmin.out.nc';

%focean_ustar = '/data/gladish/gc_output/paper_jobs/no_rot_jul4_restart_2/no_rot_jul4_restart_2.out.nc';
%fplume_ustar = '/data/gladish/gc_output/paper_jobs/no_rot_jul4_restart_2/plume.no_rot_jul4_restart_2.out.nc';
focean_ustar = '/data/gladish/gc_output/paper_jobs/ustar/last_time.aug28_ustar.out.nc';
fplume_ustar = '/data/gladish/gc_output/paper_jobs/ustar/last_time.plume.aug28_ustar.out.nc';

%focean_sgd = '/data/gladish/gc_output/paper_jobs/no_rot_jul4_restart_2/no_rot_jul4_restart_2.out.nc';
%fplume_sgd = '/data/gladish/gc_output/paper_jobs/no_rot_jul4_restart_2/plume.no_rot_jul4_restart_2.out.nc';
focean_sgd = '/data/gladish/gc_output/paper_jobs/sgd/last_time.aug28_sgd.out.nc';
fplume_sgd = '/data/gladish/gc_output/paper_jobs/sgd/last_time.plume.aug28_sgd.out.nc';


%dplume_min = nc_plume_read(fplume_min,1,205);
%dice_min =   nc_ice_read(focean_min,1,205);
%dplume_ustar = nc_plume_read(fplume_ustar,1,161);
%dice_ustar =   nc_ice_read(focean_ustar,1,161);
%dplume_sgd = nc_plume_read(fplume_sgd,1,161);
%dice_sgd =   nc_ice_read(focean_sgd,1,161);

dplume_min = nc_plume_read(fplume_min,1,1);
dice_min =   nc_ice_read(focean_min,1,1);
dplume_ustar = nc_plume_read(fplume_ustar,1,1);
dice_ustar =   nc_ice_read(focean_ustar,1,1);
dplume_sgd = nc_plume_read(fplume_sgd,1,1);
dice_sgd =   nc_ice_read(focean_sgd,1,1);

lw = 2.0;
fs = 18;

clf;
hold on

[ax1,h1,h2] = plotyy(dice_min.ygrid/1000.0, ...
                    mean(dplume_min.bmelt(:,:,end),1),...
                    dice_min.ygrid/1000.0, ...
                    mean(dice_min.lsurf(:,:,end),1));

set(ax1,'ycolor','k')
set(h1,'Linewidth',lw,'LineStyle','-','color','k');
set(h2,'Linewidth',lw,'LineStyle','-','color','b');

[ax2,h3,h4] = plotyy(dice_ustar.ygrid/1000.0, ...
                    mean(dplume_ustar.bmelt(:,:,end),1),...
                    dice_ustar.ygrid/1000.0, ...
                    mean(dice_ustar.lsurf(:,:,end),1));

set(h3,'Linewidth',lw,'LineStyle','-','color','k');
set(h4,'Linewidth',lw,'LineStyle','-','color','b');

[ax3,h5,h6] = plotyy(dice_sgd.ygrid/1000.0, ...
                    mean(dplume_sgd.bmelt(:,:,end),1),...
                    dice_sgd.ygrid/1000.0, ...
                    mean(dice_sgd.lsurf(:,:,end),1));
                

set(h5,'Linewidth',lw,'LineStyle','-','color','k');
set(h6,'Linewidth',lw,'LineStyle','-','color','b');

ymax = 24;
xmax = 62;
xmin = -2;


set(ax1(1),'xlim',[xmin xmax]);
set(ax1(2),'xlim',[xmin xmax]);
set(ax2(1),'xlim',[xmin xmax]);
set(ax2(2),'xlim',[xmin xmax]);
set(ax3(1),'xlim',[xmin xmax]);
set(ax3(2),'xlim',[xmin xmax]);

set(ax1(1),'ylim',[0 ymax]);
set(ax1(2),'ylim',[-600 0]);
set(ax2(1),'ylim',[0 ymax]);
set(ax2(2),'ylim',[-600 0]);
set(ax3(1),'ylim',[0 ymax]);
set(ax3(2),'ylim',[-600 0]);

set(ax1(1),'FontSize',fs);
set(ax1(2),'FontSize',fs);
set(ax2(1),'FontSize',fs);
set(ax2(2),'FontSize',fs);
set(ax3(1),'FontSize',fs);
set(ax3(2),'FontSize',fs);

set(ax1(1),'YTick',0:4:ymax);
set(ax1(2),'YTick',-600:100:0);
set(ax2(1),'YTick',0:4:ymax);
set(ax2(2),'YTick',-600:100:0);
set(ax3(1),'YTick',0:4:ymax);
set(ax3(2),'YTick',-600:100:0);

set(ax3(1),'YColor','k');
set(ax3(2),'YColor','b');

set(get(ax1(1),'Ylabel'),'String','basal melt rate (m/a)','FontSize',fs)
set(get(ax1(2),'Ylabel'),'String','ice draft (m)','FontSize',fs)

xlabel('Along shelf distance (km)','FontSize',fs);

%legend([h1 h3 h5 h2 h4 h6], ...
%       'basal melt (min. thickness)',...
%       'basal melt (u_* offset)', ...
%       'basal melt (sub-glacial discharge)', ...
%       'ice draft (min thickness)', ...
%       'ice draft (u_* offset)',...       
%       'ice draft (sub-glacial discharge)',...
%       'Location','NorthEast');

