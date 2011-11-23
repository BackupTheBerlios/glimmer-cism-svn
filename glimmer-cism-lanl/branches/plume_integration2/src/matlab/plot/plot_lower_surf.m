%%% plot a 3d surface showing melt rates coloring the lower surface of the shelf


ncontours = 25;
fs = 28;

ec = 'None';
astrg = 0.6;
azimuth = 17.5;
azimuth = 40;

elevation = 50;
elevation = 30;
lpos = [0 0 -5];

thermal_forcing = @(depth,salt,temp) temp-(-5.73e-2*salt-7.61d-4*depth+8.32d-2);

tforce = thermal_forcing(dplume_avg.bpos-dplume_avg.pdep-800,dplume_avg.salt,dplume_avg.temp);

fig1 = figure(1);
clf;
set(fig1,'Position',[1 1 1200 800]);
%subplot(1,2,1);
%s = surf(dice.xgrid/1000,dice.ygrid/1000, ...
%	 double(dice_avg.lsurf)',double(tforce)');
%set(s,'FaceLighting','phong','FaceColor','interp','AmbientStrength',astrg);
%camlight right;
%%light('Position',[50 20 200],'style','infinite');

%colorbar;
%set(gca,'FontSize',fs);
%set(s,'EdgeColor',ec);
%xlabel('Across shelf distance (km)','FontSize',fs);
%ylabel('Along shelf distance (km)','FontSize',fs);
%zlabel('Ice draft (m)','FontSize',fs);
%title('Ice basal surface colored by thermal forcing (deg C)','FontSize',fs);

%view([azimuth elevation]);



%subplot(1,2,2);
s = surf(dice.xgrid/1000,dice.ygrid/1000, ...
         double(dice_avg.lsurf)',double(dice_avg.bmlt)');
set(s,'FaceLighting','phong','FaceColor','interp','AmbientStrength',astrg);
camlight right;

colorbar;
set(gca,'FontSize',fs);
set(s,'EdgeColor',ec);
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);
zlabel('Ice draft (m)','FontSize',fs);
%title('Ice basal surface colored by melt rate (m/a)','FontSize',fs);
%azimuth = -135;
%elevation = -30;

%azimuth = 45;
%elevation = 20;
view([azimuth elevation]);

%subplot(1,2,2);
%hold on
%draft = squeeze(dplume.draft(:,:,end));
%su = squeeze(dplume.su(:,:,end));
%sv = squeeze(dplume.sv(:,:,end));
%contourf(x,y,-draft,40,'EdgeColor','None');
%colorbar('FontSize',fs);
%set(gca,'FontSize',fs);
%%caxis([000 650])
%scale=2.0;stride=2;
%i0 = 2;
%quiver(x(i0:stride:end)/1000,y(1:stride:end)/1000,...
%      su(i0:stride:end,1:stride:end),...
%      sv(i0:stride:end,1:stride:end),...
%      scale,'k','LineWidth',lw);
%colormap jet;
%xlabel('Across shelf distance (km)','FontSize',fs);
%ylabel('Along shelf distance (km)','FontSize',fs);
%caxis([0 600]);
%title('Contours of Ice Draft (m) with plume velocities','FontSize',fs);
%hold off
%print('-depsc',strcat([fig_dir,'/plume_draft_vel']));
