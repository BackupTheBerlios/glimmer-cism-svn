%%% plot a 3d surface showing melt rates coloring the lower surface of the shelf

ncontours = 25;
fs = 20;

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
set(fig1,'Position',[1 1 800 500]);


s = surf(dice.xgrid/1000,dice.ygrid/1000, ...
         double(dice_avg.lsurf)',double(dice_avg.bmlt)');
set(s,'FaceLighting','phong','FaceColor','interp','AmbientStrength',astrg);
camlight right;

colorbar('FontSize',fs);
set(gca,'FontSize',fs);
set(s,'EdgeColor',ec);
xlabel('Across shelf distance (km)','FontSize',fs);
xlabel('x (km)','FontSize',fs);
xlim([0 20]);
ylim([-2 45]);
ylabel('Along shelf distance (km)','FontSize',fs);
ylabel('y (km)','FontSize',fs);
zlabel('Ice draft (m)','FontSize',fs);

azimuths = 0:2:360;
for i=1:length(azimuths)
  az  = azimuth + azimuths(i);
  view([az elevation]);

set(gcf,'PaperPositionMode','auto');
print('-depsc',sprintf('/scratch/cvg222/lsurf_frame_%03d.eps',i));
end
