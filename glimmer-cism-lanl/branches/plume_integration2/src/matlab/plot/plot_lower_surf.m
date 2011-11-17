%%% plot a 3d surface showing melt rates coloring the lower surface of the shelf


ncontours = 25;
fs = 14;

ec = 'None';
astrg = 0.6;
azimuth = 17.5;
azimuth = 30;

elevation = 50;
elevation = 30;
lpos = [0 0 -5];

thermal_forcing = @(depth,salt,temp) temp-(-5.73e-2*salt-7.61d-4*depth+8.32d-2);

tforce = thermal_forcing(dplume_avg.bpos-dplume_avg.pdep-800,dplume_avg.salt,dplume_avg.temp);

fig1 = figure(1);
clf;
set(fig1,'Position',[1 1 1600 800]);
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
title('Ice basal surface colored by melt rate (m/a)','FontSize',fs);
%azimuth = -135;
%elevation = -30;

%azimuth = 45;
%elevation = 20;
view([azimuth elevation]);

