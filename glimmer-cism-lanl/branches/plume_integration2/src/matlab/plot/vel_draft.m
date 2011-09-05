figure(1);
clf;

xs = 20:2:61;
ys = 1:2:82;
t = 51;

hold on;
contourf(dice6.Xgrid(xs,ys)/1000,dice6.Ygrid(xs,ys)/1000,...
         -dice6.lsurf(xs,ys,t),30,'EdgeColor','None');colorbar; 
caxis([0 600]);
scale = 1.0;
lw = 1.5;

su = dplume6.su(xs,ys,t);
sv = dplume6.sv(xs,ys,t);

su(2,2) = 0;
sv(2,2) = 0.5;
quiver(dice6.Xgrid(xs,ys)/1000,dice6.Ygrid(xs,ys)/1000,...
       su,sv,...
       scale,'k','LineWidth',lw);
   
fs = 16;
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);
set(gca,'FontSize',fs);