function [] = plot_ice_velocity_overlay(data)

fs = 14;
%figure('Units','centimeters','Position',[0 0 40 20]);
x = data.Xgrid;
y = data.Ygrid;
[m,n] = size(x);
u = squeeze(data.uvelmean(1:m,1:n,end));
v = squeeze(data.vvelmean(1:m,1:n,end));

draft = -data.lsurf(:,:,end);


stride = 3;
stride2 = 5;
istart = 2;
jstart = 1;
scale_vel = 0.025;
scale = 0;
lw = 1.0;

figure(1);
clf;


hold on
contourf(x/1000.0,y/1000.0,draft,40,'EdgeColor','None');colorbar;
caxis([0 600])
quiver(x(istart:stride:end,jstart:stride2:end)/1000, ...
       y(istart:stride:end,jstart:stride2:end)/1000,...
       scale_vel*u(istart:stride:end,jstart:stride2:end),...
       scale_vel*(v(istart:stride:end,jstart:stride2:end)-1000.0),...
      scale,'k','LineWidth',lw,'autoscale','off');
colormap jet;
rectangle('Position',[15 36.5 4.5 2.5],'Facecolor','w');
h = quiver(19,37,0,scale_vel*50,scale,'k','LineWidth',lw,'autoscale','off');
adjust_quiver_arrowhead_size(h,4);
text(15.5,37.5,'50 m/year','FontSize',12);

xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);

hold off

set(gca,'FontSize',fs);
end

