figure(4);
clf;

c0 = -2.5e-2;
c1 = 2.5e-2;

fs = 14;

x = plume_res.x;
y = plume_res.y;
[x,y] = meshgrid(x,y);

subset = @(M) M(1:80,1:80);

subplot(2,2,1);
contourf(subset(x),subset(y),subset(plume_res.D_vel_div)','EdgeColor','None');colorbar;
caxis([c0 c1]);
title('Column stretching (m/s)','FontSize',fs);

subplot(2,2,2);
contourf(subset(x),subset(y), ...
 subset(plume_res.D_adv_x'+plume_res.D_adv_y'+plume_res.D_vel_div'), ...
 'EdgeColor','None');colorbar;
caxis([c0 c1]);

subplot(2,2,3);
contourf(subset(x),subset(y),subset(-plume_res.thk_flux_div)','EdgeColor','None');colorbar;
caxis([c0 c1]);

subplot(2,2,4);
contourf(subset(x),subset(y),subset(plume_res.train)','EdgeColor','None');colorbar;
caxis([c0 c1]);


figure(5);
clf;
hold on
%plot(flatten_field(plume_res.train),flatten_field(plume_res.train),'k.');

%plot(flatten_field(plume_res.train), ...
%     flatten_field(plume_res.D_adv_x),'r.');
xmin = min(flatten_field(plume_res.train));
xmax = max(flatten_field(plume_res.train))+0.0025;
ymin = min(flatten_field(plume_res.D_vel_div));
ymax = max(flatten_field(plume_res.D_adv_x+plume_res.D_adv_y))+0.0025;



ms = 16.0;
lw = 2.5;
plot([xmin,xmax],[0, 0],'k-','LineWidth',lw);       
plot([0,0],[ymin, ymax],'k-','LineWidth',lw);       
plot([xmin xmax],[xmin xmax],'k-','linewidth',lw);

h1 = plot(flatten_field(plume_res.train), ...
     flatten_field(plume_res.D_adv_x+plume_res.D_adv_y),'b.','MarkerSize',ms);
h2 = plot(flatten_field(plume_res.train), ...
     flatten_field(plume_res.D_vel_div),'r.','MarkerSize',ms);
h3 = plot(flatten_field(plume_res.train), ...
     flatten_field(plume_res.thk_flux_div),'g.','MarkerSize',ms);

xlim([-0.025 0.004]);

fs = 14;
legend([h1 h2 h3],'advection of thickness',...
       'column stretching', ...
       'thickness flux divergence','Location','Northwest');

title('Plume thickness equation terms','FontSize',fs)
xlabel('entrainment/detrainment rates (m/s)','FontSize',fs);
ylabel('rate of thickness increase (m/s)','FontSize',fs);
set(gca,'FontSize',fs);
%obj = legend('$y=x$','$\vec{u}\cdot\nabla D$', ...
%       '$D\nabla\cdot\vec{u}$','$\nabla\cdot\left(D\vec{u}\right)$',...
%       'Location','Northwest'); 
%set(obj,'Interpreter','latex');
%   % 'Interpreter','latex');

