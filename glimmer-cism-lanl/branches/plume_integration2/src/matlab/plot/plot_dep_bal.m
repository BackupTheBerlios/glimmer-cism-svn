figure(4);
clf;

c0 = -2.5e-2;
c1 = 2.5e-2;

subplot(2,2,1);
contourf(x,y,plume_res.D_vel_div','EdgeColor','None');colorbar;
caxis([c0 c1]);

subplot(2,2,2);
contourf(x,y,plume_res.D_adv_x'+plume_res.D_adv_y'+plume_res.D_vel_div','EdgeColor','None');colorbar;
caxis([c0 c1]);

subplot(2,2,3);
contourf(x,y,-plume_res.thk_flux_div','EdgeColor','None');colorbar;
caxis([c0 c1]);

subplot(2,2,4);
contourf(x,y,plume_res.train','EdgeColor','None');colorbar;
caxis([c0 c1]);


figure(5);
clf;
hold on
%plot(flatten_field(plume_res.train),flatten_field(plume_res.train),'k.');

%plot(flatten_field(plume_res.train), ...
%     flatten_field(plume_res.D_adv_x),'r.');
xmin = min(flatten_field(plume_res.train));
xmax = max(flatten_field(plume_res.train));
ymin = min(flatten_field(plume_res.D_vel_div));
ymax = max(flatten_field(plume_res.D_adv_x+plume_res.D_adv_y));

plot([xmin xmax],[xmin xmax],'k-','linewidth',2.0);

ms = 16.0;

plot(flatten_field(plume_res.train), ...
     flatten_field(plume_res.D_adv_x+plume_res.D_adv_y),'b.','MarkerSize',ms);
plot(flatten_field(plume_res.train), ...
     flatten_field(plume_res.D_vel_div),'r.','MarkerSize',ms);
plot(flatten_field(plume_res.train), ...
     flatten_field(plume_res.thk_flux_div),'m.','MarkerSize',ms);
 plot([xmin,xmax],[0, 0],'k-');       
 plot([0,0],[ymin, ymax],'k-');       
xlim([-0.025 0.004]);

fs = 14;
title('Plume thickness equation terms','FontSize',fs)
xlabel('entrainment/detrainment rates (m/s)','FontSize',fs);
ylabel('rate of thickness increase (m/s)','FontSize',fs);
set(gca,'FontSize',fs);
obj = legend('$y=x$','$\vec{u}\cdot\nabla D$', ...
       '$D\nabla\cdot\vec{u}$','$\nabla\cdot\left(D\vec{u}\right)$',...
       'Location','Northwest'); %,...
set(obj,'Interpreter','latex');
   % 'Interpreter','latex');