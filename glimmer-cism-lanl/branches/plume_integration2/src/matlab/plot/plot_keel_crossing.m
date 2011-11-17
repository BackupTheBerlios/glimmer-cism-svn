function plot_keel_crossing(dice,docean,times,fig_dir)

[flat_ocean,flat_ice] = flatten_gc(docean,dice,times);

fs = 14;
figure(1);
clf;


l = (flat_ocean.curvature >= 0.0);
l2 = (flat_ocean.curvature < 0.0);

xmin = min(flat_ocean.gradx(l));
xmax = max(flat_ocean.gradx(l));
ymin = min(flat_ocean.u(l));
ymax = max(flat_ocean.u(l));

xmin = min(flat_ocean.gradx);
xmax = max(flat_ocean.gradx);
ymin = min(flat_ocean.u);
ymax = max(flat_ocean.u);

fig1 = figure(1);
keelsize_x = 800;
keelsize_y = 650;
set(fig1,'Position',[1 1 keelsize_x keelsize_y]);
subplot(1,2,1);
hold on;
plot(flat_ocean.gradx(l2),flat_ocean.u(l2),'r.')
plot([xmin xmax], [0 0],'k')
plot([0 0], [ymin ymax],'k')
xlim([xmin xmax]);
ylim([ymin ymax]);
xlabel('cross-shelf draft gradient','FontSize',fs);
ylabel('cross-shelf plume transport (m^2/s)','FontSize',fs);
title('Transport in channel crests','fontsize',fs)
set(gca,'FontSize',fs);

subplot(1,2,2);
hold on
plot(flat_ocean.gradx(l),flat_ocean.u(l),'b.')
plot([xmin xmax], [0 0],'k')
plot([0 0], [ymin ymax],'k')
xlim([xmin xmax]);
ylim([ymin ymax]);
xlabel('cross-shelf draft gradient','FontSize',fs);
%ylabel('cross-shelf plume transport (m^2/s)','FontSize',fs);
title('Transport over channel keels','fontsize',fs)
set(gca,'FontSize',fs);

print('-depsc',strcat([fig_dir,'/plume_keel_crossing']));






end 
