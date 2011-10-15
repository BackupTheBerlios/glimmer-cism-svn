function plot_keel_crossing(dice,docean,tlen)

[flat_ocean,flat_ice] = flatten_gc(docean,dice,tlen);

figure(1);
clf;

hold on;

l = (flat_ocean.curvature >= 0.0);

plot(flat_ocean.gradx(l),flat_ocean.u(l),'b.')
xmin = min(flat_ocean.gradx(l));
xmax = max(flat_ocean.gradx(l));
ymin = min(flat_ocean.u(l));
ymax = max(flat_ocean.u(l));

plot([xmin xmax], [0 0],'k')
plot([0 0], [ymin ymax],'k')
xlim([xmin xmax]);
ylim([ymin ymax]);

fs = 16;
xlabel('cross-shelf draft gradient','FontSize',fs);
ylabel('cross-shelf plume transport (m^2/s)','FontSize',fs);
set(gca,'FontSize',fs);



end 
