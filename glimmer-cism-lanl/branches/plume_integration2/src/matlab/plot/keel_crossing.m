f_ice6 = '/data/gladish/gc_output/igs_jobs/may29_highres_temp_-10.0_tau_25000.0_diff_10.0_k_4.0_amp_50.0_tempbot_-0.1_pmin_10.0/may29_highres_temp_-10.0_tau_25000.0_diff_10.0_k_4.0_amp_50.0_tempbot_-0.1_pmin_10.0.out.nc';
f_plume6='/data/gladish/gc_output/igs_jobs/may29_highres_temp_-10.0_tau_25000.0_diff_10.0_k_4.0_amp_50.0_tempbot_-0.1_pmin_10.0/plume.may29_highres_temp_-10.0_tau_25000.0_diff_10.0_k_4.0_amp_50.0_tempbot_-0.1_pmin_10.0.out.nc';

%dice6 = nc_ice_read(f_ice6,10,51);
%dplume6 = nc_plume_read(f_plume6,10,51);

[flat_ocean,flat_ice] = flatten_gc(dplume6,dice6,51);
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



