amps = {'25.0','50.0'};
ks = {'0.0','1.0','2.0','3.0', ...
      '4.0','5.0','6.0','7.0', ...
      '8.0','9.0','10.0','11.0', ...
      '12.0','13.0','14.0','15.0'};

dice = struct()

for i=1:length(amps)
  for j= 1:length(ks)
 
 datafile = ['/archive/cvg222/gc_output/2011/paper_jobs/',...
 sprintf('oct30_perturb_usq.9_thk_perturb_k_%s_amp_%s/',ks{j},amps{i}),...
 sprintf('oct30_perturb_usq.9_thk_perturb_k_%s_amp_%s.out.nc',ks{j},amps{i})]
 figname = sprintf('spectrum_k_%s_amp_%s.eps',ks{j},amps{i});
jobname = sprintf('k_%i_amp_%i',strread(ks{j}),strread(amps{i}))
dice.(jobname) = nc_ice_read(datafile,-1,1,-1);
plot_spectral_evolution(dice.(jobname),figname,ks{j},amps{i});
end
end
