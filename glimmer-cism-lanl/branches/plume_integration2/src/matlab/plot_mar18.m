rootdir = '/scratch/cvg222/gcp/GC_jobs/';
prefixes = {'mar15.1','mar15_koch_central'};
ks = {'2.0','4.0'};
amps = {'10.0','25.0'};
times = [1, -1];

addpath('/home/cvg222/code/pi2/src/matlab');

for i=1:2
  for j=1:2
    for k=1:2
      for l=1:2

  jobname = strcat([char(prefixes(i)),'_k_',char(ks(j)),'_amp_',char(amps(k)), ...
			    '_upvel_1000.0_temp_0.0'])
       plume_filename = strcat([rootdir,jobname,'/plume.',jobname,'.out.nc'])
	times(l)

	data = nc_plume_read(plume_filename,times(l));

plot_channels(data)
  fname = strcat(['/home/cvg222/matlab_figs/',jobname,'_time_',sprintf('%d',times(l)),'.channels.eps'])
print('-depsc',fname);
close

plot_bmlt_train_overlay(data)
  fname = strcat(['/home/cvg222/matlab_figs/',jobname,'_time_',sprintf('%d',times(l)),'.melt_train.eps'])
print('-depsc',fname);
close

 plot_velocity_overlay(data)
  fname = strcat(['/home/cvg222/matlab_figs/',jobname,'_time_',sprintf('%d',times(l)),'.velocity.eps'])
  print('-depsc',fname);
  close
  plot_interface_pdep_overlay(data)
    fname = strcat(['/home/cvg222/matlab_figs/',jobname,'_time_',sprintf('%d',times(l)),'.interface_pdep.eps'])
  print('-depsc',fname);
  close

if (i==2)
  plot_entr_overlay(data)
    fname = strcat(['/home/cvg222/matlab_figs/',jobname,'_time_',sprintf('%d',times(l)),'.entr.eps'])
  print('-depsc',fname);
  close
end

end
end
end
end
