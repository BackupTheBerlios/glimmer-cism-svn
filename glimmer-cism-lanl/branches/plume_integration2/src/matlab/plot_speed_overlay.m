function [] = plot_bmlt_depth(nc_filename,timeslice)

  data = nc_plume_read(nc_filename,timeslice);
  figure('Units','centimeters','Position',[0 0 30 20]);

end 
