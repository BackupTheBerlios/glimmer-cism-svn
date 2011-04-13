function [] = plot_bmlt_depth(nc_filename,timeslice)

  [x,y,su,sv,u,v,bmelt,bpos] = nc_plume_read(nc_filename,timeslice);
  figure('Units','centimeters','Position',[0 0 10 20]);

  subplot(2,1,1);
  contourf(x/1000.0,y/1000.0,bpos'-650.0,20);colorbar;
  subplot(2,1,2); 
  contourf(x/1000.0,y/1000.0,bmelt',10,'k') ;colorbar;

end 
