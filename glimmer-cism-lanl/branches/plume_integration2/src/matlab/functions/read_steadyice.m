function [y_thk, thk, y_vel, vvel] = read_steadyice(ncfilename)

nc = ncload(ncfilename);

ntimes = nc.dims.time.length;
nx0 = nc.dims.x0.length;
nx1 = nc.dims.x1.length;
ny0 = nc.dims.y0.length;
ny1 = nc.dims.y1.length;

times = ncgetdata(nc, 'time');
x0 = ncgetdata(nc,'x0');
x1 = ncgetdata(nc,'x1');
y_vel = ncgetdata(nc,'y0');
y_thk = ncgetdata(nc,'y1');

thk = ncgetdata(nc,'thk');
vvel = ncgetdata(nc,'vvelhom');

thk_center = 0.5 * (thk(floor(nx1/2),:) + thk(floor(nx1/2)+1,:));
vel_center = 0.5 * (vvel(floor(nx0/2),:) + vvel(floor(nx0/2)+1,:));

%now extract the section corresponding to the final time
thk = thk_center((ny1*ntimes - (ny1-1)):(ny1*ntimes));
vvel = vel_center((ny0*ntimes - (ny0-1)):(ny0*ntimes));

end 

