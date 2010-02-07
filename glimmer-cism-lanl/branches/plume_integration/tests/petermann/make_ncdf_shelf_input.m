function make_ncdf_shelf(nx,ny,dx,dy, x_wavenumber, amplitude, gl_depth, if_depth, out_file_name)

% gl_depth and if_depth must be negative numbers

ncid = netcdf.create(out_file_name,'NC_CLOBBER');

x1_id = netcdf.defDim(ncid, 'x1', nx);
y1_id = netcdf.defDim(ncid, 'y1', ny);
time_id = netcdf.defDim(ncid, 'time', 1);

netcdf.defVar(ncid, 'x1','double',[x1_id]);
netcdf.defVar(ncid, 'y1', 'double', [y1_id]);
netcdf.defVar(ncid,'time','int',[time_id]);

thk_id = netcdf.defVar(ncid, 'thk', 'double', [x1_id,y1_id,time_id]);
kinbcmask_id = netcdf.defVar(ncid,'kinbcmask', 'int', [x1_id,y1_id,time_id]);
topg_id = netcdf.defVar(ncid,'topg','double',[x1_id,y1_id]);
%lsrf_id = netcdf.defVar(ncid,'lsrf','double', [x1_id,y1_id,time_id]);

netcdf.endDef(ncid);

xs = (0:(nx-1))*dx;
ys = (0:(ny-1))*dy;

netcdf.putVar(ncid, x1_id, xs);
netcdf.putVar(ncid, y1_id, ys);
netcdf.putVar(ncid, time_id, [1]);

% NB: Xs and Ys looks backwards in the following line,
% but it is correct because we want Xs(i,j) to increase 
% with increasing i and we want Ys(i,j) to increase with 
% increasing j.
[Ys,Xs] = meshgrid(xs,ys);

% Gaussian bump
%xmid = dx*nx / 2;
%ymid = dy*ny / 2;
%thk = exp(10*(-(Xs-xmid).^2 - (Ys-ymid).^2)/(dx*nx+dy*ny)^2);

% makes linear rising shelf with longitudinal channels
lsrf = gl_depth + (if_depth-gl_depth)*(Ys/(ys(ny)));
lsrf = lsrf + 0.5*amplitude*(1 - cos((Xs/xs(nx))*2*pi*x_wavenumber));
thk = 1030/920*abs(lsrf);

thk(1:2,:) = 0.0;
thk((nx-1):nx,:) = 0.0;
thk(:,1:2) = 0.0;
% experimenting with zeroing out the north edge to create 
% a calving front
thk(:,(ny-3):ny) = 0.0;

kinbcmask = round(zeros(nx,ny));
kinbcmask(1:3, :) = 1;
kinbcmask((nx-2):nx,:) = 1;
kinbcmask(:,1:3) = 1;

topg = zeros(nx,ny);
topg(1:nx,1:ny) = gl_depth ;
topg(1:2,:) = 0.0; % sea level
topg((nx-1):nx,:) = 0.0;
topg(:,1:2) = 0.0;

%contourf(Xs,Ys,thk);colorbar;
%contourf(thk);colorbar;

netcdf.putVar(ncid, thk_id,  thk );
netcdf.putVar(ncid, kinbcmask_id, kinbcmask);
netcdf.putVar(ncid, topg_id, topg);

% NOTE glimmer doesn't want to read in lsrf
% netcdf.putVar(ncid, lsrf_id, lsrf);

netcdf.close(ncid)

end
