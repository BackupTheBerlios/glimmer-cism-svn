function [amb_temp,amb_salt] = nc_read_amb(fname)


amb = load(fname);
depth = amb(:,1);
temp = amb(:,2);
salt = amb(:,3);
density = amb(:,4);

amb_temp = @(z) interp1(depth,temp,z,'linear');
amb_salt = @(z) interp1(depth,salt,z,'linear');


end
