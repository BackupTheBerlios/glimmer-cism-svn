function [ y_vel, y_thk,vvel_center,thk_center ] = centerline_profile(y0,y1,vvel,thck )

[m,n] = size(thck);

ice_front = 5;
end_buf = 2;
level = 1;

y_vel = y0((ice_front-1):(n-1-end_buf));
y_thk = y1(ice_front:(n-end_buf));
vvel_center = vvel(floor(m/2), (ice_front-1):(n-1-end_buf),  level);

thk_center = thck(floor(m/2), ice_front:(n-end_buf));

end

