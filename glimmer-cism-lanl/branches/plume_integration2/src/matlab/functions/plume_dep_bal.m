function [ pdep_res ] = plume_dep_bal( dplume )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

dx = dplume.x(2)-dplume.x(1);
dy = dplume.y(2)-dplume.y(1);

cen_grad_x = @(A) (A(3:end,:,:)-A(1:end-2,:,:))/(2*dx);
cen_grad_y = @(A) (A(:,3:end,:)-A(:,1:end-2,:))/(2*dy);
right_diff_x = @(A) (A(2:end,:,:)-A(1:end-1,:,:))/dx;
right_diff_y = @(A) (A(:,2:end,:)-A(:,1:end-1,:))/dx;

U=dplume.u(1:end-1,:,:);
su=dplume.su(1:end-1,:,:);
V=dplume.v(:,1:end-1,:);
sv=dplume.sv(:,1:end-1,:);
salt = dplume.salt;
temp = dplume.temp;
ambsurf = dplume.bpos - dplume.pdep;
depth = abs(dplume.draft);
rhop = dplume.rhop;
pdep = dplume.pdep;

pdep_res.x = (dplume.x - 1375)/1000;
pdep_res.y = (dplume.y - 875)/1000;

pdep_res.thk_flux_div = ...
    pad_edge(1,0, ...
     (0.5*(pdep(3:end,:,:)+pdep(2:end-1,:,:)).*su(2:end,:,:) - ...
      0.5*(pdep(2:end-1,:,:)+pdep(1:end-2,:,:)).*su(1:end-1,:,:)))/dx + ...
    pad_edge(0,1, ...
     (0.5*(pdep(:,3:end,:)+pdep(:,2:end-1,:)).*sv(:,2:end,:) - ...
      0.5*(pdep(:,2:end-1,:)+pdep(:,1:end-2,:)).*sv(:,1:end-1,:)))/dy;      
  
pdep_res.D_vel_div  = pdep .* (Dx_xgrid(su,dx) + Dy_ygrid(sv,dy));

pdep_adv_y_derv = right_diff_y(pdep) .* sv;
pdep_res.D_adv_y = pad_edge(0,1, 0.5*(pdep_adv_y_derv(:,2:end,:)+pdep_adv_y_derv(:,1:end-1,:)));

pdep_adv_x_derv = right_diff_x(pdep) .* su;
pdep_res.D_adv_x = pad_edge(1,0, 0.5*(pdep_adv_x_derv(2:end,:,:)+pdep_adv_x_derv(1:end-1,:,:)));

pdep_res.D_adv_derv = pdep_res.D_adv_y+pdep_res.D_adv_x;

pdep_res.train = dplume.train/ (3600*24*365.25);
pdep_res.melt = dplume.bmelt;

end

