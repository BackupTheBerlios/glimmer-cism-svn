function [] = plume_momentum_bal(dplume,amb_t,amb_s)

% evaluate the terms in the momentum equations on the scalar
% grid with a 2 cell border omitted around the edges


g = 9.81;
rho0 = 1030.0;
betaS = 7.86e-4;
betaT = 3.87e-5;
f = 2*(2*pi)/(3600.0*24.0)*sin(80 * pi/180);

cd = 2.5e-3;
diff = 10.0;
dx = dplume.x(2)-dplume.x(1);
dy = dplume.y(2)-dplume.y(1);

cen_div_func = @(flux_x,flux_y) ...
  pad_edge((flux_x(3:end,2:end-1,:)-flux_x(1:end-2,2:end-1,:)) / ...
	   (dplume.x(3)-dplume.x(1)) + ...
	   (flux_y(2:end-1,3:end,:)-flux_y(2:end-1,1:end-2,:)) / ...
	   (dplume.y(3)-dplume.y(1)));

cen_grad_x = @(A) (A(3:end,:,:)-A(1:end-2,:,:))/dx;
cen_grad_y = @(A) (A(:,3:end,:)-A(:,1:end-2,:))/dy;

u_trans_adv = zeros(size(dplume.u,1) ,size(dplume.u,2) ,size(dplume.u,3));
v_trans_adv = zeros(size(dplume.v,1) ,size(dplume.v,2) ,size(dplume.v,3));

u_diff = zeros(size(dplume.u,1),size(dplume.u,2),size(dplume.u,3));
v_diff = zeros(size(dplume.u,1),size(dplume.u,2),size(dplume.u,3));

den_grad_x = zeros(size(dplume.u,1),size(dplume.u,2),size(dplume.u,3));
den_grad_y = zeros(size(dplume.u,1),size(dplume.u,2),size(dplume.u,3));

isopyc_grad_x = zeros(size(dplume.u,1),size(dplume.u,2),size(dplume.u,3));
isopyc_grad_y = zeros(size(dplume.u,1),size(dplume.u,2),size(dplume.u,3));

drag_x = zeros(size(dplume.u,1),size(dplume.u,2),size(dplume.u,3));
drag_y = zeros(size(dplume.u,1),size(dplume.u,2),size(dplume.u,3));

cor_x = zeros(size(dplume.u,1),size(dplume.u,2),size(dplume.u,3));
cor_y = zeros(size(dplume.u,1),size(dplume.u,2),size(dplume.u,3));

U=dplume.u(1:end-1,:,:);
su=dplume.su(1:end-1,:,:);
V=dplume.v(:,1:end-1,:);
sv=dplume.sv(:,1:end-1,:);
salt = dplume.salt;
temp = dplume.temp;
ambsurf = dplume.bpos - dplume.pdep;
depth = abs(dplume.draft);
rhop = dplume.rhop(:,:,:);
pdep = dplume.pdep(:,:,:);



% D_x(U*su)+D_y(U*sy)
u_trans_adv = Dx_xgrid(U.*su)/dx + Dy_xgrid(U.*y_to_x(sv))/dy;

% D_x(V*su)+D_y(V*sy)
v_trans_adv = Dx_ygrid(V.*x_to_y(su))/dx + Dy_ygrid(V.*sv)/dy;

u_diff = pad_edge(1,0,cen_grad_x(diff*pdep.*Dx_xgrid(su)))/(2*dx*dx) + ...
         pad_edge(0,1,cen_grad_y(diff*pdep.*Dy_xgrid(su)))/(2*dy*dy);
v_diff = pad_edge(1,0,cen_grad_x(diff*pdep.*Dx_ygrid(sv)))/(2*dx*dx) + ...
         pad_edge(0,1,cen_grad_y(diff*pdep.*Dy_ygrid(sv)))/(2*dy*dy);

total_diff = sqrt(u_diff .^2 + v_diff.^2);
trans_adv =  u_trans_adv.^2+v_trans_adv.^2;

cor_x = -2*(2*pi/(3600.0*24.0))*sin(80.0*pi/180)* ...
  	       pad_edge(0,1,0.5*(V(:,2:end,:)+V(:,1:end-1,:)));
cor_y = 2*(2*pi/(3600.0*24.0))*sin(80.0*pi/180)* ...
               pad_edge(1,0,0.5*(U(2:end,:,:)+U(1:end-1,:,:)));

den_grad_x = (g/(2*rho0))*pad_edge(1,0,cen_grad_x(rhop)/dx);
den_grad_y = (g/(2*rho0))*pad_edge(0,1,cen_grad_y(rhop)/dy);

isopyc_grad_x = g*(betaS*(amb_s(ambsurf)-salt)-betaT*(amb_t(ambsurf)-temp)) .* ...
		  pdep.*pad_edge(1,0,cen_grad_x(ambsurf));
isopyc_grad_y = g*(betaS*(amb_s(ambsurf)-salt)-betaT*(amb_t(ambsurf)-temp)).* ...
		  pdep.*pad_edge(0,1,cen_grad_y(ambsurf));

drag_x = 0.5*cd*pad_edge(0,0,sqrt(sum_square_xy_grid(su,sv))) ... 
   	      .*pad_edge(1,0,U(2:end,:,:)+U(1:end-1,:,:));

drag_y = 0.5*cd*pad_edge(0,0,sqrt(sum_square_xy_grid(su,sv))) ...
	       .*pad_edge(0,1,V(:,2:end,:)+V(:,1:end-1,:));

drag = sqrt(drag_x.^2+drag_y.^2);

reynolds = trans_adv ./ total_diff;
rossby = trans_adv./sqrt(cor_x.^2+cor_y.^2);
speed = sum_square_xy_grid(su,sv);

figure(1);
clf;
subplot(2,2,1);
contourf(dplume.xgrid,dplume.ygrid,u_trans_adv,'EdgeColor','None');colorbar;
caxis([-1e-3 1e-3]);
subplot(2,2,2);
contourf(dplume.xgrid,dplume.ygrid,v_trans_adv,'EdgeColor','None');colorbar;
caxis([-0.01 0.01]);
subplot(2,2,3);
contourf(dplume.xgrid,dplume.ygrid,sqrt(cor_x.^2+cor_y.^2),'EdgeColor','None');colorbar;
subplot(2,2,4);
contourf(dplume.xgrid,dplume.ygrid,log(rossby),'EdgeColor','None');colorbar;

figure(2);
plot(flatten_field(rossby),flatten_field(dplume.bmelt(:,:,1)),'.');

figure(3);
subplot(2,2,1);
contourf(dplume.xgrid,dplume.ygrid,u_diff,'EdgeColor','None');colorbar;
subplot(2,2,2);
contourf(dplume.xgrid,dplume.ygrid,v_diff,'EdgeColor','None');colorbar;
subplot(2,2,3);
contourf(dplume.xgrid,dplume.ygrid,reynolds,'EdgeColor','None');colorbar;
subplot(2,2,4);
contourf(dplume.xgrid,dplume.ygrid,cd*(speed.^2)./trans_adv,'EdgeColor','None');colorbar;


end
