function [] = plume_momentum_bal(dplume,amb_t,amb_s)

% evaluate the terms in the momentum equations on the scalar
% grid with a 2 cell border omitted around the edges

  fs = 14;

g = 9.81;
rho0 = 1030.0;
betaS = 7.86e-4;
betaT = 3.87e-5;
f = 2*(2*pi)/(3600.0*24.0)*sin(80 * pi/180);

cd = 2.5e-3;
diff = 10.0;
dx = dplume.x(2)-dplume.x(1);
dy = dplume.y(2)-dplume.y(1);

%cen_div_func = @(flux_x,flux_y) ...
%  pad_edge((flux_x(3:end,2:end-1,:)-flux_x(1:end-2,2:end-1,:)) / ...
%	   (2*dx) + ...
%	   (flux_y(2:end-1,3:end,:)-flux_y(2:end-1,1:end-2,:)) / ...
%	   (2*dy);

cen_grad_x = @(A) (A(3:end,:,:)-A(1:end-2,:,:));
cen_grad_y = @(A) (A(:,3:end,:)-A(:,1:end-2,:));

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


vel_div = Dx_xgrid(su)/dx + Dy_ygrid(sv)/dy;
flux_div = Dx_xgrid(U)/dx + Dy_ygrid(V)/dy;

delta_w = vel_div .* pdep;

 trim_edge = @(M) M(4:end-3,4:end-3,:);

% D_x(U*su)+D_y(U*sy)
u_trans_adv = Dx_xgrid(U.*su)/dx + Dy_xgrid(U.*y_to_x(sv))/dy;

% D_x(V*su)+D_y(V*sy)
v_trans_adv = Dx_ygrid(V.*x_to_y(su))/dx + Dy_ygrid(V.*sv)/dy;

u_diff = pad_edge(1,0,cen_grad_x(diff*pdep.*Dx_xgrid(su)))/(2*dx*dx) + ...
         pad_edge(0,1,cen_grad_y(diff*pdep.*Dy_xgrid(su)))/(2*dy*dy);
v_diff = pad_edge(1,0,cen_grad_x(diff*pdep.*Dx_ygrid(sv)))/(2*dx*dx) + ...
         pad_edge(0,1,cen_grad_y(diff*pdep.*Dy_ygrid(sv)))/(2*dy*dy);


grad_pdep = sqrt(pad_edge(1,0,cen_grad_x(pdep)/(2*dx)).^2 + ...
                 pad_edge(0,1,cen_grad_y(pdep)/(2*dy)).^2);


adv_derv_pdep = pad_edge(1,0, ...
                         0.5*(su(2:end,:,:)+su(1:end-1,:,:)).* cen_grad_x(pdep)/(2*dx) ...
                ) + ...
                pad_edge(0,1, ...
                         0.5*(sv(:,2:end,:)+sv(:,1:end-1,:)).* cen_grad_y(pdep)/(2*dy) ...
                );


total_diff = sqrt(u_diff .^2 + v_diff.^2);
trans_adv =  u_trans_adv.^2+v_trans_adv.^2;

cor_x = -2*(2*pi/(3600.0*24.0))*sin(80.0*pi/180)* ...
  	       pad_edge(0,1,0.5*(V(:,2:end,:)+V(:,1:end-1,:)));
cor_y = 2*(2*pi/(3600.0*24.0))*sin(80.0*pi/180)* ...
               pad_edge(1,0,0.5*(U(2:end,:,:)+U(1:end-1,:,:)));

den_grad_x = (g/(2*rho0))*pad_edge(1,0,cen_grad_x(rhop)/(2*dx));
den_grad_y = (g/(2*rho0))*pad_edge(0,1,cen_grad_y(rhop)/(2*dy));

isopyc_grad_x = g*(betaS*(amb_s(ambsurf)-salt)-betaT*(amb_t(ambsurf)-temp)) .* ...
  pdep.*pad_edge(1,0,cen_grad_x(ambsurf)/(2*dx));
isopyc_grad_y = g*(betaS*(amb_s(ambsurf)-salt)-betaT*(amb_t(ambsurf)-temp)).* ...
  pdep.*pad_edge(0,1,cen_grad_y(ambsurf)/(2*dy));



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

x = dplume.xgrid/1000.0;
y = dplume.ygrid/1000.0;

contourf(x,y,u_trans_adv,'EdgeColor','None');colorbar;
%caxis([-1e-3 1e-3]);
title('u convection term','FontSize',fs);

subplot(2,2,2);
contourf(x,y,v_trans_adv,'EdgeColor','None');colorbar;
%caxis([-0.01 0.01]);
title('v convection term','FontSize',fs)

subplot(2,2,3);
contourf(x,y,sqrt(cor_x.^2+cor_y.^2),'EdgeColor','None');colorbar;
title('coriolis accel','FontSize',fs);

subplot(2,2,4);
contourf(x,y,rossby,'EdgeColor','None');colorbar;
title('Rossby Number','FontSize',fs);

figure(2);
plot(flatten_field(rossby),flatten_field(dplume.bmelt(:,:,1)),'.');

figure(3);
subplot(2,2,1);
contourf(x,y,u_diff,'EdgeColor','None');colorbar;
title('u diffusion term','FontSize',fs);

subplot(2,2,2);
contourf(x,y,v_diff,'EdgeColor','None');colorbar;
title('v diffusion term','FontSize',fs);

subplot(2,2,3);
contourf(x,y,log(reynolds),'EdgeColor','None');colorbar;
title('log(Reynolds number)','FontSize',fs);

subplot(2,2,4);
contourf(x,y,cd*(speed.^2),'EdgeColor','None');
colorbar;
%./trans_adv,'EdgeColor','None');colorbar;
title('non-dimensional drag','FontSize',fs);


figure(4);
clf;
%contourf(x,y,delta_w,'EdgeColor','None');
%colorbar;
subplot(2,2,1);
hold on
plot(trim_edge(delta_w),trim_edge(dplume.train)/(3600.0*24.0*365.25),'k.');
plot(delta_w,delta_w,'b.');
xlabel('delta w (m/s)','FontSize',fs);
ylabel('trainment rate (m/s)','FontSize',fs);
hold off
title('Velocity divergence in plume','FontSize',fs);

subplot(2,2,2);
hold on
plot(trim_edge(adv_derv_pdep),trim_edge(dplume.train)/(3600.0*24.0*365.25),'k.');
plot(adv_derv_pdep,adv_derv_pdep,'b.');
xlabel('adv derivative of pdep (m/s)','FontSize',fs);
ylabel('trainment rate (m/s)','FontSize',fs);
hold off
title('Advective derivative of pdep','FontSize',fs);

subplot(2,2,3);
hold on
plot(trim_edge(flux_div),trim_edge(dplume.train)/(3600.0*24.0*365.25),'k.');
plot(flux_div,flux_div,'b.');
xlabel('flux div (m/s)','FontSize',fs);
ylabel('trainment rate (m/s)','FontSize',fs);
hold off
title('Total flux divergnece','FontSize',fs);

subplot(2,2,4);
hold on
bmlt =  dplume.bmelt/(3600.0*24*365.25);
%plot(bmlt,dplume.train/(3600.0*24.0*365.25),'k.');
%plot(bmlt,bmlt,'b.');
plot(trim_edge(flux_div),trim_edge(delta_w+adv_derv_pdep),'r.');
plot(flux_div,flux_div,'b.');
xlabel('flux_div (m/s)','FontSize',fs);
ylabel('sum of two terms (m/s)','FontSize',fs);
hold off
title('bmlt','FontSize',fs);


%figure(6)
%  contourf(grad_pdep);colorbar;
% contourf(x,y,adv_derv_pdep);colorbar;

end
