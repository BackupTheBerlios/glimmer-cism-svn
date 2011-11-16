function [res] = plume_momentum_bal2(dplume,amb_t,amb_s)

% evaluate the terms in the momentum equations on the scalar
% grid with a 2 cell border omitted around the edges

g = 9.81;
rho0 = 1030.0;
betaS = 7.86e-4;
betaT = 3.87e-5;
lat = 80.0;
f = 2*(2*pi)/(3600.0*24.0)*sin(lat * pi/180);

cd = 2.5e-3;
diff = 10.0;

dx = dplume.x(2)-dplume.x(1);
dy = dplume.y(2)-dplume.y(1);

cen_grad_x = @(A) (A(3:end,:,:)-A(1:end-2,:,:))/(2*dx);
cen_grad_y = @(A) (A(:,3:end,:)-A(:,1:end-2,:))/(2*dy);

U=dplume.u(1:end-1,:,:);
su=dplume.su(1:end-1,:,:);
V=dplume.v(:,1:end-1,:);
sv=dplume.sv(:,1:end-1,:);
salt = dplume.salt;
temp = dplume.temp;
ambsurf = dplume.bpos - dplume.pdep;
rhop = dplume.rhop;
pdep = dplume.pdep;
bmlt = dplume.bmelt/(3600.0*24.0*365.25);
train = dplume.train/(3600.0*24.0*365.25);

res.x = (dplume.x-1375)/1000;
res.y = (dplume.y - 875)/1000;

res.su_adv_derv = 0.5*(pad_edge(1,0,su(2:end,:,:)+su(1:end-1,:,:)) ) ...
                 .*Dx_xgrid(su,dx) + ...
                 0.5*(pad_edge(0,1,sv(:,2:end,:)+sv(:,1:end-1,:))) ...
                 .*Dy_xgrid(su,dy);
res.sv_adv_derv = 0.5*(pad_edge(1,0,su(2:end,:,:)+su(1:end-1,:,:)) ) ...
                 .*Dx_ygrid(sv,dx) + ...
                 0.5*(pad_edge(0,1,sv(:,2:end,:)+sv(:,1:end-1,:))) ...
                 .*Dy_ygrid(sv,dy);

res.u_boundary = -(1.0*cd*sqrt(sum_square_xy_grid(su,sv))+1.0*bmlt+1.0*0.5*(train+abs(train))).* ...
  pad_edge(1,0,0.5*(su(2:end  ,:,:).^2./U(2:end  ,:,:) + ...
                    su(1:end-1,:,:).^2./U(1:end-1,:,:)));

res.v_boundary = -(1.0*cd*sqrt(sum_square_xy_grid(su,sv))+1.0*bmlt+1.0*0.5*(train+abs(train))).* ...
  pad_edge(0,1,0.5*(sv(:,2:end  ,:).^2./V(:,2:end  ,:) + ...
                    sv(:,1:end-1,:).^2./V(:,1:end-1,:)));
              
res.u_diff = pad_edge(1,0,cen_grad_x(diff*Dx_xgrid(su,dx))) + ...
             pad_edge(0,1,cen_grad_y(diff*Dy_xgrid(su,dy)));
res.v_diff = pad_edge(1,0,cen_grad_x(diff*Dx_ygrid(sv,dx))) + ...
             pad_edge(0,1,cen_grad_y(diff*Dy_ygrid(sv,dy)));


res.cor_x = -f*pad_edge(0,1,0.5*(sv(:,    2:end,:)+sv(:,      1:end-1,:)));
res.cor_y =  f*pad_edge(1,0,0.5*(su(2:end,:,    :)+su(1:end-1,:,      :)));

res.den_grad_x = (g/(2*rho0))*pdep.*pad_edge(1,0,cen_grad_x(rhop));
res.den_grad_y = (g/(2*rho0))*pdep.*pad_edge(0,1,cen_grad_y(rhop));

res.g_prime = g*(betaS*(amb_s(ambsurf)-salt)-betaT*(amb_t(ambsurf)-temp));

res.isopyc_grad_x = res.g_prime .* pad_edge(1,0,cen_grad_x(ambsurf));
res.isopyc_grad_y = res.g_prime .* pad_edge(0,1,cen_grad_y(ambsurf));

res.isopyc_grad = sqrt(res.isopyc_grad_x.^2+res.isopyc_grad_y.^2);


res.u_diff2 = (diff) .* (Dx_xgrid(su,dx).*pad_edge(1,0,cen_grad_x(log(pdep))) + ...
			 Dy_xgrid(su,dy).*pad_edge(0,1,cen_grad_y(log(pdep))));
res.v_diff2 = (diff./pdep) .* (Dx_ygrid(sv,dx).*pad_edge(1,0,cen_grad_x(log(pdep))) + ...
			       Dy_ygrid(sv,dy).*pad_edge(0,1,cen_grad_y(log(pdep))));

end
