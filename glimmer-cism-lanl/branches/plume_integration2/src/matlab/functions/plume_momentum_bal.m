function [res] = plume_momentum_bal(dplume,amb_t,amb_s)

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
depth = abs(dplume.draft);
rhop = dplume.rhop;
pdep = dplume.pdep;

res.vel_div  = Dx_xgrid(su,dx) + Dy_ygrid(sv,dy);
res.flux_div = Dx_xgrid(U, dx) + Dy_ygrid(V ,dy);

res.delta_w = res.vel_div .* pdep;

res.u_flux_div = Dx_xgrid(U.*su,dx)         + Dy_xgrid(U.*y_to_x(sv),dy);
res.v_flux_div = Dx_ygrid(V.*x_to_y(su),dx) + Dy_ygrid(V.*sv,dy);

res.u_adv_derv = 0.5*(pad_edge(1,0,su(2:end,:,:)+su(1:end-1,:,:)) ) ...
                 .*Dx_xgrid(U,dx) + ...
                 0.5*(pad_edge(0,1,sv(:,2:end,:)+sv(:,1:end-1,:))) ...
                 .*Dy_xgrid(U,dy);
res.v_adv_derv = 0.5*(pad_edge(1,0,su(2:end,:,:)+su(1:end-1,:,:)) ) ...
                 .*Dx_ygrid(V,dx) + ...
                 0.5*(pad_edge(0,1,sv(:,2:end,:)+sv(:,1:end-1,:))) ...
                 .*Dy_ygrid(V,dy);


res.u_diff = pad_edge(1,0,cen_grad_x(diff*pdep.*Dx_xgrid(su,dx))) + ...
             pad_edge(0,1,cen_grad_y(diff*pdep.*Dy_xgrid(su,dy)));
res.v_diff = pad_edge(1,0,cen_grad_x(diff*pdep.*Dx_ygrid(sv,dx))) + ...
             pad_edge(0,1,cen_grad_y(diff*pdep.*Dy_ygrid(sv,dy)));


res.grad_pdep = sqrt(pad_edge(1,0,cen_grad_x(pdep).^2) + ...
                     pad_edge(0,1,cen_grad_y(pdep).^2));

res.adv_derv_pdep = ...
pad_edge(1,0,0.5*(su(2:end,:,:)+su(1:end-1,:,:)).* cen_grad_x(pdep)) + ...
pad_edge(0,1,0.5*(sv(:,2:end,:)+sv(:,1:end-1,:)).* cen_grad_y(pdep));

res.total_diff     =  sqrt(res.u_diff .^2 + res.v_diff.^2);
res.trans_flux_div =  sqrt(res.u_flux_div.^2+res.v_flux_div.^2);

res.cor_x = -f*pad_edge(0,1,0.5*(V(:,2:end,:)+V(:,1:end-1,:)));
res.cor_y =  f*pad_edge(1,0,0.5*(U(2:end,:,:)+U(1:end-1,:,:)));

res.den_grad_x = (g/(2*rho0))*pdep.^2.*pad_edge(1,0,cen_grad_x(rhop));
res.den_grad_y = (g/(2*rho0))*pdep.^2.*pad_edge(0,1,cen_grad_y(rhop));

res.isopyc_grad_x = g*(betaS*(amb_s(ambsurf)-salt)-betaT*(amb_t(ambsurf)-temp)) .* ...
                    pdep.*pad_edge(1,0,cen_grad_x(ambsurf));
res.isopyc_grad_y = g*(betaS*(amb_s(ambsurf)-salt)-betaT*(amb_t(ambsurf)-temp)).* ...
                     pdep.*pad_edge(0,1,cen_grad_y(ambsurf));

res.isopyc_grad = sqrt(res.isopyc_grad_x.^2+res.isopyc_grad_y.^2);

res.drag_x = -0.5*cd*pad_edge(0,0,sqrt(sum_square_xy_grid(su,sv))) ... 
   	          .*pad_edge(1,0,su(2:end,:,:)+su(1:end-1,:,:));

res.drag_y = -0.5*cd*pad_edge(0,0,sqrt(sum_square_xy_grid(su,sv))) ...
   	          .*pad_edge(0,1,sv(:,2:end,:)+sv(:,1:end-1,:));

res.drag = sqrt(res.drag_x.^2+res.drag_y.^2);

res.reynolds = res.trans_flux_div ./ res.total_diff;
res.rossby =   res.trans_flux_div ./ sqrt(res.cor_x.^2+res.cor_y.^2);
res.speed = sum_square_xy_grid(su,sv);

end
