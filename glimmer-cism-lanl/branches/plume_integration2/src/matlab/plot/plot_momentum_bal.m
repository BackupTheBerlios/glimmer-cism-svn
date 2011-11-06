trim_edge = @(M) M(4:end-3,4:end-3,:);
fs = 14;

subset = @(M) M(1:40,1:50,:);

x = dplume.xgrid/1000.0;
y = dplume.ygrid/1000.0;

%[x,y] = meshgrid(x,y);

figure(1);
clf;

cmax = 5e-4;

subplot(2,3,1);
contourf(subset(x),subset(y), ...
         subset(res.u_flux_div),'EdgeColor','None');colorbar;
caxis([0 cmax ]);
title('transport flux divergence','FontSize',fs);

subplot(2,3,2);
contourf(subset(x),subset(y), ...
	 subset(res.u_diff),'EdgeColor','None');colorbar;
caxis([0 cmax]);
title('total transport diffusion','FontSize',fs)

subplot(2,3,3);
contourf(subset(x),subset(y), ...
	 subset(res.cor_x),'EdgeColor','None');colorbar;
title('coriolis accel','FontSize',fs);
caxis([0 cmax]);

subplot(2,3,4);
contourf(subset(x),subset(y), ...
	 subset(res.isopyc_grad_x),...
	 'EdgeColor','None');colorbar;
title('isopycnal flattening','FontSize',fs);
caxis([0 cmax]);

subplot(2,3,5);
contourf(subset(x),subset(y), ...
	 subset(res.den_grad_x), ...
	 'EdgeColor','None');colorbar;
%contourf(subset(x),subset(y), ...
%	 subset(sqrt(res.u_adv_derv.^2+res.v_adv_derv.^2)), ...
%	 'EdgeColor','None');colorbar;
title('Density gradient accel.','FontSize',fs);
caxis([0 cmax]);


subplot(2,3,6);
contourf(subset(x),subset(y), ...
	 subset(res.drag_x),'EdgeColor','None');colorbar;
title('Drag term','FontSize',fs);
caxis([0 cmax]);

%figure(2);
%plot(flatten_field(res.rossby),flatten_field(dplume.bmelt(:,:,1)),'.');

figure(3);
subplot(2,2,1);
contourf(x,y,res.u_diff,'EdgeColor','None');colorbar;
title('u diffusion term','FontSize',fs);

subplot(2,2,2);
contourf(x,y,res.v_diff,'EdgeColor','None');colorbar;
title('v diffusion term','FontSize',fs);

subplot(2,2,3);
%contourf(x,y,log(res.reynolds),'EdgeColor','None');colorbar;
%title('log(Reynolds number)','FontSize',fs);
contourf(x,y,res.v_flux_div,'EdgeColor','None');colorbar;
title('u flux divergence term','FontSize',fs);

subplot(2,2,4);
%contourf(x,y,res.drag,'EdgeColor','None');colorbar;
%title('non-dimensional drag','FontSize',fs);
contourf(x,y,res.v_flux_div,'EdgeColor','None');colorbar;
title('u flux divergence term','FontSize',fs);


figure(4);
clf;
%contourf(x,y,res.delta_w,'EdgeColor','None');
%colorbar;
subplot(2,2,1);
hold on
plot(trim_edge(res.delta_w),trim_edge(mean(dplume.train,3))/(3600.0*24.0*365.25),'k.');
plot(res.delta_w,res.delta_w,'b.');
xlabel('delta w (m/s)','FontSize',fs);
ylabel('trainment rate (m/s)','FontSize',fs);
hold off
title('Velocity divergence in plume','FontSize',fs);

subplot(2,2,2);
hold on
plot(trim_edge(res.adv_derv_pdep),trim_edge(mean(dplume.train,3))/(3600.0*24.0*365.25),'k.');
plot(res.adv_derv_pdep,res.adv_derv_pdep,'b.');
xlabel('adv derivative of pdep (m/s)','FontSize',fs);
ylabel('trainment rate (m/s)','FontSize',fs);
hold off
title('Advective derivative of pdep','FontSize',fs);

subplot(2,2,3);
hold on
plot(trim_edge(res.flux_div),trim_edge(mean(dplume.train,3))/(3600.0*24.0*365.25),'k.');
plot(res.flux_div,res.flux_div,'b.');
xlabel('flux div (m/s)','FontSize',fs);
ylabel('trainment rate (m/s)','FontSize',fs);
hold off
title('Total flux divergnece','FontSize',fs);

subplot(2,2,4);
hold on
bmlt =  dplume.bmelt/(3600.0*24*365.25);
%plot(bmlt,dplume.train/(3600.0*24.0*365.25),'k.');
%plot(bmlt,bmlt,'b.');
plot(trim_edge(res.flux_div),trim_edge(res.delta_w+res.adv_derv_pdep),'r.');
plot(res.flux_div,res.flux_div,'b.');
xlabel('flux_div (m/s)','FontSize',fs);
ylabel('sum of two terms (m/s)','FontSize',fs);
hold off
title('bmlt','FontSize',fs);


