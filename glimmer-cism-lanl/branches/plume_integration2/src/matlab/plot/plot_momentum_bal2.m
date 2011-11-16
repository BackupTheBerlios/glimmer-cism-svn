trim_edge = @(M) M(4:end-3,4:end-3,:);
fs = 14;

selectset = @(M) M(1:76,1:155,:);
selectset = @(M) M(1:76,1:80);

subset = @(M)selectset(four_by_four_avg(M));

x = res.x; %dplume.xgrid/1000.0 - 1.375;
y = res.y; %rdplume.ygrid/1000.0 - 875;
[y,x] = meshgrid(y,x);


figure(1);
clf;
ncontours = 20;

cmax = 1.0e-4;

subplot(2,3,1);
contourf(subset(x),subset(y), ...
         subset(res.su_adv_derv),ncontours,'EdgeColor','None');colorbar;
caxis([-cmax cmax ]);
title('x acceleration','FontSize',fs);
set(gca,'FontSize',fs);

subplot(2,3,2);
contourf(subset(x),subset(y), ...
	 subset(res.u_diff+res.u_diff2),ncontours,'EdgeColor','None');colorbar;
caxis([-cmax cmax]);
title('diffusion terms','FontSize',fs)
set(gca,'FontSize',fs);

subplot(2,3,3);
contourf(subset(x),subset(y), ...
	 -subset(res.cor_x),ncontours,'EdgeColor','None');colorbar;
title('coriolis acceleration','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);

subplot(2,3,4);
contourf(subset(x),subset(y), ...
	 subset(res.isopyc_grad_x),ncontours,...
	 'EdgeColor','None');colorbar;
title('interface flattening term','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);

subplot(2,3,5);
contourf(subset(x),subset(y), ...
	 subset(res.den_grad_x), ncontours,...
	 'EdgeColor','None');colorbar;
%contourf(subset(x),subset(y), ...
%	 subset(sqrt(res.u_adv_derv.^2+res.v_adv_derv.^2)), ...
%	 'EdgeColor','None');colorbar;
title('density gradient term','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);


subplot(2,3,6);
%contourf(subset(x),subset(y), ...
%	 subset(res.drag_x),'EdgeColor','None');colorbar;
%title('Drag term','FontSize',fs);
contourf(subset(x),subset(y), ...
	 subset(res.u_boundary),ncontours,'EdgeColor','None');colorbar;
title('boundary term','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);

%%%%%%%% same thing for y direction
figure(2);
clf;
ncontours = 20;
cmax = 1.0e-4;

subplot(2,3,1);
contourf(subset(x),subset(y), ...
         subset(res.sv_adv_derv),ncontours,'EdgeColor','None');colorbar;
caxis([-cmax cmax ]);
title('transport flux divergence','FontSize',fs);
title('y acceleration','FontSize',fs);
set(gca,'FontSize',fs);

subplot(2,3,2);
contourf(subset(x),subset(y), ...
	 subset(res.v_diff+res.v_diff2),ncontours,'EdgeColor','None');colorbar;
caxis([-cmax cmax]);
title('diffusion terms','FontSize',fs)
set(gca,'FontSize',fs);

subplot(2,3,3);
contourf(subset(x),subset(y), ...
	 -subset(res.cor_y),ncontours,'EdgeColor','None');colorbar;
title('coriolis accel','FontSize',fs);
title('coriolis acceleration','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);

subplot(2,3,4);
contourf(subset(x),subset(y), ...
	 subset(res.isopyc_grad_y),ncontours,...
	 'EdgeColor','None');colorbar;
title('isopycnal flattening','FontSize',fs);
title('interface flattening term','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);

subplot(2,3,5);
contourf(subset(x),subset(y), ...
	 subset(res.den_grad_y),ncontours, ...
	 'EdgeColor','None');colorbar;
title('Density gradient accel.','FontSize',fs);
title('density gradient term','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);

subplot(2,3,6);
contourf(subset(x),subset(y), ...
	      subset(res.v_boundary),ncontours,'EdgeColor','None');colorbar;
title('boundary term','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);


%%%%%%%%%%%%%%%%
%  More scatter plots showing balances in x and y direction
%%%%%%%%%%%%%%%

figure(4);
clf;

subplot(1,2,1);
hold on
a = flatten_field(subset(1.0*res.isopyc_grad_x + 1.0*res.den_grad_x + ...
                         1.0*res.u_diff+1.0*res.u_diff2+...
                         1.0*res.u_boundary));
b= flatten_field(subset(1.0*res.su_adv_derv+1.0*res.cor_x));

limlo = -0.5e-4;
limhi = 0.5e-4;
plot([limlo, limhi],[limlo, limhi],'r-');
%xlim([limlo limhi]);
plot(a,b,'b.');
m = corr([a,b]);
%title('x balance','FontSize',fs);
xlabel('interface flattening term (m^2/s^2)','FontSize',fs);
ylabel('Coriolis acceleration (m^2/s^2)','FontSize',fs);
set(gca,'FontSize',fs);
title(strcat(['correlation =',sprintf('%3.2f',m(1,2))]),'FontSize',fs);

subplot(1,2,2);
hold on
a = flatten_field(subset(1.0*res.isopyc_grad_y + 1.0*res.den_grad_y + ...
                         1.0*res.v_diff+1.0*res.v_diff2+...
                         1.0*res.v_boundary));
b= flatten_field(subset(1.0*res.sv_adv_derv+1.0*res.cor_y));

m = corr([a b]);
plot(a,b,'b.');
limlo = -0.5e-4;
limhi = 0.5e-4;
plot([limlo, limhi],[limlo, limhi],'r-');
%xlim([limlo limhi]);
%ylim([limlo limhi]);

set(gca,'FontSize',fs);
xlabel('interface flattening term + detrainment term (m^2/s^2)','FontSize',fs);
ylabel('non-Coriolis inertial terms (m^2/s^2)','FontSize',fs);
title(strcat(['correlation =',sprintf('%3.2f',m(1,2))]),'FontSize',fs);

figure(5);
clf;
contourf(subset(x),subset(y),subset(res.su_adv_derv+res.cor_x - res.u_diff-res.u_diff2-res.den_grad_x-res.isopyc_grad_x-res.u_boundary),'EdgeColor','None');colorbar;

figure(6);
clf;
contourf(subset(x),subset(y),subset(res.sv_adv_derv+res.cor_y - res.v_diff-res.v_diff2-res.den_grad_y-res.isopyc_grad_y-res.v_boundary),'EdgeColor','None');colorbar;
