trim_edge = @(M) M(4:end-3,4:end-3,:);
fs = 14;
fs_label = 16;
label_x = -2.5;
label_y = 21;

selectset = @(M) M(2:79,2:80);

%subset = @(M)selectset(four_by_four_avg(M));
%subset = @(M)selectset(two_by_two_avg(M));
subset = @(M)selectset(three_by_three_avg(M));

x = res.x; 
y = res.y; 
[y,x] = meshgrid(y,x);

figure(1);
clf;
ncontours = 20;
cmax = 1.5e-3;

subplot(2,3,1);
contourf(subset(x),subset(y), ...
         subset(res.u_flux_div),ncontours,'EdgeColor','None');colorbar;
caxis([-cmax cmax ]);
title('flux divergence','FontSize',fs);
set(gca,'FontSize',fs);
text(label_x,label_y,'a','color','k','FontSize',fs_label);

subplot(2,3,3);
contourf(subset(x),subset(y), ...
	 subset(res.u_diff),ncontours,'EdgeColor','None');colorbar;
caxis([-cmax cmax]);
title('diffusion term','FontSize',fs)
set(gca,'FontSize',fs);
text(label_x,label_y,'c','color','k','FontSize',fs_label);

subplot(2,3,2);
contourf(subset(x),subset(y), ...
	 -subset(res.cor_x),ncontours,'EdgeColor','None');colorbar;
title('coriolis acceleration','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);
text(label_x,label_y,'b','color','k','FontSize',fs_label);

subplot(2,3,4);
contourf(subset(x),subset(y), ...
	 subset(res.isopyc_grad_x),ncontours,...
	 'EdgeColor','None');colorbar;
title('interface flattening term','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);
text(label_x,label_y,'d','color','k','FontSize',fs_label);

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
text(label_x,label_y,'e','color','k','FontSize',fs_label);

subplot(2,3,6);
%contourf(subset(x),subset(y), ...
%	 subset(res.drag_x),'EdgeColor','None');colorbar;
%title('Drag term','FontSize',fs);
contourf(subset(x),subset(y), ...
	 subset(res.u_trans_detrain),ncontours,'EdgeColor','None');colorbar;
title('detrainment term','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);
text(label_x,label_y,'f','color','k','FontSize',fs_label);

%%%%%%%% same thing for y direction
figure(2);
clf;
ncontours = 20;
cmax = 1.25e-3;

subplot(2,3,1);
contourf(subset(x),subset(y), ...
         subset(res.v_flux_div),ncontours,'EdgeColor','None');colorbar;
caxis([-cmax cmax ]);
title('transport flux divergence','FontSize',fs);
title('flux divergence','FontSize',fs);
set(gca,'FontSize',fs);
text(label_x,label_y,'a','color','k','FontSize',fs_label);

subplot(2,3,3);
contourf(subset(x),subset(y), ...
	 subset(res.v_diff),ncontours,'EdgeColor','None');colorbar;
caxis([-cmax cmax]);
title('total transport diffusion','FontSize',fs)
title('diffusion term','FontSize',fs)
set(gca,'FontSize',fs);
text(label_x,label_y,'c','color','k','FontSize',fs_label);

subplot(2,3,2);
contourf(subset(x),subset(y), ...
	 -subset(res.cor_y),ncontours,'EdgeColor','None');colorbar;
title('coriolis accel','FontSize',fs);
title('coriolis acceleration','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);
text(label_x,label_y,'b','color','k','FontSize',fs_label);

subplot(2,3,4);
contourf(subset(x),subset(y), ...
	 subset(res.isopyc_grad_y),ncontours,...
	 'EdgeColor','None');colorbar;
title('isopycnal flattening','FontSize',fs);
title('interface flattening term','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);
text(label_x,label_y,'d','color','k','FontSize',fs_label);

subplot(2,3,5);
contourf(subset(x),subset(y), ...
	 subset(res.den_grad_y),ncontours, ...
	 'EdgeColor','None');colorbar;
title('Density gradient accel.','FontSize',fs);
title('density gradient term','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);
text(label_x,label_y,'e','color','k','FontSize',fs_label);

subplot(2,3,6);
contourf(subset(x),subset(y), ...
	      subset(res.drag_y),ncontours,'EdgeColor','None');colorbar;
contourf(subset(x),subset(y), ...
	 subset(res.v_trans_detrain+0.0*res.isopyc_grad_y+ ...
             0.0*res.v_diff-0.0*res.cor_y+0.0*res.drag_y), ...
	 ncontours,'EdgeColor','None');colorbar;
%title('Detrainment term','FontSize',fs);
%title('Drag term','FontSize',fs);
title('detrainment term','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);
text(label_x,label_y,'f','color','k','FontSize',fs_label);

%%%%%%%%%%%%%%%%
%  More scatter plots showing balances in x and y direction
%%%%%%%%%%%%%%%

figure(4);
clf;

subplot(2,2,1);
hold on
a = flatten_field(subset(res.isopyc_grad_x+0.0*res.den_grad_x));
b = flatten_field(subset(res.cor_x));
limlo = -2.0e-3;
limhi = 0.5e-3;
plot([limlo, limhi],[limlo, limhi],'r-');
xlim([limlo limhi]);
plot(a,b,'b.');
m = corr([a,b]);
%title('x balance','FontSize',fs);
xlabel('g''DA_x (m^2/s^2)','FontSize',fs);
ylabel('-fDV (m^2/s^2)','FontSize',fs);
set(gca,'FontSize',fs);
title(strcat(['correlation =',sprintf('%3.2f',m(1,2))]),'FontSize',fs);
text(-19e-4,3.0e-4,'a','color','k','FontSize',fs_label);

subplot(2,2,2);
hold on
a = flatten_field(subset(1.0*res.isopyc_grad_y+1.0*res.v_trans_detrain));
b= flatten_field(subset(1.0*res.v_flux_div));
m = corr([a b]);
plot(a,b,'b.');

limlo = -1.0e-3;
limhi = 1.0e-3;
plot([limlo, limhi],[limlo, limhi],'r-');
xlim([limlo limhi]);
ylim([limlo limhi]);

set(gca,'FontSize',fs);
xlabel('g''DA_y - d''V (m^2/s^2)','FontSize',fs);
ylabel('(DUU)_x+(DUV)_y (m^2/s^2)','FontSize',fs);
title(strcat(['correlation =',sprintf('%3.2f',m(1,2))]),'FontSize',fs);
text(-0.9e-3,0.8e-3,'b','color','k','FontSize',fs_label);

subplot(2,2,4);
hold on
a = flatten_field(subset(1.0*res.isopyc_grad_y./res.pdep + res.entr_drag));
b= flatten_field(subset(1.0*res.sv_adv_derv));
m = corr([a b]);
plot(a,b,'b.');
limlo = -5.0e-5;
limhi = 5.0e-5;
plot([limlo, limhi],[limlo, limhi],'r-');
xlim([limlo limhi]);
ylim([limlo limhi]);

set(gca,'FontSize',fs);
xlabel('g'' A_y - d''V/D (m/s^2)','FontSize',fs);
ylabel('UV_x + VV_y  (m/s^2)','FontSize',fs);
title(strcat(['correlation =',sprintf('%3.2f',m(1,2))]),'FontSize',fs);
text(-4.5e-5,4.25e-5,'d','color','k','FontSize',fs_label);

subplot(2,2,3);
hold on
a = flatten_field(subset(res.isopyc_grad_x./res.pdep));
b= flatten_field(subset(res.cor_y./res.pdep));
m = corr([a b]);
plot(a,b,'b.');
limlo = -5.0e-5;
limhi = 5.0e-5;
plot([limlo, limhi],[limlo, limhi],'r-');
xlim([limlo limhi]);
ylim([limlo limhi]);

set(gca,'FontSize',fs);
xlabel('g'' A_x (m/s^2)','FontSize',fs);
ylabel('-fV  (m/s^2)','FontSize',fs);
title(strcat(['correlation =',sprintf('%3.2f',m(1,2))]),'FontSize',fs);
text(-4.5e-5,4.25e-5,'c','color','k','FontSize',fs_label);

if (false)
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
end

if (false)

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

end
