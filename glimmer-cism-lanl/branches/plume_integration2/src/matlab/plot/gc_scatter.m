function [] = gc_scatter(flat_ocean,flat_ice,flat_ocean_transient,flat_ice_transient,amb_t_func,fig_dir)

doTransient = false;
scatter_xsize = 1000;
scatter_ysize = 600;

fs = 16;
fn = 'Helvetica';

c1 = 'm';
c1 = 'k';
c2 = 'k';
c_transient = 'b.';
c_steady = 'r.';
c_steady = 'b.';
c5 = 'r-';

thermal_forcing = @(depth,salt,temp) temp-(-5.73e-2*salt-7.61d-4*depth+8.32d-2);

if (doTransient)

flattforcing = thermal_forcing(abs(flat_ocean_transient.draft), ...
                               flat_ocean_transient.salt, ...
                               flat_ocean_transient.temp);

figure(100);
clf;
set(figure(100),'Position',[1 1 scatter_xsize scatter_ysize]);

subplot(1,3,1);
hold on
plot(flat_ocean_transient.speed,flat_ocean_transient.bmlt,c_transient);

subplot(1,3,2);
hold on
plot(flattforcing,flat_ocean_transient.bmlt,c_transient);

subplot(1,3,3);
hold on
plot(flat_ocean_transient.train,flat_ocean_transient.bmlt,c_transient);

figure(101);
%clf;

subplot(1,3,1);
hold on
plot(flat_ocean_transient.grad,flat_ocean_transient.bmlt,c_transient);

subplot(1,3,2);
hold on
plot(flat_ice_transient.y_adv,flat_ocean_transient.bmlt,c_transient);

subplot(1,3,3);
hold on
plot(-flat_ocean_transient.draft,flat_ocean_transient.bmlt,c_transient);

end

flattforcing = thermal_forcing(abs(flat_ocean.draft),flat_ocean.salt,flat_ocean.temp);
flatambtemp = amb_t_func(-flat_ocean.draft);
flat_exp_var = flat_ocean.gradx.^(1/4);
flat_exp_var = flat_ocean.speed.*flattforcing;

corr_mat = corrcoef([flat_ocean.bmlt,flat_ocean.speed,flattforcing,flat_ocean.train, ...
                     flat_ocean.grad,flat_ocean.gradx,flat_ocean.grady,flat_ice.y_adv,-flat_ocean.draft,...
			     flattforcing.*flat_ocean.grad,flatambtemp,flatambtemp.*flat_ocean.grad,flat_ocean.gradx,...
                       flat_exp_var]);

speed_cor = corr_mat(1,2);
temp_cor =  corr_mat(1,3);
train_cor =  corr_mat(1,4);
grad_cor =  corr_mat(1,5);
draft_cor=  corr_mat(1,9);
yadv_cor =  corr_mat(1,8);
tforc_grad_cor = corr_mat(1,10);
ambt_cor = corr_mat(1,11);
ambt_grad_cor = corr_mat(1,12);
gradx_cor = corr_mat(1,13);
exp_var_cor = corr_mat(1,14);

figure(100);
if (not(doTransient))
  clf;
  set(figure(100),'Position',[1 1 scatter_xsize scatter_ysize]);
end

subplot(1,3,1);
plot(flat_ocean.speed,flat_ocean.bmlt,c_steady);
title(strcat(['correlation = ', sprintf('%4.3f',speed_cor)]),'FontSize',fs,'Color',c1,'FontName',fn);
set(gca,'FontSize',fs);
xlabel('plume speed (m/s)','FontSize',fs,'FontName',fn);
ylabel('melt rate (m/year)','FontSize',fs,'FontName',fn);
hold off

subplot(1,3,2);
plot(flattforcing,flat_ocean.bmlt,c_steady);
title(strcat(['correlation = ', sprintf('%4.3f',temp_cor)]),'FontSize',fs,'Color',c2,'FontName',fn);
set(gca,'FontSize',fs);
xlabel('thermal forcing (C)','FontSize',fs,'FontName',fn);
ylabel('melt rate (m/year)','FontSize',fs,'FontName',fn);
hold off

subplot(1,3,3);
plot(flat_ocean.train,flat_ocean.bmlt,c_steady);
%plot(abs(flat_ocean.sv),flat_ocean_transient.bmlt,c_steady);

title(strcat(['correlation = ', sprintf('%4.3f',train_cor)]),'FontSize',fs,'Color',c2,'FontName',fn);
set(gca,'FontSize',fs);
xlabel('entrainment (km/year)','FontSize',fs,'FontName',fn);
xlim([-500 100]);
%xlim([0 0.4]);
ylabel('melt rate (m/year)','FontSize',fs,'FontName',fn);
hold off

print('-depsc',strcat([fig_dir,'/plume_scatter']));

figure(101);
if (not(doTransient))
  clf;
end
set(figure(101),'Position',[1 1 scatter_xsize scatter_ysize]);

subplot(1,3,1);
plot(flat_ocean.grad,flat_ocean.bmlt,c_steady);
set(gca,'FontSize',fs);
title(strcat(['correlation = ', sprintf('%4.3f',grad_cor)]),'FontSize',fs,'Color',c2,'FontName',fn);
xlabel('ice gradient (non-dimensional)','FontSize',fs,'FontName',fn);
ylabel('melt rate (m/year)','FontSize',fs,'FontName',fn);
hold off

subplot(1,3,2);
set(gca,'FontSize',fs);
if(not(doTransient))
  hold on;
end
plot(flat_ice.y_adv,flat_ocean.bmlt,c_steady);
plot(0:0.5:75,0:0.5:75,c5,'linewidth',2.0);
title(strcat(['correlation = ', sprintf('%4.3f',yadv_cor)]),'FontSize',fs,'Color',c1,'FontName',fn);
xlim([-5 80]);
xlabel('v_0 H_y (m/year)','FontSize',fs,'FontName',fn);
ylabel('melt rate (m/year)','FontSize',fs,'FontName',fn);
hold off


subplot(1,3,3);
plot(-flat_ocean.draft,flat_ocean.bmlt,c_steady);
title(strcat(['correlation = ', sprintf('%4.3f',draft_cor)]),'FontSize',fs,'Color',c2,'FontName',fn);
set(gca,'FontSize',fs);
xlabel('ice draft (m)','FontSize',fs,'FontName',fn);
ylabel('melt rate (m/year)','FontSize',fs,'FontName',fn);
xlim([0 600]);hold off

print('-depsc',strcat([fig_dir,'/ice_scatter']));

if (false)
figure(102);
if (not(doTransient))
  clf;
end

subplot(1,3,1);
plot(flat_ocean.grad.*flattforcing,flat_ocean.bmlt,c_steady);
set(gca,'FontSize',fs);
title(strcat(['correlation = ', sprintf('%4.3f',tforc_grad_cor)]),'FontSize',fs,'Color',c2,'FontName',fn);
xlabel('ice gradient * thermal forcing','FontSize',fs,'FontName',fn);
ylabel('melt rate (m/year)','FontSize',fs,'FontName',fn);
hold off

subplot(1,3,2);
plot(flatambtemp,flat_ocean.bmlt,c_steady);
title(strcat(['correlation = ', sprintf('%4.3f',ambt_cor)]),'FontSize',fs,'Color',c1,'FontName',fn);
set(gca,'FontSize',fs);
xlabel('ambient water temp','FontSize',fs,'FontName',fn);
ylabel('melt rate (m/year)','FontSize',fs,'FontName',fn);
hold off

subplot(1,3,3);
plot(flat_exp_var,flat_ocean.bmlt,c_steady);
title(strcat(['correlation = ', sprintf('%4.3f',exp_var_cor)]),'FontSize',fs,'Color',c2,'FontName',fn);
set(gca,'FontSize',fs);
xlabel('exp var','FontSize',fs,'FontName',fn);
ylabel('melt rate (m/year)','FontSize',fs,'FontName',fn);
hold off

end

coeffs = corrcoef([flat_ice.ice_def, flat_ocean.gradx, flat_ocean.grady, flat_ocean.grad]);
gradx_cor = coeffs(1,2);
grady_cor = coeffs(1,3);
grad_cor = coeffs(1,4);

if (false)
figure(103);
clf;
subplot(1,2,1);
plot(flat_ocean.gradx,flat_ice.ice_def,c_steady);
set(gca,'FontSize',fs);
title(strcat(['correlation = ', sprintf('%4.3f',gradx_cor)]),'FontSize',fs,'Color',c2,'FontName',fn);
ylabel('ice deformation (m/year)','FontSize',fs,'FontName',fn);
xlabel('ice basal surface gradient, x-direction','FontSize',fs,'FontName',fn);
hold off

subplot(1,2,2);
plot(flat_ocean.grady,flat_ice.ice_def,c_steady);
title(strcat(['correlation = ', sprintf('%4.3f',grady_cor)]),'FontSize',fs,'Color',c1,'FontName',fn);
set(gca,'FontSize',fs);
ylabel('ice deformation (m/year)','FontSize',fs,'FontName',fn);
xlabel('ice basal surface gradient, y-direction','FontSize',fs,'FontName',fn);
hold off

%subplot(1,3,3);
%plot(flat_ocean.gradx,-flat_ice.thk_div,c_steady);
%title(strcat(['correlation = ', sprintf('%4.3f',grad_cor)]),'FontSize',fs,'Color',c2,'FontName',fn);
%set(gca,'FontSize',fs);
%ylabel('ice divergence (m/year)','FontSize',fs,'FontName',fn);
%xlabel('ice thickness gradient, x-direction','FontSize',fs,'FontName',fn);
%hold off

end

end
