clear all;
jobs = getenv('GC_JOBS');

if (false)
fig1 = figure(1);
clf;
set(fig1,'Position',[1 1 1200 1200]);

fig2 = figure(2);
clf;
set(fig2,'Position',[1 1 1200 1200]);
end 

fig3 = figure(3);
clf;
set(fig3,'Position',[1 1 1000 700]);

xtext = -3.5;
ytext = 21.5;


istart = -1;
iend = -1;
stride = 1;
minslices = 1;
control_job = 'central_final_nov17';
f_ice_control = strcat([jobs,'/',control_job,'/',control_job,'.out.nc']);
f_ocean_control = strcat([jobs,'/',control_job,'/plume.',control_job,'.out.nc']);
dice_control = nc_ice_read(f_ice_control,  istart,stride,iend,minslices);
dplume_control = nc_plume_read(f_ocean_control, istart,stride,iend);

lsurf_control = dice_control.lsurf(:,:,end);

precontrol_jobs = {'oct30_perturb_usq.11_sgd_flux_0.5', ...
		   'oct30_perturb_usq.11_sgd_flux_1.0', ...
		   'oct30_perturb_usq.11_sgd_flux_1.5'};

control_jobs = {'oct30_perturb_usq.8_sgd_flux_0.5',...
	       'no_tangle_oct30_perturb_usq_sgd_flux_1.0', ...
		'oct30_perturb_usq.8_sgd_flux_1.5'};

control_jobs_no_homotopy = {'oct30_perturb_usq.13_sgd_flux_0.5',...
			    'oct30_perturb_usq.13_sgd_flux_1.0',...
			    'oct30_perturb_usq.13_sgd_flux_1.5'};

control_jobs_sgd_type1 = {'oct30_perturb_usq.12_sgd_flux_0.5',...
			  'oct30_perturb_usq.12_sgd_flux_1.0',...
			  'oct30_perturb_usq.12_sgd_flux_1.5'};

for i=3:3

  pcj = precontrol_jobs{i}
  cj = control_jobs{i}
  cjnh = control_jobs_no_homotopy{i}
  cjst1 = control_jobs_sgd_type1{i}
  pre_f_ice = strcat([jobs,'/',pcj,'/',pcj,'.out.nc']);
  pre_f_ocean = strcat([jobs,'/',pcj,'/plume.',pcj,'.out.nc']);
  f_ice = strcat([jobs,'/',cj,'/',cj,'.out.nc']);
  f_ocean = strcat([jobs,'/',cj,'/plume.',cj,'.out.nc']);

  f_ice_nh = strcat([jobs,'/',cjnh,'/',cjnh,'.out.nc']);
  f_ocean_nh = strcat([jobs,'/',cjnh,'/plume.',cjnh,'.out.nc']);
  f_ice_st1 = strcat([jobs,'/',cjst1,'/',cjst1,'.out.nc']);
  f_ocean_st1 = strcat([jobs,'/',cjst1,'/plume.',cjst1,'.out.nc']);

pre_dice = nc_ice_read(pre_f_ice,  istart,stride,iend,minslices);
pre_dplume = nc_plume_read(pre_f_ocean, istart,stride,iend);
dice = nc_ice_read(f_ice,  istart,stride,iend,minslices);
dplume = nc_plume_read(f_ocean, istart,stride,iend);
dice_nh = nc_ice_read(f_ice_nh,  istart,stride,iend,minslices);
dplume_nh = nc_plume_read(f_ocean_nh, istart,stride,iend);
dice_st1 = nc_ice_read(f_ice_st1,  istart,stride,iend,minslices);
dplume_st1 = nc_plume_read(f_ocean_st1, istart,stride,iend);

x = dice.xgrid/1000.0;
y = dice.ygrid/1000.0;

pre_lsrf = pre_dice.lsurf(:,:,end);
pre_thk = pre_dice.thk(:,:,end);
lsrf = dice.lsurf(:,:,end);
thk = dice.thk(:,:,end);
lsrf_nh = dice_nh.lsurf(:,:,end);
thk_nh = dice_nh.thk(:,:,end);
lsrf_st1 = dice_st1.lsurf(:,:,end);
thk_st1 = dice_st1.thk(:,:,end);

pre_bmlt = pre_dice.bmlt(:,:,end);
bmlt = dice.bmlt(:,:,end);
bmlt_nh = dice_nh.bmlt(:,:,end);
bmlt_st1 = dice_st1.bmlt(:,:,end);

ncontours = 30;
fs = 14;
fslabel= 16;
jstart = 3;
jfinal = 83;
subset = @(M) M(:,jstart:jfinal);

x = x;
y = y(jstart:jfinal);

discharges = {'0.5','1.0','1.5'};
labels1 = {'a','b','c'};
labels2 = {'d','e','f'};
labels3 = {'g','h','i'};
labels4 = {'j','k','l'};


if (false)

subplot(4,3,i);
contourf(x,y,subset(pre_lsrf-0.0*lsurf_control)',ncontours,'EdgeColor','None');colorbar('FontSize',fs)
%caxis([-100 100]);
%xlabel('Across shelf (km)','FontSize',fs);
%ylabel('Along shelf (km)','FontSize',fs);
set(gca,'FontSize',fs);
title([' ',discharges{i},' km^3/a'],'FontSize',fs);
text(xtext,ytext,labels1{i},'FontSize',fslabel);

subplot(4,3,i+3);
contourf(x,y,subset(lsrf-0.0*lsurf_control)',ncontours,'EdgeColor','None');colorbar('FontSize',fs)
%caxis([-100 100])
%xlabel('Across shelf (km)','FontSize',fs);
%ylabel('Along shelf (km)','FontSize',fs);
set(gca,'FontSize',fs);
title([' ',discharges{i},' km^3/a'],'FontSize',fs);
text(xtext,ytext,labels2{i},'FontSize',fslabel);

subplot(4,3,i+6);
contourf(x,y,subset(lsrf_nh-0.0*lsurf_control)',ncontours,'EdgeColor','None');colorbar('FontSize',fs)
%caxis([-100 100]);
%xlabel('Across shelf (km)','FontSize',fs);
%ylabel('Along shelf (km)','FontSize',fs);
set(gca,'FontSize',fs);
title([' ',discharges{i},' km^3/a'],'FontSize',fs);
text(xtext,ytext,labels3{i},'FontSize',fslabel);

subplot(4,3,i+9);
contourf(x,y,subset(lsrf_st1-0.0*lsurf_control)',ncontours,'EdgeColor','None');colorbar('FontSize',fs)
%caxis([-100 100])
%xlabel('Across shelf (km)','FontSize',fs);
%ylabel('Along shelf (km)','FontSize',fs);
set(gca,'FontSize',fs);
title([' ',discharges{i},' km^3/a'],'FontSize',fs);
text(xtext,ytext,labels4{i},'FontSize',fslabel);
	 


subplot(4,3,i);
contourf(x,y,subset(pre_bmlt)',ncontours,'EdgeColor','None');colorbar('FontSize',fs)
caxis([0 110]);
%xlabel('Across shelf (km)','FontSize',fs);
%ylabel('Along shelf (km)','FontSize',fs);
set(gca,'FontSize',fs);
title([' ',discharges{i},' km^3/a'],'FontSize',fs);
text(xtext,ytext,labels1{i},'FontSize',fslabel);

subplot(4,3,i+3);
contourf(x,y,subset(bmlt)',ncontours,'EdgeColor','None');colorbar('FontSize',fs)
caxis([0 110]);
%xlabel('Across shelf (km)','FontSize',fs);
%ylabel('Along shelf (km)','FontSize',fs);
set(gca,'FontSize',fs);
title([' ',discharges{i},' km^3/a'],'FontSize',fs);
text(xtext,ytext,labels2{i},'FontSize',fslabel);

subplot(4,3,i+6);
contourf(x,y,subset(bmlt_nh)',ncontours,'EdgeColor','None');colorbar('FontSize',fs)
caxis([0 110]);
%xlabel('Across shelf (km)','FontSize',fs);
%ylabel('Along shelf (km)','FontSize',fs);
set(gca,'FontSize',fs);
title([' ',discharges{i},' km^3/a'],'FontSize',fs);
text(xtext,ytext,labels3{i},'FontSize',fslabel);

subplot(4,3,i+9);
contourf(x,y,subset(bmlt_st1)',ncontours,'EdgeColor','None');colorbar('FontSize',fs)
caxis([0 110]);
%xlabel('Across shelf (km)','FontSize',fs);
%ylabel('Along shelf (km)','FontSize',fs);
set(gca,'FontSize',fs);
title([' ',discharges{i},' km^3/a'],'FontSize',fs);
text(xtext,ytext,labels4{i},'FontSize',fslabel);

end

end



%%%% still have i==3

subplot(2,3,1);
contourf(x,y,subset(pre_bmlt)',ncontours,'EdgeColor','None');colorbar('FontSize',fs)
caxis([0 110]);
%xlabel('Across shelf (km)','FontSize',fs);
%ylabel('Along shelf (km)','FontSize',fs);
set(gca,'FontSize',fs);
%title('melt rates, pre-control','FontSize',fs);
title('Melt rates, pre-control','FontSize',fs);
text(xtext,ytext,'a','FontSize',fslabel);

subplot(2,3,3);
contourf(x,y,subset(bmlt)',ncontours,'EdgeColor','None');colorbar('FontSize',fs)
caxis([0 110]);
%xlabel('Across shelf (km)','FontSize',fs);
%ylabel('Along shelf (km)','FontSize',fs);
set(gca,'FontSize',fs);
title('Melt rates, control','FontSize',fs);
text(xtext,ytext,'c','FontSize',fslabel);

%subplot(2,3,3);
%contourf(x,y,subset(bmlt_nh)',ncontours,'EdgeColor','None');colorbar('FontSize',fs)
%caxis([0 110]);
%%xlabel('Across shelf (km)','FontSize',fs);
%%ylabel('Along shelf (km)','FontSize',fs);
%set(gca,'FontSize',fs);
%title('from control','FontSize',fs);
%text(xtext,ytext,'c.','FontSize',fslabel);

subplot(2,3,2);
contourf(x,y,subset(bmlt_st1)',ncontours,'EdgeColor','None');colorbar('FontSize',fs)
caxis([0 110]);
%xlabel('Across shelf (km)','FontSize',fs);
%ylabel('Along shelf (km)','FontSize',fs);
set(gca,'FontSize',fs);
title('Melt rates, pre-control, jets','FontSize',fs);
text(xtext,ytext,'b','FontSize',fslabel);

subplot(2,3,4);
contourf(x,y,subset(pre_thk)',ncontours,'EdgeColor','None');colorbar('FontSize',fs)
caxis([0 600]);
%xlabel('Across shelf (km)','FontSize',fs);
%ylabel('Along shelf (km)','FontSize',fs);
set(gca,'FontSize',fs);
title('Ice thickness, pre-control','FontSize',fs);
text(xtext,ytext,'d','FontSize',fslabel);

subplot(2,3,6);
contourf(x,y,subset(thk)',ncontours,'EdgeColor','None');colorbar('FontSize',fs)
caxis([0 600]);
%xlabel('Across shelf (km)','FontSize',fs);
%ylabel('Along shelf (km)','FontSize',fs);
set(gca,'FontSize',fs);
title('Ice thickness, control','FontSize',fs);
text(xtext,ytext,'f','FontSize',fslabel);

%subplot(2,3,6);
%contourf(x,y,subset(thk_nh)',ncontours,'EdgeColor','None');colorbar('FontSize',fs)
%caxis([0 600]);
%%xlabel('Across shelf (km)','FontSize',fs);
%%ylabel('Along shelf (km)','FontSize',fs);
%set(gca,'FontSize',fs);
%%title([' ',discharges{i},' km^3/a'],'FontSize',fs);
%text(xtext,ytext,'f.','FontSize',fslabel);

subplot(2,3,5);
contourf(x,y,subset(thk_st1)',ncontours,'EdgeColor','None');colorbar('FontSize',fs)
caxis([0 600]);
%xlabel('Across shelf (km)','FontSize',fs);
%ylabel('Along shelf (km)','FontSize',fs);
set(gca,'FontSize',fs);
title('Ice thickness, pre-control, jets','FontSize',fs);
text(xtext,ytext,'e','FontSize',fslabel);

