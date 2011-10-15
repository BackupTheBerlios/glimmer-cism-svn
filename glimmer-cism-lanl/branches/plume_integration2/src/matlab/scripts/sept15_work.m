jobDir = '/scratch/cvg222/gcp/GC_jobs_paper/';
plume_files = {};
plume_files{1} = [jobDir,'v2_sept_perturb_acab_0.0/plume.v2_sept_perturb_acab_0.0.out.nc'];
plume_files{2} = [jobDir,'v2_sept_perturb_acab_-2.0/plume.v2_sept_perturb_acab_-2.0.out.nc'];
plume_files{3} = [jobDir,'v2_sept_perturb_bottom_-0.15/plume.v2_sept_perturb_bottom_-0.15.out.nc'];
plume_files{4} = [jobDir,'v2_sept_perturb_bottom_-0.35/plume.v2_sept_perturb_bottom_-0.35.out.nc'];
plume_files{5} = [jobDir,'v2_sept_perturb_cdb_0.00125/plume.v2_sept_perturb_cdb_0.00125.out.nc'];
plume_files{6} = [jobDir,'v2_sept_perturb_cdb_0.00375/plume.v2_sept_perturb_cdb_0.00375.out.nc'];
plume_files{7} = [jobDir,'v2_sept_perturb_entype2_5/plume.v2_sept_perturb_entype2_5.out.nc'];
plume_files{8} = [jobDir,'central_final/plume.central_final.out.nc'];
plume_files{9} = [jobDir,'v2_sept_perturb_phi_60.0/plume.v2_sept_perturb_phi_60.0.out.nc'];
plume_files{10} = [jobDir,'v2_sept_perturb_plume_min_thickness_10.0/plume.v2_sept_perturb_plume_min_thickness_10.0.out.nc'];
plume_files{11} = [jobDir,'v2_sept_perturb_plume_min_thickness_15.0/plume.v2_sept_perturb_plume_min_thickness_15.0.out.nc'];
plume_files{12} = [jobDir,'v2_sept_perturb_plume_min_thickness_5.0/plume.v2_sept_perturb_plume_min_thickness_5.0.out.nc'];
plume_files{13} = [jobDir,'v2_sept_perturb_tauxy0_-12500.0/plume.v2_sept_perturb_tauxy0_-12500.0.out.nc'];
plume_files{14} = [jobDir,'v2_sept_perturb_tauxy0_-50000/plume.v2_sept_perturb_tauxy0_-50000.out.nc'];
plume_files{15} = [jobDir,'v2_sept_perturb_temp_interior_-15.0/plume.v2_sept_perturb_temp_interior_-15.0.out.nc'];
plume_files{16} = [jobDir,'v2_sept_perturb_temp_interior_-5.0/plume.v2_sept_perturb_temp_interior_-5.0.out.nc'];
plume_files{17} = [jobDir,'v2_sept_perturb_top_-0.85/plume.v2_sept_perturb_top_-0.85.out.nc'];
plume_files{18} = [jobDir,'v2_sept_perturb_top_-1.35/plume.v2_sept_perturb_top_-1.35.out.nc'];
plume_files{19} = [jobDir,'v2_sept_perturb_top_-1.85/plume.v2_sept_perturb_top_-1.85.out.nc'];
plume_files{20} = [jobDir,'v2_sept_perturb_train_time_const_1800.0/plume.v2_sept_perturb_train_time_const_1800.0.out.nc'];
plume_files{21} = [jobDir,'v2_sept_perturb_train_time_const_5400.0/plume.v2_sept_perturb_train_time_const_5400.0.out.nc'];
plume_files{22} = [jobDir,'v2_sept_perturb_visc_15.0/plume.v2_sept_perturb_visc_15.0.out.nc'];


ice_files = {};
ice_files{1} = [jobDir,'v2_sept_perturb_acab_0.0/v2_sept_perturb_acab_0.0.out.nc'];
ice_files{2} = [jobDir,'v2_sept_perturb_acab_-2.0/v2_sept_perturb_acab_-2.0.out.nc'];
ice_files{3} = [jobDir,'v2_sept_perturb_bottom_-0.15/v2_sept_perturb_bottom_-0.15.out.nc'];
ice_files{4} = [jobDir,'v2_sept_perturb_bottom_-0.35/v2_sept_perturb_bottom_-0.35.out.nc'];
ice_files{5} = [jobDir,'v2_sept_perturb_cdb_0.00125/v2_sept_perturb_cdb_0.00125.out.nc'];
ice_files{6} = [jobDir,'v2_sept_perturb_cdb_0.00375/v2_sept_perturb_cdb_0.00375.out.nc'];
ice_files{7} = [jobDir,'v2_sept_perturb_entype2_5/v2_sept_perturb_entype2_5.out.nc'];
ice_files{8} = [jobDir,'central_final/central_final.out.nc'];
ice_files{9} = [jobDir,'v2_sept_perturb_phi_60.0/v2_sept_perturb_phi_60.0.out.nc'];
ice_files{10} = [jobDir,'v2_sept_perturb_plume_min_thickness_10.0/v2_sept_perturb_plume_min_thickness_10.0.out.nc'];
ice_files{11} = [jobDir,'v2_sept_perturb_plume_min_thickness_15.0/v2_sept_perturb_plume_min_thickness_15.0.out.nc'];
ice_files{12} = [jobDir,'v2_sept_perturb_plume_min_thickness_5.0/v2_sept_perturb_plume_min_thickness_5.0.out.nc'];
ice_files{13} = [jobDir,'v2_sept_perturb_tauxy0_-12500.0/v2_sept_perturb_tauxy0_-12500.0.out.nc'];
ice_files{14} = [jobDir,'v2_sept_perturb_tauxy0_-50000/v2_sept_perturb_tauxy0_-50000.out.nc'];
ice_files{15} = [jobDir,'v2_sept_perturb_temp_interior_-15.0/v2_sept_perturb_temp_interior_-15.0.out.nc'];
ice_files{16} = [jobDir,'v2_sept_perturb_temp_interior_-5.0/v2_sept_perturb_temp_interior_-5.0.out.nc'];
ice_files{17} = [jobDir,'v2_sept_perturb_top_-0.85/v2_sept_perturb_top_-0.85.out.nc'];
ice_files{18} = [jobDir,'v2_sept_perturb_top_-1.35/v2_sept_perturb_top_-1.35.out.nc'];
ice_files{19} = [jobDir,'v2_sept_perturb_top_-1.85/v2_sept_perturb_top_-1.85.out.nc'];
ice_files{20} = [jobDir,'v2_sept_perturb_train_time_const_1800.0/v2_sept_perturb_train_time_const_1800.0.out.nc'];
ice_files{21} = [jobDir,'v2_sept_perturb_train_time_const_5400.0/v2_sept_perturb_train_time_const_5400.0.out.nc'];
ice_files{22} = [jobDir,'v2_sept_perturb_visc_15.0/v2_sept_perturb_visc_15.0.out.nc'];

n_cases = length(plume_files);

%for i=1:n_cases
%  i
%  dice(i)   = nc_ice_read(ice_files{i},10,-1);
%  dplume(i) = nc_plume_read(plume_files{i},10,-1);
%end

mgrid = size(dice(1).thk,1);
ngrid = size(dice(1).thk,2);

mean_thk = zeros(mgrid,ngrid);
mean_bmlt = zeros(mgrid,ngrid);

thk_anoms = zeros(mgrid,ngrid,n_cases);
bmlt_anoms = zeros(mgrid,ngrid,n_cases);

for i=1:n_cases
  mean_thk = mean_thk + (1.0/n_cases)*dice(i).thk(:,:,end);
  mean_bmlt = mean_bmlt+(1.0/n_cases)*dice(i).bmlt(:,:,end);
end

for i=1:n_cases
  thk_anoms(:,:,i) = dice(i).thk(:,:,end) - mean_thk;
  bmlt_anoms(:,:,i) = dice(i).bmlt(:,:,end) - mean_bmlt;
end
  
%X_thk = reshape(thk_anoms,[mgrid*ngrid n_cases]);
%[U_thk,S_thk,V_thk] = svd(X_thk);

%X_bmlt = reshape(bmlt_anoms,[mgrid*ngrid n_cases]);
%[U_bmlt,S_bmlt,V_bmlt] = svd(X_bmlt);

%Us_thk = reshape(U_thk(:,1:n_cases),[mgrid ngrid n_cases]);
%Us_bmlt = reshape(U_bmlt(:,1:n_cases),[mgrid ngrid n_cases]);

xs = dice(1).xgrid/1000.0;
ys = dice(1).ygrid/1000.0;

fs = 16;

S_thk_frac = diag(S_thk).*diag(S_thk) / sum(diag(S_thk).*diag(S_thk));
S_bmlt_frac = diag(S_bmlt).*diag(S_bmlt) / sum(diag(S_bmlt).*diag(S_bmlt));

figure(1);
subplot(1,3,1);
contourf(xs,ys,mean_thk',25,'EdgeColor','None');colorbar;
title('mean ice thickness','FontSize',fs);
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);

subplot(1,3,2);
hold on
contourf(xs,ys,Us_thk(:,:,1)','EdgeColor','None');colorbar;
%caxis([min(min(Us_thk(:,:,1))) 0]);
contour(xs,ys,Us_thk(:,:,1)',[0 0]);
title(['First EOF (',sprintf('%3.1f',100*S_thk_frac(1)),' percent of var.)'],'FontSize',fs);
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);

subplot(1,3,3);
hold on
contourf(xs,ys,Us_thk(:,:,2)','EdgeColor','None');colorbar;
contour(xs,ys,Us_thk(:,:,2)',[0 0]);
title(['Second EOF (',sprintf('%3.1f',100*S_thk_frac(2)),' percent of var.)'],'FontSize',fs);
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);

figure(2);
subplot(1,3,1);
contourf(xs,ys,mean_bmlt',25,'EdgeColor','None');colorbar;
title('mean bmlt rate (m/year)','FontSize',fs);
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);

subplot(1,3,2);
hold on
contourf(xs,ys,Us_bmlt(:,:,1)','EdgeColor','None');colorbar;
%caxis([min(min(Us_bmlt(:,:,1))) 0]);
contour(xs,ys,Us_bmlt(:,:,1)',[0 0]);
title(['First EOF (',sprintf('%3.1f',100*S_bmlt_frac(1)),' percent of var.)'],'FontSize',fs);
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);

subplot(1,3,3);
hold on
contourf(xs,ys,Us_bmlt(:,:,2)','EdgeColor','None');colorbar;
contour(xs,ys,Us_bmlt(:,:,2)',[0 0]);
title(['Second EOF (',sprintf('%3.1f',100*S_bmlt_frac(2)),' percent of var.)'],'FontSize',fs);
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);
