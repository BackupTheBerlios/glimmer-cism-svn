jname = 'oct30_perturb_usq.2_entype2_5_ci_10.0_cn_10.0_cs_10.0';
jname = 'oct30_perturb_usq.2_entype2_7_m_5.0';
jname = 'oct30_perturb_usq.2_visc_100.0';

f = sprintf('%s/%s/%s.out.nc',getenv('GC_JOBS'),jname,jname)
fp = sprintf('%s/%s/plume.%s.out.nc',getenv('GC_JOBS'),jname,jname)

  dice = nc_ice_read(f,-1,1,-1);
  dplume = nc_plume_read(fp,-1,1,-1);

xs = dice.Xgrid/1000;
ys = dice.Ygrid/1000;

draft = dice.lsurf(:,:,end);

contourf(xs,ys,-draft,30,'EdgeColor','None');colorbar;

fs = 16;
xlabel('Across shelf distance (km)','FontSize',fs)
ylabel('Along shelf distance (km)','FontSize',fs)
  set(gca,'FontSize',fs);
  title('Ice draft (m)','FontSize',fs);


