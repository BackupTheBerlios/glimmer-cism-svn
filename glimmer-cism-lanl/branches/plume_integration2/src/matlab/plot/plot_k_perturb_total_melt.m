clear all;

jobs = getenv('GC_JOBS');

usq2_amp25_ks = [0,1,2,3,4,5,6,7,...
		  8,9,10,11,12,13,14,15];
usq2_amp50_ks = [0,1,2,3,4,5, ...
		  6,7,8,9,10,11,12,13,14,15];

usq3_amp25_ks = [6,8,9,10];
usq3_amp50_ks = [7,8,9,10,11,12,13,14,15];

usq2_amp25_jobs = cell(length(usq2_amp25_ks),1);
usq2_amp50_jobs = cell(length(usq2_amp50_ks),1);
usq3_amp25_jobs = cell(length(usq3_amp25_ks),1);
usq3_amp50_jobs = cell(length(usq3_amp50_ks),1);


m_2_25 = zeros(length(usq2_amp25_ks),1);
m_2_50 = zeros(length(usq2_amp50_ks),1);
m_3_25 = zeros(length(usq3_amp25_ks),1);
m_3_50 = zeros(length(usq3_amp50_ks),1);

u_2_25 = zeros(length(usq2_amp25_ks),1);
u_2_50 = zeros(length(usq2_amp50_ks),1);
u_3_25 = zeros(length(usq3_amp25_ks),1);
u_3_50 = zeros(length(usq3_amp50_ks),1);

for i=1:length(usq2_amp25_ks)
  usq2_amp25_jobs{i} = strcat(['oct30_perturb_usq.2_thk_perturb_k_',sprintf('%.1f',usq2_amp25_ks(i)),'_amp_25.0']);
jname = usq2_amp25_jobs{i};
f_ice = strcat([jobs,'/',jname,'/',jname,'.out.nc']);
f_ocean = strcat([jobs,'/',jname,'/plume.',jname,'.out.nc']);
istart = -1;
iend = -1;
stride = 1;
dice = nc_ice_read(f_ice,  istart,stride,iend);
[melt,f_in,f_out,acab,unsteady] = mass_balance(dice);
m_2_25(i) = melt;
u_2_25(i) = unsteady;
end

for i=1:length(usq2_amp50_ks)
  usq2_amp50_jobs{i} = strcat(['oct30_perturb_usq.2_thk_perturb_k_',sprintf('%.1f',usq2_amp50_ks(i)),'_amp_50.0']);
jname = usq2_amp50_jobs{i};
f_ice = strcat([jobs,'/',jname,'/',jname,'.out.nc']);
f_ocean = strcat([jobs,'/',jname,'/plume.',jname,'.out.nc']);
istart = -1;
iend = -1;
stride = 1;
dice = nc_ice_read(f_ice,  istart,stride,iend);
[melt,f_in,f_out,acab,unsteady] = mass_balance(dice);
m_2_50(i) = melt;
u_2_50(i) = unsteady;
end

for i=1:length(usq3_amp25_ks)
  usq3_amp25_jobs{i} = strcat(['oct30_perturb_usq.3_thk_perturb_k_',sprintf('%.1f',usq3_amp25_ks(i)),'_amp_25.0']);
jname = usq3_amp25_jobs{i};
f_ice = strcat([jobs,'/',jname,'/',jname,'.out.nc']);
f_ocean = strcat([jobs,'/',jname,'/plume.',jname,'.out.nc']);
istart = -1;
iend = -1;
stride = 1;
dice = nc_ice_read(f_ice,  istart,stride,iend);
[melt,f_in,f_out,acab,unsteady] = mass_balance(dice);
m_3_25(i) = melt;
u_3_25(i) = unsteady;
end

for i=1:length(usq3_amp50_ks)
  usq3_amp50_jobs{i} = strcat(['oct30_perturb_usq.3_thk_perturb_k_',sprintf('%.1f',usq3_amp50_ks(i)),'_amp_50.0']);
jname = usq3_amp50_jobs{i};
f_ice = strcat([jobs,'/',jname,'/',jname,'.out.nc']);
f_ocean = strcat([jobs,'/',jname,'/plume.',jname,'.out.nc']);
istart = -1;
iend = -1;
stride = 1;
dice = nc_ice_read(f_ice,  istart,stride,iend);
[melt,f_in,f_out,acab,unsteady] = mass_balance(dice);
m_3_50(i) = melt;
u_3_50(i) = unsteady;
end



figure(1);
clf;
ms = 15.0;
fs = 15;

ax1Color = 'b';
ax2Color = 'r';
amp25Mark = '.';
amp50Mark = 'x';

hl1 = line(usq2_amp25_ks,100*m_2_25/f_in);
ax1 = gca;

ax2 = axes('Position',get(ax1,'Position'),...
	   'XAxisLocation','bottom',...
           'YAxislocation','right',...
            'Color','None',...
            'XColor','k','YColor',ax2Color);
%hold(ax2);
%hold(ax1);

hl2 = line(usq2_amp25_ks,100*u_2_25/f_in);

set(ax1,'XColor','k','Ycolor',ax1Color);
set(ax2,'XColor','k','Ycolor',ax2Color);

set(hl1,'MarkerSize',ms,'LineStyle',amp25Mark,'Color',ax1Color);
set(hl2,'MarkerSize',ms,'LineStyle',amp25Mark,'Color',ax2Color);

set(ax2,'Ydir','reverse');
set(ax2,'ylim',[-50,5]);
set(ax1,'ylim',[50,105]);
set(ax1,'xlim',[-1 16]);
set(ax2,'xlim',[-1 16]);

axes(ax1);
hl3 = line(usq2_amp50_ks,100*m_2_50/f_in);
  set(hl3,'MarkerSize',ms,'LineStyle',amp50Mark,'Color',ax1Color);
axes(ax2);
hl4 = line(usq2_amp50_ks,100*u_2_50/f_in);
set(hl4,'MarkerSize',ms,'LineStyle',amp50Mark,'Color',ax2Color);

%plot(usq3_amp25_ks,100*m_3_25/f_in,'b+','MarkerSize',ms);
%plot(usq3_amp25_ks,100*u_3_25/f_in,'c+','MarkerSize',ms);

%plot(usq3_amp50_ks,100*m_3_50/f_in,'r+','MarkerSize',ms);
%plot(usq3_amp50_ks,100*u_3_50/f_in,'m+','MarkerSize',ms);

set(get(ax1,'ylabel'),'String','Total basal melting (% of influx)');
set(get(ax1,'ylabel'),'FontSize',fs);
set(get(ax2,'ylabel'),'String','Volume rate of change (% of influx)');
set(get(ax2,'ylabel'),'FontSize',fs);

set(get(ax1,'xlabel'),'String','Cross-shelf perturbation wavenumber (1/km)');
set(get(ax1,'xlabel'),'FontSize',fs);

hold off
