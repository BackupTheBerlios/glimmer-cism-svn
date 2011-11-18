clear all;

jobs = getenv('GC_JOBS');

baseName = 'oct30_perturb_usq.8';

usq_amp25_ks = [0,1,2,3,4,5,6,7,...
		  8,9,10,11,12,13,14,15];
usq_amp50_ks = [0,1,2,3,4,5, ...
		  6,7,8,9,10,11,12,13,14,15];

usq_amp25_jobs = cell(length(usq_amp25_ks),1);
usq_amp50_jobs = cell(length(usq_amp50_ks),1);


m_25 = zeros(length(usq_amp25_ks),1);
m_50 = zeros(length(usq_amp50_ks),1);

u_25 = zeros(length(usq_amp25_ks),1);
u_50 = zeros(length(usq_amp50_ks),1);

for i=1:length(usq_amp25_ks)
  usq_amp25_jobs{i} = strcat([baseName,'_thk_perturb_k_',sprintf('%.1f',usq_amp25_ks(i)),'_amp_25.0']);
jname = usq_amp25_jobs{i};
f_ice = strcat([jobs,'/',jname,'/',jname,'.out.nc']);
f_ocean = strcat([jobs,'/',jname,'/plume.',jname,'.out.nc']);
istart = -1;
iend = -1;
stride = 1;
dice = nc_ice_read(f_ice,  istart,stride,iend);
[melt,f_in,f_out,acab,unsteady] = mass_balance(dice);
m_25(i) = melt;
u_25(i) = unsteady;
end

for i=1:length(usq_amp50_ks)
  usq_amp50_jobs{i} = strcat([baseName,'_thk_perturb_k_',sprintf('%.1f',usq_amp50_ks(i)),'_amp_50.0']);
jname = usq_amp50_jobs{i};
f_ice = strcat([jobs,'/',jname,'/',jname,'.out.nc']);
f_ocean = strcat([jobs,'/',jname,'/plume.',jname,'.out.nc']);
istart = -1;
iend = -1;
stride = 1;
dice = nc_ice_read(f_ice,  istart,stride,iend);
[melt,f_in,f_out,acab,unsteady] = mass_balance(dice);
m_50(i) = melt;
u_50(i) = unsteady;
end

figure(1);
clf;
ms = 15.0;
fs = 15;

ax1Color = 'b';
ax2Color = 'r';
amp25Mark = '.';
amp50Mark = 'x';

hl1 = line(usq_amp25_ks,100*m_25/f_in);
ax1 = gca;

ax2 = axes('Position',get(ax1,'Position'),...
	   'XAxisLocation','bottom',...
           'YAxislocation','right',...
            'Color','None',...
            'XColor','k','YColor',ax2Color);
%hold(ax2);
%hold(ax1);

hl2 = line(usq_amp25_ks,100*u_25/f_in);

set(ax1,'XColor','k','Ycolor',ax1Color);
set(ax2,'XColor','k','Ycolor',ax2Color);

set(hl1,'MarkerSize',ms,'LineStyle',amp25Mark,'Color',ax1Color);
set(hl2,'MarkerSize',ms,'LineStyle',amp25Mark,'Color',ax2Color);

set(ax2,'Ydir','reverse');
set(ax2,'ylim',[-50,10]);
set(ax1,'ylim',[50,105]);
set(ax1,'xlim',[-1 16]);
set(ax2,'xlim',[-1 16]);

axes(ax1);
hl3 = line(usq_amp50_ks,100*m_50/f_in);
  set(hl3,'MarkerSize',ms,'LineStyle',amp50Mark,'Color',ax1Color);
axes(ax2);
hl4 = line(usq_amp50_ks,100*u_50/f_in);
set(hl4,'MarkerSize',ms,'LineStyle',amp50Mark,'Color',ax2Color);


set(get(ax1,'ylabel'),'String','Total basal melting (% of influx)');
set(get(ax1,'ylabel'),'FontSize',fs);
set(get(ax2,'ylabel'),'String','Volume rate of change (% of influx)');
set(get(ax2,'ylabel'),'FontSize',fs);

set(get(ax1,'xlabel'),'String','Cross-shelf perturbation wavenumber (1/km)');
set(get(ax1,'xlabel'),'FontSize',fs);

hold off
