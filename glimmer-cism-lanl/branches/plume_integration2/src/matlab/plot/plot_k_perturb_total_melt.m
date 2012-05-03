%clear all;

recalculate = true;
%recalculate = false;

if (recalculate)
  clear all;
jobs = getenv('GC_JOBS');

baseName = 'oct30_perturb_usq.9';

usq_amp25_ks = [0,1,2,3,4,5,6,7,...
		  8,9,10,11,12,13,14,15];
usq_amp50_ks = [0,1,2,3,4,5, ...
		  6,7,8,9,10,11,12,13,14,15];

%usq_amp25_ks = [1,2,3,4,5,6,7,...
%		 8,9,10,11,12,13,14,15];
%usq_amp50_ks = [2, ...
%		  6,7,8,9,10,11,12,13,14,15];

usq_amp25_jobs = cell(length(usq_amp25_ks),1);
usq_amp50_jobs = cell(length(usq_amp50_ks),1);


m_25 = zeros(length(usq_amp25_ks),1);
m_50 = zeros(length(usq_amp50_ks),1);

u_25 = zeros(length(usq_amp25_ks),1);
u_50 = zeros(length(usq_amp50_ks),1);

for i=1:length(usq_amp25_ks)
  usq_amp25_jobs{i} = strcat([baseName,'_thk_perturb_k_',sprintf('%.1f',usq_amp25_ks(i)),'_amp_25.0']);
jname = usq_amp25_jobs{i}
f_ice = strcat([jobs,'/',jname,'/',jname,'.out.nc']);
f_ocean = strcat([jobs,'/',jname,'/plume.',jname,'.out.nc']);
istart = -1;
iend = -1;
stride = 1;
dice = nc_ice_read(f_ice,  istart,stride,iend);
dplume = nc_plume_read(f_ocean, istart,stride,iend);

[melt,melt_applied,f_in,f_out,acab,unsteady] = mass_balance(dice,dplume,-1);

m_25(i) = melt;
ma_25(i) = melt_applied;
u_25(i) = unsteady;

end

for i=1:length(usq_amp50_ks)
  usq_amp50_jobs{i} = strcat([baseName,'_thk_perturb_k_',sprintf('%.1f',usq_amp50_ks(i)),'_amp_50.0']);
jname = usq_amp50_jobs{i}
f_ice = strcat([jobs,'/',jname,'/',jname,'.out.nc']);
f_ocean = strcat([jobs,'/',jname,'/plume.',jname,'.out.nc']);
istart = -1;
iend = -1;
stride = 1;
dice = nc_ice_read(f_ice,  istart,stride,iend);
dplume = nc_plume_read(f_ocean, istart,stride,iend);
[melt,melt_applied,f_in,f_out,acab,unsteady] = mass_balance(dice,dplume,-1);

m_50(i) = melt;
ma_50(i) = melt_applied;
u_50(i) = unsteady;
end

end
fig1 = figure(1);
set(fig1,'Position',[1 1 800 600]);
clf;
ms = 15.0;
fs = 18;
fs = 18;

ax1Color = 'b';
ax2Color = 'k';
amp25Mark = '-';
amp50Mark = 'x';

which_ks = 1:5;
which_ks = 1:16;
melt_color = 'b';
hl1 = line(usq_amp25_ks(which_ks),100*m_25(which_ks)/f_in);
%hl1a = line(usq_amp25_ks,100*ma_25/f_in);

ax1 = gca;


ax2 = axes('Position',get(ax1,'Position'),...
	   'XAxisLocation','bottom',...
       'YAxislocation','right',...
       'Color','None',...
	   'XColor','k','YColor',ax2Color);

set(ax1,'FontSize',fs);
set(ax2,'FontSize',fs);

hl2 = line(usq_amp25_ks,100*u_25/f_in);

set(ax1,'XColor','k','Ycolor',ax1Color);
set(ax2,'XColor','k','Ycolor',ax2Color);
lw = 2.0;
set(hl1,'MarkerSize',ms,'LineStyle',amp25Mark,'LineWidth',lw,'Color',melt_color);
%set(hl1a,'MarkerSize',ms,'LineStyle',amp25Mark,'LineWidth',lw,'Color','b');
set(hl2,'MarkerSize',ms,'LineStyle',amp25Mark,'LineWidth',lw,'Color',ax2Color);

set(ax2,'Ydir','reverse');
set(ax2,'ylim',[-25,2]);
set(ax2,'ytick',-6:2:2);
set(ax1,'ylim',[55,120]);
set(ax1,'ytick',70:10:120);
set(ax1,'xlim',[-1 16]);
set(ax2,'xlim',[-1 16]);

axes(ax1);
hl3 = line(usq_amp50_ks(which_ks),100*m_50(which_ks)/f_in);
%hl3a = line(usq_amp50_ks,100*ma_50/f_in);

set(hl3,'MarkerSize',ms,'LineStyle',amp50Mark,'Color',melt_color);
%set(hl3a,'MarkerSize',ms,'LineStyle',amp50Mark,'Color','b');

axes(ax2);
hl4 = line(usq_amp50_ks,100*u_50/f_in);
set(hl4,'MarkerSize',ms,'LineStyle',amp50Mark,'Color',ax2Color);


set(get(ax1,'ylabel'),'String','Total basal melting (% of influx)');
set(get(ax1,'ylabel'),'FontSize',fs);
set(get(ax2,'ylabel'),'String','Volume rate of change (% of influx)');
set(get(ax2,'ylabel'),'FontSize',fs);

set(get(ax1,'xlabel'),'String','Cross-shelf perturbation wavenumber');
set(get(ax1,'xlabel'),'FontSize',fs);

hold off

%legend([hl1,hl1a,hl2,hl3,hl3a,hl4], ...
legend([hl1,hl2,hl3,hl4], ...
       '25 m - m''',... %       '25 m - \gamma m''',...
        '25 m - unsteadiness',...
        '50 m - m''',... %       '50 m - \gamma m''', ...
       '50 m - unsteadiness');
