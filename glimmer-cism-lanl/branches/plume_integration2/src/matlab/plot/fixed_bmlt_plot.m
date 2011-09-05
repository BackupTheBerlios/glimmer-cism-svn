%fixed bmlt plot

jname = '/data/gladish/gc_output/paper_jobs/may30_fixed_bmlt_temp_-20.0_tau_25000.0/may30_fixed_bmlt_temp_-20.0_tau_25000.0.out.nc';

%data = nc_ice_read(jname,10, 200);

fs = 16;
lw = 3.0;

figure(1);
clf;
bmlt = data.bmlt(1,:,end);
xmid = floor((length(data.xgrid)+1)/2);
xmid2 = ceil((length(data.xgrid)+1)/2);
%center_thk = 0.5*(data.thk(xmid,:,end)+data.thk(xmid2,:,end));
%[ax,h1,h2] = plotyy(data.ygrid/1000.0,bmlt,data.ygrid/1000.0,center_thk);
center_draft = 0.5*(data.lsurf(xmid,:,end)+data.lsurf(xmid2,:,end));
[ax,h1,h2] = plotyy(data.ygrid/1000.0,bmlt,data.ygrid/1000.0,center_draft);
set(get(ax(1),'Ylabel'),'String','basal melt rate (m/a)','FontSize',fs)
set(get(ax(2),'Ylabel'),'String','ice draft (m)','FontSize',fs)
set(h1,'Linewidth',lw,'Linestyle','-');
set(h2,'Linewidth',lw,'Linestyle','-');
set(ax(1),'xlim',[-2 55]);
set(ax(2),'xlim',[-2 55]);
set(ax(1),'ylim',[0 27]);
set(ax(2),'ylim',[-550 0]);
set(ax(1),'ytick',0:5:25);
set(ax(2),'ytick',-500:100:0);
set(ax(1),'FontSize',fs);
set(ax(2),'FontSize',fs);

xlabel('along shelf distance (km)','FontSize',fs)

figure(2);
clf;
ice_deformation = -data.flux_div + data.y_adv;
ice_deformation(:,1:2,:) = 0.0;
ice_deformation(:,end-1:end,:) = 0.0;

subplot(1,3,1);
contourf(data.Xgrid/1000,data.Ygrid/1000,ice_deformation(:,:,end),30,'EdgeColor','None');colorbar;
caxis([-5 1]);
subplot(1,3,2);
contourf(data.Xgrid/1000,data.Ygrid/1000,-data.y_adv(:,:,end),30,'EdgeColor','None');colorbar;
caxis([0 25]);
subplot(1,3,3);
contourf(data.Xgrid/1000,data.Ygrid/1000,data.bmlt(:,:,end),30,'EdgeColor','None');colorbar;
caxis([0 25]);

figure(3);
clf;
center_ice_deformation = 0.5*(ice_deformation(xmid,:,end)+ice_deformation(xmid2,:,end));
center_y_adv = 0.5*(-data.y_adv(xmid,:,end)-data.y_adv(xmid2,:,end));

%[ax,h1,h2] = plotyy(data.ygrid/1000.0,bmlt,data.ygrid/1000,center_ice_deformation);
hold on
ybegin = 3;
ystop = 1;
h = plot(data.ygrid(ybegin:end-ystop)/1000.0,bmlt(ybegin:end-ystop),'k');
set(h,'linestyle','-','linewidth',lw);
h = plot(data.ygrid(ybegin:end-ystop)/1000,center_y_adv(ybegin:end-ystop),'r');
set(h,'linestyle','-','linewidth',lw);
h = plot(data.ygrid(ybegin:end-ystop)/1000,center_ice_deformation(ybegin:end-ystop),'b');
set(h,'linestyle','-','linewidth',lw);
set(gca,'fontsize',fs);
plot([0 70],[0 0],'k')
xlim([-1 55])
%set(h1,'Linewidth',lw,'Linestyle','-');
%set(h2,'Linewidth',lw,'Linestyle','-');
%set(h3,'Linewidth',lw,'Linestyle','-');
xlabel('along shelf distance (km)','FontSize',fs);
ylabel('ice thickness tendency (m/a)','FontSize',fs);
%legend('melt rate', ...
%       '-v_0 \partial_y H',...
%       '-\partial_x (uH) - \partial_y((v-v_0)H)'),legend('boxoff');
xlabel('along shelf distance (km)','FontSize',fs)

