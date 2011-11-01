function [] = plot_draft_sections(ice_data,y_sections)

k = length(y_sections);

x = ice_data.xgrid;
y = ice_data.ygrid;
tend = size(ice_data.lsurf,3);

draft = ice_data.lsurf(:,:,tend);
colors='krbgm';
fs = 18;
fs2 = 20;

%figure(1);
%clf;

hold on

%locs = []

for i=1:k
    colr = colors(mod(i-1,length(colors))+1);
    plot(x/1000.0,draft(:,y_sections(i)),[colr,'-'],'LineWidth',2.0);
    set(gca,'FontSize',fs);
    %locs = [locs,sprintf('%3.1f,',y(y_sections(i)))];
text(21,-650+i*90,strcat([sprintf('%4.1f',(y(y_sections(i))-100.0)/1000.0),' km']),'FontSize',fs,'Color',colr);
end
xlabel('Cross-shelf distance (km)','FontSize',fs);
ylabel('Ice Draft (m)','FontSize',fs);
title('Cross-shelf sections of ice draft','FontSize',fs);

%text(19,-600+7.25*50,'distance from gr. ln.','FontSize',16,'Color','k');

%clf;
%contourf(x/1000.0,y/1000.0,draft',30);colorbar;
end
