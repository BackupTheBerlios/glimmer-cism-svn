function [] = plot_bmlt_draft_overlay(data)

x = data.x;
y = data.y;
bmelt = data.bmelt;
pdep = data.pdep;
draft = data.draft;
ambdepth = data.ambdepth;
train = data.train;

%figure('Units','centimeters','Position',[0 0 40 20]);
fs = 9;

%subplot(1,2,1);
hold on
contourf(x/1000.0,y/1000.0,bmelt',20,'EdgeColor','None') ;colorbar;
C =contour(x/1000.0,y/1000.0,draft',-600:25:0,'w');
caxis([0, max(max(bmelt))]);
title('Basal melt rate (m/y) with ice draft contours','FontSize',fs);
xlabel('Cross shelf distance (km)','FontSize',14);
ylabel('Along shelf distance (km)','FontSize',14);
%clabel(C,'FontSize',14,'Color','w');
title('Basal melt rate (m/year) with ice draft contours','FontSize',14);

%subplot(1,2,2);
%hold on
%contourf(data.x/1000.0,data.y/1000.0,train',...
%         20,'EdgeColor','None') ;colorbar;
%contour(data.x/1000.0,data.y/1000.0,train',[0 0],'k')
%contour(data.x/1000.0,data.y/1000.0,data.draft',20,'w');
%caxis([min(min(data.train)), max(max(train))]);
%title('train rate (m/year) with draft contours','FontSize',fs);
%hold off

end 
