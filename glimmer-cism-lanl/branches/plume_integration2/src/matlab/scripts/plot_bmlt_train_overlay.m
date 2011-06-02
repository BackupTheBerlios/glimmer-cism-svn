function [] = plot_bmlt_train_overlay(data)

x = data.x;
y = data.y;
bmelt = data.bmelt;
pdep = data.pdep;
draft = data.draft;
ambdepth = data.ambdepth;
train = data.train;

figure('Units','centimeters','Position',[0 0 40 20]);
fs = 9;

subplot(1,2,1);
hold on
contourf(x/1000.0,y/1000.0,bmelt',20,'EdgeColor','None') ;colorbar;
contour(x/1000.0,y/1000.0,draft',20,'w');
caxis([0, max(max(bmelt))]);
title('Basal melt rate (m/y) with draft contours','FontSize',fs);

subplot(1,2,2);
hold on
contourf(data.x/1000.0,data.y/1000.0,train',...
         20,'EdgeColor','None') ;colorbar;
contour(data.x/1000.0,data.y/1000.0,train',[0 0],'k')
contour(data.x/1000.0,data.y/1000.0,data.draft',20,'w');
caxis([min(min(data.train)), max(max(train))]);
title('train rate (m/year) with draft contours','FontSize',fs);
hold off

end 
