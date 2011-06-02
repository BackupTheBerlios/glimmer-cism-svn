function [] = plot_interface_pdep_overlay(data)

  p_data = data;
  figure('Units','centimeters','Position',[0 0 40 20]);

fs = 14;
subplot(1,2,1);
hold on
contourf(p_data.x/1000.0,p_data.y/1000.0,p_data.ambdepth', ...
        20,'EdgeColor','None') ;colorbar;
contour(p_data.x/1000.0,p_data.y/1000.0,p_data.draft',20,'w')
caxis([min(min(p_data.ambdepth)) max(max(p_data.ambdepth))]);
title('Ambient interface depth with ice draft contours','FontSize',fs);

subplot(1,2,2);
hold on
contourf(p_data.x/1000.0,p_data.y/1000.0,p_data.pdep',...
        20,'EdgeColor','None') ;colorbar;
contour(p_data.x/1000.0,p_data.y/1000.0,p_data.draft',20,'w')
caxis([0 max(max(p_data.pdep))]);
title('Plume depth with ice draft contours','FontSize',fs);
end 
