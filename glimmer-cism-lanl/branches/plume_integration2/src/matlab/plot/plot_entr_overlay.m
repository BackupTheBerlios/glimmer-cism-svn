function [] = plot_channels(data)

x = data.x;
y = data.y;
su = data.su;
sv = data.sv;
u = data.u;
v = data.v;
entr = data.entr;
pdep = data.pdep;
draft = data.draft;
ambdepth = data.ambdepth;

  figure('Units','centimeters','Position',[0 0 20 20]);

fs = 9;

  hold on
  contourf(x/1000.0,y/1000.0,entr',20,'EdgeColor','None') ;colorbar;
  contour(x/1000.0,y/1000.0,draft',20,'w')
  caxis([min(min(entr)) max(max(entr))]);
  title('entrainment rate with draft contours','FontSize',fs);

end 
