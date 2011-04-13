function [] = plot_channels(data)

x = data.x;
y = data.y;
su = data.su;
sv = data.sv;
u = data.u;
v = data.v;
bmelt = data.bmelt;
pdep = data.pdep;
draft = data.draft;
ambdepth = data.ambdepth;

  figure('Units','centimeters','Position',[0 0 20 20]);

fs = 14;

[m,n] = size(draft);
channel_amp = zeros(m,n);
mean_draft = zeros(m,n);
max_draft = zeros(m,n);

for i=1:n
  channel_amp(:,i) = draft(:,i) - min(draft(:,i));
  mean_draft(:,i) = sum(draft(:,i))/m;
  max_draft(:,i) = min(draft(:,i));
end

subplot(1,1,1);
  hold on
  contourf(x/1000.0,y/1000.0,channel_amp',20,'EdgeColor','None') ;colorbar;
%  contour(x/1000.0,y/1000.0,mean_draft',20,'w')
  contour(x/1000.0,y/1000.0,max_draft',20,'w')
  caxis([min(min(channel_amp)) max(max(channel_amp))]);
  title('Channel depth with deepest-draft contours','FontSize',fs);

end 
