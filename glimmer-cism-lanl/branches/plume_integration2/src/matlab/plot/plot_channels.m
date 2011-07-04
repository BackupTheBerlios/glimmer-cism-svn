function [] = plot_channels(data, timeslice)

[m,n,k] = size(data.bmelt);
if (timeslice < 0) 
    timeslice = k;
end

x = data.x;
y = data.y;
su = data.su(:,:,timeslice);
sv = data.sv(:,:,timeslice);
u = data.u(:,:,timeslice);
v = data.v(:,:,timeslice);
bmelt = data.bmelt(:,:,timeslice);
pdep = data.pdep(:,:,timeslice);
draft = data.draft(:,:,timeslice);
ambdepth = data.ambdepth(:,:,timeslice);

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

%subplot(1,1,1);
  hold on
  contourf(x/1000.0,y/1000.0,channel_amp',20,'EdgeColor','None') ;colorbar;
%  contour(x/1000.0,y/1000.0,mean_draft',20,'w')
  [C,h] = contour(x/1000.0,y/1000.0,max_draft',10,'w');
  clabel(C,h,'manual');
  caxis([min(min(channel_amp)) max(max(channel_amp))]);
  title('Channel depth with deepest-draft contours','FontSize',fs);

end 
