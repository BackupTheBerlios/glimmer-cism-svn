%%%% plot plume cross section 


fig1 = figure(1);
set(fig1,'Position',[1 1 800 600]);
clf;

lw = 2.0;
fs = 16;


x = (dplume_avg.x - 1375)/1000.0;

yslice0 = 1;
yslices = [0,8,16,24,32,14*4,40*4]+yslice0;

hold on
for i=1:length(yslices)
    yslice = yslices(i);
    plot(x,dplume_avg.bpos(:,yslice)-800,'b-','LineWidth',lw);
plot(x,dplume_avg.bpos(:,yslice)-dplume_avg.pdep(:,yslice)-800,'r-','LineWidth',lw);

end
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Depth (m)','FontSize',fs);


set(gca,'FontSize',fs);
