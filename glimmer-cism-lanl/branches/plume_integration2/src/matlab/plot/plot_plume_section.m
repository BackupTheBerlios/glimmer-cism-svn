%%%% plot plume cross section 


figure(1);
clf;

lw = 2.0;
fs = 14;



x = (dplume_avg.x - 1375)/1000.0;


yslice0 = 1;
yslices = [0,8,16,24,32,160]+yslice0;

%subplot(2,1,1);
hold on
for i=1:length(yslices)
    yslice = yslices(i);
    plot(x,dplume_avg.bpos(:,yslice)-800,'b-','LineWidth',lw);
plot(x,dplume_avg.bpos(:,yslice)-dplume_avg.pdep(:,yslice)-800,'r-','LineWidth',lw);

end
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Depth (m)','FontSize',fs);

legend('Ice base','Plume lower surface','Location','SouthEast');

title('Sections across plume at 1km,3km,5km,7km,9km and 40km downstream','FontSize',fs);

if (false)
subplot(2,1,2);
yslices = [21];
for i=1:length(yslices)
    yslice = yslices(i);
    plotyy(x,dplume_avg.bpos(:,yslice)-800, ...
           x,res.g_prime(:,yslice));

end
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Depth (m)','FontSize',fs);

legend('Ice base','reduced gravity','Location','SouthEast');

title('Sections across plume at 1km,3km,5km and 7km downstream','FontSize',fs);
end