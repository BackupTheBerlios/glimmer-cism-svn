%%% plot plume temp and salt at depths

flatsalt = flatten_field(dplume_avg.salt);
flattemp = flatten_field(dplume_avg.temp);
flatdraft = flatten_field(dplume_avg.bpos-800.0);
flatmiddepth = flatten_field(dplume_avg.bpos-0.5*dplume_avg.pdep-800.0);
uspeed = pad_edge(1,0,dplume_avg.su(2:end-1,:,:)+dplume_avg.su(1:end-2,:,:))/0.5;
vspeed = pad_edge(0,1,dplume_avg.sv(:,2:end-1,:)+dplume_avg.sv(:,1:end-2,:))/0.5;
flatspeed = flatten_field( sqrt( uspeed.^2+vspeed.^2));
       
fs = 14;
fs_label = 18;

fig1 = figure(1);
set(fig1,'Position',[1 1 1200 800]);
clf;

subplot(1,2,1);
hold on
plot(flattemp,flatmiddepth,'r.');
plot(amb_t(abs(flatmiddepth)),flatmiddepth,'k.');
xlabel('Temperature (^\circ C)','FontSize',fs);
ylabel('Depth (m)','FontSize',fs);
set(gca,'FontSize',fs);
xlim([-2 0.15]);
text(-2.25,15,'a','color','k','FontSize',fs_label)
hold off

subplot(1,2,2);
hold on;
plot(flatsalt,flatmiddepth,'b.');
plot(amb_s(abs(flatmiddepth)),flatmiddepth,'k.');
xlabel('Salinity (psu)','FontSize',fs);
ylabel('Depth (m)','FontSize',fs);
set(gca,'FontSize',fs);
xlim([34.6 34.76]);
text(34.5825,15,'b','color','k','FontSize',fs_label)
hold off;
