h_diff_t = dice.H_diff_t(:,:,end);
lsurf = dice.lsurf(:,:,end);
bmlt = dice.bmlt(:,:,end);
thk = dice.thk(:,:,end);
uflx_conv = dice.uflx_conv(:,:,end);

x = dice.Xgrid/1000;
y = dice.Ygrid/1000;

%contourf(x,y,h_diff_t);colorbar;

flatlsurf = flatten_field(lsurf);
flat_h_diff = flatten_field(h_diff_t);
flatbmlt = flatten_field(bmlt);
flatthk = flatten_field(thk);
flatuflxconv = flatten_field(uflx_conv);

fig1 = figure(1);
clf;
hold on
plot(flat_h_diff,flatbmlt,'k.');
lw = 2.0;
plot([0 10],[0 10],'r-','LineWidth',lw);
plot([0 -10],[0 10],'r-','LineWidth',lw);
plot([0 10],[0 20],'g-','LineWidth',lw);
plot([0 -10],[0 20],'g-','LineWidth',lw);
xlim([-10 10]);
ylim([-5 75]);
fs = 14;
xlabel('thickening from artificial ice diffusion (m/a)','FontSize',fs);
ylabel('basal melt rate (m/a)','FontSize',fs);

figure(2);
hold on;
plot(flatuflxconv,flat_h_diff,'k.');
plot([-10 10],[-10 10],'r-','LineWidth',lw);
plot([10 -10],[-10 10],'r-','LineWidth',lw);
hold off;
