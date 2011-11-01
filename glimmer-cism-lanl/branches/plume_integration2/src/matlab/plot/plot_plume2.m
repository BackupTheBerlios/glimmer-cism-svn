function [] = plot_plume2(dplume,dice,timeslice,flabel)

fig_dir = '/home/cvg222/matlab_figs/';

fs = 18;
fs3 = 14;
fs6 = 12;
fs2 = 14;

solo_fig_x_size = 800;
solo_fig_y_size = 800;
three_panel_fig_x_size = 1250;
three_panel_fig_y_size = 500;
four_panel_fig_x_size = 1200;
four_panel_fig_y_size = 800;
six_panel_fig_x_size = 1500;
six_panel_fig_y_size = 1000;


thermal_forcing = @(depth,salt,temp) temp-(-5.73e-2*salt-7.61d-4*depth+8.32d-2);

if (timeslice < 0) 
    timeslice = size(dplume.su, 3);
end 

x = dice.xgrid;
y = dice.ygrid;

ice_u = squeeze(dice.uvelmean(:,:,timeslice));
ice_v = squeeze(dice.vvelmean(:,:,timeslice));

su = squeeze(dplume.su(:,:,timeslice));
sv = squeeze(dplume.sv(:,:,timeslice));
u = squeeze(dplume.u(:,:,timeslice));
v = squeeze(dplume.v(:,:,timeslice));
train = squeeze(dplume.train(:,:,timeslice));
bmelt = squeeze(dplume.bmelt(:,:,timeslice));
temp = squeeze(dplume.temp(:,:,timeslice));
salt = squeeze(dplume.salt(:,:,timeslice));
pdep = squeeze(dplume.pdep(:,:,timeslice));
draft = squeeze(dplume.draft(:,:,timeslice));
thk_t = squeeze(dice.thk_t(:,:,timeslice));
T_forcing = thermal_forcing(abs(draft),salt,temp);
grad = squeeze(dplume.grad(:,:,timeslice));
speed = sqrt(su.*su + sv.*sv);

j_start = 3;
j_stop = 118;

ice_def = - dice.flux_div(:,j_start:j_stop,timeslice) + ...
            dice.y_adv(:,j_start:j_stop,timeslice);
y_adv = - dice.y_adv(:,j_start:j_stop,timeslice);
ice_adv = -dice.flux_div(:,j_start:j_stop,timeslice) + ...
           dice.vel_div(:,j_start:j_stop,timeslice);   

[m,n] = size(draft);
channel_amp = zeros(m,n);
mean_draft = zeros(m,n);
max_draft = zeros(m,n);

for i=1:n
  channel_amp(:,i) = draft(:,i) - min(draft(:,i));
  mean_draft(:,i) = sum(draft(:,i))/m;
  max_draft(:,i) = min(draft(:,i));
end

stride = 2;
stride2 = 4;

scale = 2.0;
lw = 0.5;

switch1 = false;
switch2 = false;
switch3= true;

if (switch1) 
  plot_solo_panels;
end

if (switch2)
plot_3_panel_ice_geometry;
end

if (switch2)
plot_4_panel_melting;
end

if (switch2)
plot_draft_and_plumevel;
end

if (switch3)
plot_6_panel_plume;
end

if (switch2)
plot_cross_shelf_avg;
end

end
