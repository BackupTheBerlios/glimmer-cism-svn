% read in netcdf model output

% To use this script you need:

% to know where the netcdf files containing the model output are:

model_prefix = '/home/gladish/research/gcp_resources/nc_files/last_time.ssj1_';
model_suffix = '_sept16.1.out.nc';

% to know where the quasi-exact solution .mat files are:
exact_prefix = '/home/gladish/research/gcp_resources/mat_files/ssj1_';
exact_suffix = '_exact.sept16.1000m.mat';

tau = [0,10,25,50];
%tau = [0];
acab = [-2.0, 0.0, 2.0];
acab = [0.0];

markers = 'brgk';

for k=1:length(acab)

figure;

for j=1:length(tau)

modelmarks = markers(j);
tauxy0 = sprintf('%i', tau(j));
acab_str = sprintf('%.1f', acab(k));
model_output_filename = strcat(model_prefix, ...
			                   tauxy0,'_kPa_', ...
			                   acab_str,'_acab', ...
                               model_suffix);
                           
[x0,model_yvel,x1,model_ythk,model_thk,uvel,model_vvel] = ...
    nc_read(model_output_filename, -1);

nx1 = length(x1);
nx0 = length(x0);
model_thk = 0.5 * (model_thk(floor(nx1/2),:) + model_thk(floor(nx1/2)+1,:));
model_vvel = 0.5 * (model_vvel(floor(nx0/2),:,1) + model_vvel(floor(nx0/2)+1,:,1));

%[model_ythk,model_thk,model_yvel,model_vvel] = ...
%    read_steadyice(model_output_filename);                           

% load the exact solution from a .mat file
exact_sol_filename = strcat(exact_prefix, ...
                            tauxy0,'_kPa_', ...
                            acab_str, '_acab', ...
                            exact_suffix);
                       
load (exact_sol_filename,  'exact_y', 'exact_thk', 'exact_vvel','exact_thk_pat');

fs = 11;
lwidth = 2;
lwidth_pat = 1;
lwidth_dot = 1;

subplot(2,1,1);
hold on
plot(  exact_y / 1000.0, exact_vvel, 'k','LineWidth',lwidth);  % exact solution
plot(  model_yvel(4:end) / 1000.0, model_vvel(4:end), ...
       strcat(modelmarks,'*'),'LineWidth',lwidth_dot);   % model output

xlabel('y (km)','FontSize',fs);
ylabel('velocity','FontSize',fs);
set(gca,'FontSize',fs);
title(strcat('Steady velocity profile (',...
	     'acab = ', acab_str, ' m/year )'),'FontSize',fs);


subplot(2,1,2);
hold on
plot( exact_y / 1000.0, exact_thk, 'k','LineWidth',lwidth);    %exact solution
plot( exact_y / 1000.0, exact_thk_pat, 'm','LineWidth',lwidth_pat); % Paterson solution
plot( model_ythk(5:end) / 1000.0, model_thk(5:end), ...
      strcat(modelmarks,'*'),'LineWidth',lwidth_dot);   % model output
xlabel('y (km)','FontSize',fs);
ylabel('height','FontSize',fs);
title(strcat('Steady height profile (',...
	     'acab = ', acab_str, ' m/year )'),'FontSize',fs);
set(gca,'FontSize',fs);

end 
subplot(2,1,1);
legend('quasi-analytic (0 kPa)', ...
       'model output (0 kPa)',  ...     
       'quasi-analytic (10 kPa)', ...
       'model output (10 kPa)', ...
       'quasi-analytic (25 kPa)', ...
       'model output (25 kPa)', ... 
       'quasi-analytic (50 kPa)', ...
       'model output (50 kPa)', ...
       'Location','BestOutside') ;

subplot(2,1,2);
legend('quasi-analytic (0 kPa)', ...
       'Paterson sol (0 kPa)', ...
       'model output (0 kPa)',  ...     
       'quasi-analytic (10 kPa)', ...
       'Paterson sol (10 kPa)', ...
       'model output (10 kPa)', ...
       'quasi-analytic (25 kPa)', ...
       'Paterson sol (25 kPa)', ...
       'model output (25 kPa)', ... 
       'quasi-analytic (50 kPa)', ...
       'Paterson sol (50 kPa)', ...
       'model output (50 kPa)', ...
       'Location','BestOutside') ;

end 
