

% read in netcdf model output
ncfilename = '/Users/carl/computation/jobs/job_data/ssj2_50_kPa_0.0_acab/ssj2_50_kPa_0.0_acab.out.nc'
[y_thk, thk, y_vel, vel, thk_paterson] = read_steadyice(ncfilename)

% load the exact solution from a .mat file
exact_sole_filename = 'ssj2_50_kpa_0.0_acab_exact.m'
load exact_sol_filename exact_y exact_thk exact_vvel exact_thk_paterson

fs = 20;

figure;


subplot(2,1,1);
hold on
plot(  exact_y / 1000.0, exact_vvel, 'r*');  % exact solution
plot(  y_vel / 1000.0, vel, 'k*');   % model output

xlabel('y (km)','FontSize',fs);
ylabel('velocity','FontSize',fs);
set(gca,'FontSize',fs);
title('Centerline velocities (tauxy0 = 50 kPa)','FontSize',fs);
legend('quasi-analytic', ...
       'model output','Location','BestOutside') ;

hold off

subplot(2,1,2);
hold on
plot( exact_y / 1000.0, exact_thk, 'b*');    %exact solution
plot( y_thk / 1000.0, thk, 'k*');   % model output
plot( exact_y / 1000.0, thk_paterson, 'g*'); % Paterson solution
xlabel('y (km)','FontSize',fs);
ylabel('height','FontSize',fs);
title('Centerline thickness (tauxy0 = 50 kPa)','FontSize',fs);
set(gca,'FontSize',fs);
legend('quasi-analytic','model output','Paterson','Location','BestOutside') ;
hold off

