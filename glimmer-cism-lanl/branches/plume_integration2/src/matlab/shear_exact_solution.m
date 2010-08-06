% read in netcdf model output

jtype = ['1','2'];
jtype =['2'];
tau = [0,10,25,50];
tau = [0,10,25,50];
acab = [-2.0, 0.0, 2.0];

markers = 'brgk';
for i=1:length(jtype)
for k=1:length(acab)

figure;

for j=1:length(tau)

modelmarks = markers(j);
tauxy0 = sprintf('%i', tau(j));
acab_str = sprintf('%.1f', acab(k));
prefix = strcat('ssj', jtype(i),'_');
model_output_filename = strcat(prefix, ...
			       tauxy0,'_kPa_', ...
			       acab_str,'_acab_model.mat');
load(model_output_filename, 'model_ythk', 'model_thk', 'model_vvel', 'model_yvel');

% load the exact solution from a .mat file
exact_sol_filename = strcat(prefix, ...
	                    tauxy0,'_kpa_', ...
			    acab_str, '_acab_exact.mat');
load (exact_sol_filename,  'exact_y', 'exact_thk', 'exact_vvel',  'exact_thk_pat');

fs = 20;


subplot(2,1,1);
hold on
plot(  exact_y / 1000.0, exact_vvel, 'k');  % exact solution
plot(  model_yvel(4:end) / 1000.0, model_vvel(4:end), ...
       strcat(modelmarks,'*'));   % model output

xlabel('y (km)','FontSize',fs);
ylabel('velocity','FontSize',fs);
set(gca,'FontSize',fs);
title(strcat('Steady velocity profile (',...
	     'acab = ', acab_str, ' m/year )'),'FontSize',fs);


subplot(2,1,2);
hold on
plot( exact_y / 1000.0, exact_thk, 'k');    %exact solution
plot( exact_y / 1000.0, exact_thk_pat, 'm'); % Paterson solution
plot( model_ythk(5:end) / 1000.0, model_thk(5:end), ...
      strcat(modelmarks,'*'));   % model output
xlabel('y (km)','FontSize',fs);
ylabel('height','FontSize',fs);
title(strcat('Steady height profile (',...
	     'acab = ', acab_str, ' m/year )'),'FontSize',fs);
set(gca,'FontSize',fs);



end 
subplot(2,1,1);
legend('quasi-analytic', ...
       'model output (0 kPa)',       
       '','model output (10 kPa)',
       '','model output (25 kPa)', ... 
       '','model output (50 kPa)', ...
       'Location','Outside') ;

subplot(2,1,2);
legend('quasi-analytic', ...
       'Paterson solution', ...
       'model output (0 kPa)', ...
       '','', 'model output (10 kPa)', ...
       '','', 'model output (25 kPa)', ...
       '','', 'model output (50 kPa)', ...
       'Location','Outside') ;

end 
end
