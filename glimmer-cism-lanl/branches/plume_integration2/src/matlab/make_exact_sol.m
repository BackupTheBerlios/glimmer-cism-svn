

file_prefix = '/home/gladish/research/gcp_resources/mat_files/ssj1_';
file_suffix = '_exact.sept14.400m.mat';

m = 50;
n = 100;

hx = 20000.0 / m;
hy = 40000.0 / n;

upvel = -1000.0;
upthk = 600.0;
yin = (n - 2.5) * hy;  % correct for velocity
%yin = (n - 2.0) * hy;   % correct for thickness
yout = 4.0 * hy;  % correct for thickness

tauxy0 = [0,10, 25, 50];
acab = [-2.0, 0.0, 2.0];

L = (m-1)/2.0 * hx;  % half-width between outermost velocity points


for i=1:length(tauxy0)
 for j=1:length(acab)

[exact_y, exact_thk, exact_vvel, exact_thk_pat] = ...
    quasi_exact_shear_sol(upvel, upthk, yin, ...
                          yout, acab(j), tauxy0(i)*1000.0, ...
                          L, ny);

filename = strcat(file_prefix, ...
		  sprintf('%i',tauxy0(i)),'_kPa_', ...
		  sprintf('%2.1f',acab(j)),'_acab', ...
		  file_suffix);

save(filename, 'exact_y', 'exact_thk', 'exact_vvel', 'exact_thk_pat');

  end
end



