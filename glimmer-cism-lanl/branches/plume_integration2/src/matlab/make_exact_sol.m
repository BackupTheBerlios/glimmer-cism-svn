
upvel = -1000.0;
upthk = 600.0;
yin = 37.5 * 1000.0;  % correct for velocity
yin = 38.0 * 1000.0;   % correct for thickness
yout = 4.0 * 1000.0;  % correct for thickness
tauxy0 = [0,10, 25, 50];
acab = [-2.0, 0.0, 2.0];
jtype = ['1','2'];

L = 19.0/2.0 * 1000.0;  % half-width between outermost velocity points
ny = 40;

for k=1:length(jtype)
  for i=1:length(tauxy0)
    for j=1:length(acab)

[exact_y, exact_thk, exact_vvel, exact_thk_pat] = ...
    quasi_exact_shear_sol(upvel, upthk, yin, ...
                          yout, acab(j), tauxy0(i)*1000.0, ...
                          L, ny);

filename = strcat('ssj',jtype(k),'_', ...
		  sprintf('%i',tauxy0(i)),'_kPa_', ...
		  sprintf('%2.1f',acab(j)),'_acab_exact.mat');
save(filename, 'exact_y', 'exact_thk', 'exact_vvel', 'exact_thk_pat');

end 
end
end



