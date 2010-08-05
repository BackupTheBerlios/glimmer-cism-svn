
upvel = -1000.0;
upthk = 600.0;
yin = 40.0 * 1000.0;
yout = 5.0 * 1000.0;
acab = 0.0;
tauxy0 = 50.0 * 1000.0;
L = 10.0 * 1000.0;
ny = 40;


[exact_y, exact_thk, exact_vvel, exact_thk_pat] = ...
    quasi_exact_shear_sol(upvel, upthk, yin, ...
                          yout, acab, tauxy0, ...
                          L, ny);


save 'ssj2_50_kPa_0.0_acab_exact.mat' exact_y exact_thk exact_vvel ...
    exact_thk_pat;



