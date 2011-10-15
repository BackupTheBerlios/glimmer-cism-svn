
fice = '/scratch/cvg222/gcp/GC_jobs_paper/may29_highres_original/ice.out.nc'
dice = nc_ice_read(fice, 1, -1);

tlen = length(dice.time);
xlen = length(dice.xgrid);
ylen = length(dice.ygrid);

xs = double(dice.thk_t);
%xs = mean(double(dice.efvs),3);

xs = reshape(xs, xlen*ylen, tlen);

xs_nomean = xs - repmat(sum(xs, 2) / tlen, [1 tlen]);

xs_trend_a = zeros(xlen*ylen,1);
xs_trend_b = zeros(xlen*ylen,1);

for i=1:xlen*ylen
  p = polyfit(1:tlen, xs(i,:), 1);
  a = p(1); 
  b = p(2);
  xs(i,:) = xs(i,:) - a*(1:tlen) - b;
  xs_trend_a(i,1) = a;
  xs_trend_b(i,1) = b;
end

n_eigs = 10

corr_mat = xs_nomean * xs_nomean' / tlen;
total_variance = trace(corr_mat)

[u_eofs,vars] = eigs(corr_mat, n_eigs);
[v_eofs,vars] = eigs(xs_nomean'* xs_nomean  / tlen, n_eigs);

scores = v_eofs * vars;
u_eofs = reshape(u_eofs,xlen,ylen,n_eigs);
xs_trend_a = reshape(xs_trend_a,xlen,ylen);
xs_trend_b = reshape(xs_trend_b,xlen,ylen);


fs = 16

subplot(2,2,1);
contourf(dice.xgrid/1000,dice.ygrid/1000,squeeze(u_eofs(:,:,1))');colorbar;
title('first EOF of variations in thk\_t','fontsize',fs)
set(gca,'fontsize',fs);

subplot(2,2,2);
contourf(dice.xgrid/1000,dice.ygrid/1000,squeeze(u_eofs(:,:,2))');colorbar;
title('second EOF of variations in thk\_t','fontsize',fs)
set(gca,'fontsize',fs);

subplot(2,2,3);
plot((1:tlen)/20.0, v_eofs(:,1),'b');
title(strcat(['time series of first EOF (', sprintf('%03.2f',vars(1,1)/total_variance), ...
					  ' fraction of variance)']),'fontsize',fs);
set(gca,'fontsize',fs);

subplot(2,2,4);                     
plot((1:tlen)/20.0, v_eofs(:,2),'b');
title(strcat(['time series of second EOF (', sprintf('%03.2f',vars(2,2)/total_variance), ...
					  ' fraction of variance)']),'fontsize',fs);
set(gca,'fontsize',fs);

