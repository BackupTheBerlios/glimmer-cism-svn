function [spec,norm2] = plot_spectral_evolution(dice,fig_dir)

[m,n,k] = size(dice.thk);
fs = 16;

meanthk = repmat(mean(dice.thk,1),[m,1,1]);
dev_thk = dice.thk-meanthk;

norm2 = zeros(1,n,k);

for j=1:n
  norm2(1,j,:) = sum(squeeze(dev_thk(:,j,:))' * squeeze(dev_thk(:,j,:)));
end


spec = fft(dev_thk,[],1) / sqrt(size(dev_thk,1));

kmax = 20;

ks = 0:(length(dice.xgrid)-1);

times = [1,k];
times = [1];


for t=1:length(times)


%subplot(1,2,1);
%contourf(dice.xgrid/1000,dice.ygrid/1000,dice.thk(:,:,times(t))',20,'EdgeColor','None');
%colorbar;
%xlabel('Across shelf distance (km)','FontSize',fs);
%ylabel('Along shelf distance (km)','FontSize',fs);
%title('Ice thickness (m)','FontSize',fs);
%set(gca,'FontSize',fs);

%  subplot(2,length(times),t+length(times));

%subplot(1,2,2);
fig1 = figure(1);
clf;
set(fig1,'Position',[1 1 800 600]);

z = spec(1:kmax,:,times(t)).* conj(spec(1:kmax,:,times(t)));
z = log(max(z,1.0));
contourf(ks(1:kmax),dice.ygrid/1000,z',40,'EdgeColor','None');
colorbar;
set(gca,'FontSize',fs);
%caxis([-5 10]);

xlabel('Across shelf wavenumber (1/km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);
title('log(spectral density) for cross-shelf variations','FontSize',fs);
%title('spectral density of cross-shelf variations (m^2)','FontSize',fs);

end

print('-depsc',strcat([fig_dir,'/plume_spectral_plot']));

end

