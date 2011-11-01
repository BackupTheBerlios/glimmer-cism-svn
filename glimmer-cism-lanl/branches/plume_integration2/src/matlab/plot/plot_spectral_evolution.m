function [spec,norm2] = plot_spectral_evolution(dice)

[m,n,k] = size(dice.thk);

meanthk = repmat(mean(dice.thk,1),[m,1,1]);
dev_thk = dice.thk-meanthk;

norm2 = zeros(1,n,k);

for j=1:n
  norm2(1,j,:) = sum(squeeze(dev_thk(:,j,:))' * squeeze(dev_thk(:,j,:)));
end


spec = fft(dev_thk,[],1) / sqrt(size(dev_thk,1));

kmax = 30;

ks = 0:(length(dice.xgrid)-1);

times = [1,k];
times = [1];

fs = 12;

for t=1:length(times)
%  subplot(2,length(times),t);
subplot(1,2,1);
  contourf(dice.xgrid/1000,dice.ygrid/1000,dice.thk(:,:,times(t))',20,'EdgeColor','None');colorbar;
xlabel('Across shelf distance (km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);
title('Ice thickness (m)','FontSize',fs);

%  subplot(2,length(times),t+length(times));
subplot(1,2,2);

contourf(ks(1:kmax),dice.ygrid/1000,log(abs(spec(1:kmax,:,times(t))'.* ...
                    		       conj(spec(1:kmax,:,times(t)))')),30,'EdgeColor','None');
  colorbar;

caxis([-5 10]);

if (t==1)

else
%caxis([0 6.0e4]);
end

xlabel('Across shelf wavenumber (1/km)','FontSize',fs);
ylabel('Along shelf distance (km)','FontSize',fs);
title('log spectral density of cross-shelf variations (m^2)','FontSize',fs);


end

fs2 = 16;
%subplot(2,2,1);text(-5,45,'initial state: steady ice, imposed melting','FontSize',fs2);
%subplot(2,2,2);text(-5,45,'final state: steady ice and ocean, dynamic melting','FontSize',fs2);

end

