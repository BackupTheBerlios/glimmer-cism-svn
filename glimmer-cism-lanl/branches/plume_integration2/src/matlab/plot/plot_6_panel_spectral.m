fig1 = figure(1);
clf;

jobs = {'k_1_amp_25',...
	'k_5_amp_25',...
	'k_6_amp_25',...
	'k_1_amp_50',...
	'k_5_amp_50',...
	'k_6_amp_50'};

for i=1:length(jobs)

  d = dice.(jobs{i});

[m,n,k] = size(d.thk);
fs = 16;

meanthk = repmat(mean(d.thk,1),[m,1,1]);
dev_thk = d.thk-meanthk;

norm2 = zeros(1,n,k);

for j=1:n
  norm2(1,j,:) = sum(squeeze(dev_thk(:,j,:))' * squeeze(dev_thk(:,j,:)));
end

spec = fft(dev_thk,[],1) / sqrt(size(dev_thk,1));

kmax = 20;

ks = 0:(length(d.xgrid)-1);

times = [1,k];
times = [1];


for t=1:length(times)

subplot(2,3,i);
set(fig1,'Position',[1 1 800 600]);

z = spec(1:kmax,:,times(t)).* conj(spec(1:kmax,:,times(t)));
z = log(max(z,1.0));
contourf(ks(1:kmax),d.ygrid/1000,z',40,'EdgeColor','None');
%colorbar;
set(gca,'FontSize',fs);
caxis([0 10.5]);

textks = {'1','5','6'};
if (i<4)
  text(2,45,sprintf('k=%s',textks{i}),'FontSize',24);
end

if (i==1)
   text(-10,5,'amp = 25 m','FontSize',24,'rotation',90);
%set(h,'rotation',90);
end
if (i==4)
  text(-10,5,'amp = 50 m','FontSize',24,'rotation',90);
end
xlabel('k ','FontSize',fs);
ylabel('y (km)','FontSize',fs);
set(gca,'xtick',0:5:20)
  if(mod(i,3)==0)
   h = colorbar('East','FontSize',fs);
end
%title(sprintf('log(spectral density) of cross-shelf variations: k=%s amp=%s',..
 % k_pert,amp),'FontSize',fs);

end
end

%set(gcf,'PaperPositionMode','auto');
%print('-depsc',figfile)

%end

