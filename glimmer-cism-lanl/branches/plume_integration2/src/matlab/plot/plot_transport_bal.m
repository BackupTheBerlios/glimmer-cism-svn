
fs = 16;

selectset = @(M) M(1:76,1:155,:);
selectset = @(M) M(2:79,2:80);

%subset = @(M)selectset(four_by_four_avg(M));
subset = @(M)selectset(three_by_three_avg(M));

x = res.x; 
y = res.y; 
[y,x] = meshgrid(y,x);


fig1 = figure(1);
set(fig1,'Position',[1 1 800 600]);
clf;
ncontours = 25;

cmax = 1.4e-3;
fs_label = 16;
label_x = 2.5;
label_y = 17.5;

subplot(2,4,1);
contourf(subset(x),subset(y), ...
         subset(res.su_adv_derv),ncontours,'EdgeColor','None');
caxis([-cmax cmax ]);
title('Flux divergence','FontSize',fs);
set(gca,'FontSize',fs);
text(label_x,label_y,'a','color','k','FontSize',fs_label);

%subplot(2,3,3);
%contourf(subset(x),subset(y), ...
%	 subset(res.u_diff+res.u_diff2),ncontours,'EdgeColor','None');colorbar;
%caxis([-cmax cmax]);
%title('diffusion','FontSize',fs)
%set(gca,'FontSize',fs);
%text(label_x,label_y,'c','color','k','FontSize',fs_label);

subplot(2,4,2);
contourf(subset(x),subset(y), ...
	 subset(res.cor_x),ncontours,'EdgeColor','None');
title('Coriolis','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);
text(label_x,label_y,'b','color','k','FontSize',fs_label);

subplot(2,4,3);
contourf(subset(x),subset(y), ...
	 subset(res.isopyc_grad_x),ncontours,...
	 'EdgeColor','None');
title('Interface flattening','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);
text(label_x,label_y,'d','color','k','FontSize',fs_label);

%subplot(2,3,5);
%contourf(subset(x),subset(y), ...
%	 subset(res.u_boundary), ... %res.den_grad_x),ncontours, ... %bpos-800), ncontours,...
%         ncontours,'EdgeColor','None');colorbar;
%title('density gradient term','FontSize',fs);
%title('boundary term','FontSize',fs);
%title('ice basal surface (m)','FontSize',fs);
%caxis([-cmax cmax]);
%caxis([-600 0]);
%set(gca,'FontSize',fs);
%text(label_x,label_y,'e','color','k','FontSize',fs_label);

%subplot(2,3,6);
%contourf(subset(x),subset(y), ...
%	 subset(res.u_drag), ... % res.u_boundary),
%         ncontours,'EdgeColor','None');colorbar;
%title('boundary term','FontSize',fs);
%title('drag','FontSize',fs);
%caxis([-cmax cmax]);
%set(gca,'FontSize',fs);
%text(label_x,label_y,'f','color','k','FontSize',fs_label);

