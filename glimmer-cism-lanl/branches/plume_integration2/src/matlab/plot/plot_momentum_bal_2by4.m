trim_edge = @(M) M(4:end-3,4:end-3,:);
fs = 16;
fs_label = 16;
label_x = 2.5;
label_y = 17.5;

selectset = @(M) M(2:79,2:80);

%subset = @(M)selectset(four_by_four_avg(M));
%subset = @(M)selectset(two_by_two_avg(M));
subset = @(M)selectset(three_by_three_avg(M));

x = res.x; 
y = res.y; 
[y,x] = meshgrid(y,x);

fig1 = figure(1);
set(fig1,'Position',[1 1 800 600]);
clf;
ncontours = 20;
cmax = 1.45e-3;

subplot(2,4,1);
contourf(subset(x),subset(y), ...
         subset(res.u_flux_div),ncontours,'EdgeColor','None');
caxis([-cmax cmax ]);
title('Inertial','FontSize',fs);
set(gca,'FontSize',fs);
text(label_x,label_y,'a','color','k','FontSize',fs_label);

%subplot(2,3,3);
%contourf(subset(x),subset(y), ...
%	 subset(res.u_diff),ncontours,'EdgeColor','None');colorbar;
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
title('Interface','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);
text(label_x,label_y,'c','color','k','FontSize',fs_label);

%subplot(2,3,5);
%contourf(subset(x),subset(y), ...
%	 subset(res.den_grad_x), ncontours,...
%	 'EdgeColor','None');colorbar;
%contourf(subset(x),subset(y), ...
%	 subset(sqrt(res.u_adv_derv.^2+res.v_adv_derv.^2)), ...
%	 'EdgeColor','None');colorbar;
%title('density gradient','FontSize',fs);
%caxis([-cmax cmax]);
%set(gca,'FontSize',fs);
%text(label_x,label_y,'e','color','k','FontSize',fs_label);

subplot(2,4,4);
contourf(subset(x),subset(y), ...
	 subset(res.u_trans_detrain),ncontours,'EdgeColor','None');colorbar;
title('Detrainment','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);
text(label_x,label_y,'d','color','k','FontSize',fs_label);

%%%%%%%% same thing for y direction
cmax = 1.25e-3;

subplot(2,4,5);
contourf(subset(x),subset(y), ...
         subset(res.v_flux_div),ncontours,'EdgeColor','None');
caxis([-cmax cmax ]);
title('Inertial','FontSize',fs);
set(gca,'FontSize',fs);
text(label_x,label_y,'e','color','k','FontSize',fs_label);

subplot(2,4,6);
contourf(subset(x),subset(y), ...
	 subset(res.cor_y),ncontours,'EdgeColor','None');
title('Coriolis','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);
text(label_x,label_y,'f','color','k','FontSize',fs_label);

subplot(2,4,7);
contourf(subset(x),subset(y), ...
	 subset(res.isopyc_grad_y),ncontours,...
	 'EdgeColor','None');
title('Interface','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);
text(label_x,label_y,'g','color','k','FontSize',fs_label);

subplot(2,4,8);
contourf(subset(x),subset(y), ...
	      subset(res.drag_y),ncontours,'EdgeColor','None');colorbar;
contourf(subset(x),subset(y), ...
	 subset(res.v_trans_detrain+0.0*res.isopyc_grad_y+ ...
             0.0*res.v_diff-0.0*res.cor_y+0.0*res.drag_y), ...
	 ncontours,'EdgeColor','None');colorbar;
title('Detrainment','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);
text(label_x,label_y,'h','color','k','FontSize',fs_label);
