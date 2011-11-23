%trim_edge = @(M) M(4:end-3,4:end-3,:);

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

cmax = 4.0e-5;
fs_label = 16;
label_x = 2.5;
label_y = 17.5;

subplot(2,3,1);
contourf(subset(x),subset(y), ...
         subset(res.su_adv_derv),ncontours,'EdgeColor','None');colorbar;
caxis([-cmax cmax ]);
title('inertial term','FontSize',fs);
set(gca,'FontSize',fs);
text(label_x,label_y,'a','color','k','FontSize',fs_label);

subplot(2,3,3);
contourf(subset(x),subset(y), ...
	 subset(res.u_diff+res.u_diff2),ncontours,'EdgeColor','None');colorbar;
caxis([-cmax cmax]);
title('diffusion','FontSize',fs)
set(gca,'FontSize',fs);
text(label_x,label_y,'c','color','k','FontSize',fs_label);

subplot(2,3,2);
contourf(subset(x),subset(y), ...
	 subset(res.cor_x),ncontours,'EdgeColor','None');colorbar;
title('coriolis','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);
text(label_x,label_y,'b','color','k','FontSize',fs_label);

subplot(2,3,4);
contourf(subset(x),subset(y), ...
	 subset(res.isopyc_grad_x),ncontours,...
	 'EdgeColor','None');colorbar;
title('interface flattening','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);
text(label_x,label_y,'d','color','k','FontSize',fs_label);

subplot(2,3,5);
contourf(subset(x),subset(y), ...
	 subset(res.u_boundary), ... %res.den_grad_x),ncontours, ... %bpos-800), ncontours,...
         ncontours,'EdgeColor','None');colorbar;
title('density gradient term','FontSize',fs);
title('boundary term','FontSize',fs);
%title('ice basal surface (m)','FontSize',fs);
caxis([-cmax cmax]);
%caxis([-600 0]);
set(gca,'FontSize',fs);
text(label_x,label_y,'e','color','k','FontSize',fs_label);

subplot(2,3,6);
contourf(subset(x),subset(y), ...
	 subset(res.u_drag), ... % res.u_boundary),
         ncontours,'EdgeColor','None');colorbar;
title('boundary term','FontSize',fs);
title('drag','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);
text(label_x,label_y,'f','color','k','FontSize',fs_label);

%%%%%%%% same thing for y direction
fig2 = figure(2);
set(fig2,'Position',[1 1 800 600]);
clf;
ncontours = 20;
cmax = 1.0e-4;
cmax = 4.0e-5;

subplot(2,3,1);
contourf(subset(x),subset(y), ...
         subset(res.sv_adv_derv),ncontours,'EdgeColor','None');colorbar;
caxis([-cmax cmax ]);
title('inertial term','FontSize',fs);
set(gca,'FontSize',fs);
text(label_x,label_y,'a','color','k','FontSize',fs_label);

subplot(2,3,2);
contourf(subset(x),subset(y), ...
	 subset(res.cor_y),ncontours,'EdgeColor','None');colorbar;
title('coriolis','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);
text(label_x,label_y,'b','color','k','FontSize',fs_label);

subplot(2,3,3);
contourf(subset(x),subset(y), ...
	 subset(res.v_diff+res.v_diff2),ncontours,'EdgeColor','None');colorbar;
caxis([-cmax cmax]);
title('diffusion','FontSize',fs)
set(gca,'FontSize',fs);
text(label_x,label_y,'c','color','k','FontSize',fs_label);

subplot(2,3,4);
contourf(subset(x),subset(y), ...
	 subset(res.isopyc_grad_y),ncontours,...
	 'EdgeColor','None');colorbar;
title('interface flattening','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);
text(label_x,label_y,'d','color','k','FontSize',fs_label);

subplot(2,3,5);
contourf(subset(x),subset(y), ...
         subset(res.v_boundary),ncontours, ... 
	 'EdgeColor','None');colorbar;
title('boundary term','FontSize',fs);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);
text(label_x,label_y,'e','color','k','FontSize',fs_label);

subplot(2,3,6);
contourf(subset(x),subset(y), ...
	 subset(res.v_drag), ... %_den_grad_y),ncontours,... %res.bpos-800.0),ncontours, ...
         ncontours,'EdgeColor','None');colorbar;
title('density gradient term','FontSize',fs);
title('drag','FontSize',fs);
%title('ice basal surface (m)','FontSize',fs);
%caxis([-600 0]);
caxis([-cmax cmax]);
set(gca,'FontSize',fs);
text(label_x,label_y,'f','color','k','FontSize',fs_label);

%%%%%%%%%%%%%%%%
%  More scatter plots showing balances in x and y direction
%%%%%%%%%%%%%%%

figure(4);
clf;

subplot(1,2,1);
hold on
lhs = 1.0*res.su_adv_derv+1.0*res.cor_x;
rhs =  1.0*res.den_grad_x + 1.0*res.u_diff-1.0*res.u_diff2 + ...
       1.0*res.u_boundary + 1.0*res.u_drag + 1.0*res.isopyc_grad_x;

a = flatten_field(subset(rhs));
b= flatten_field(subset(lhs));

limlo = -0.5e-4;
limhi = 0.5e-4;

xlim([limlo limhi]);
ylim([limlo limhi]);
plot(a,b,'b.');
plot([limlo, limhi],[limlo, limhi],'r-');
m = corr([a,b]);
%title('x balance','FontSize',fs);
xlabel('RHS terms (m/s^2)','FontSize',fs);
ylabel('LHS acceleration terms (m/s^2)','FontSize',fs);
set(gca,'FontSize',fs);
title(strcat(['correlation =',sprintf('%3.2f',m(1,2))]),'FontSize',fs);

subplot(1,2,2);
hold on
lhs = 1.0*res.sv_adv_derv+1.0*res.cor_y;
rhs = 1.0*res.den_grad_y + 1.0*res.v_diff-1.0*res.v_diff2 + ...
      1.0*res.v_boundary + 1.0*res.v_drag + 1.0*res.isopyc_grad_y;
a = flatten_field(subset(rhs));
b= flatten_field(subset(lhs));


m = corr([a b]);
plot(a,b,'b.');
limlo = -0.5e-4;
limhi = 0.5e-4;
plot([limlo, limhi],[limlo, limhi],'r-');
xlim([limlo limhi]);
ylim([limlo limhi]);

set(gca,'FontSize',fs);
xlabel('RHS terms (m/s^2)','FontSize',fs);
ylabel('LHS acceleration terms (m/s^2)','FontSize',fs);
title(strcat(['correlation =',sprintf('%3.2f',m(1,2))]),'FontSize',fs);

if (false)
figure(5);
clf;
contourf(subset(x),subset(y),subset(res.su_adv_derv+res.cor_x - res.u_diff-res.u_diff2-res.den_grad_x-res.isopyc_grad_x-res.u_boundary),'EdgeColor','None');colorbar;

figure(6);
clf;
contourf(subset(x),subset(y),subset(res.sv_adv_derv+res.cor_y - res.v_diff-res.v_diff2-res.den_grad_y-res.isopyc_grad_y-res.v_boundary),'EdgeColor','None');colorbar;
end 
