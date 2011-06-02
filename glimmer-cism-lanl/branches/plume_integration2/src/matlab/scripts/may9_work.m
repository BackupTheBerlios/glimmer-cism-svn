
[ug,vg,uag,vag,agfrac] = geostrophic_vel(d_plume_trial_90);

figure(1);
clf;
fs = 16;

subplot(2,2,1);
contourf(d_plume_trial_90.x/1000.0,d_plume_trial_90.y/1000.0, ...
         ug(:,:,800)',25,'EdgeColor','None');colorbar;
xlabel('Across shelf distance','FontSize',fs);
ylabel('Along shelf distance','FontSize',fs);
title('u geostrophic','FontSize',fs);

subplot(2,2,2);
contourf(d_plume_trial_90.x/1000.0,d_plume_trial_90.y/1000.0, ...
         vg(:,:,800)',25,'EdgeColor','None');colorbar;
xlabel('Across shelf distance','FontSize',fs);
ylabel('Along shelf distance','FontSize',fs);
title('v geostrophic','FontSize',fs);

subplot(2,2,3);
contourf(d_plume_trial_90.x/1000.0,d_plume_trial_90.y/1000.0, ...
         uag(:,:,800)',25,'EdgeColor','None');colorbar;
xlabel('Across shelf distance','FontSize',fs);
ylabel('Along shelf distance','FontSize',fs);
title('u ageostrophic','FontSize',fs);

subplot(2,2,4);
contourf(d_plume_trial_90.x/1000.0,d_plume_trial_90.y/1000.0, ...
         vag(:,:,800)',25,'EdgeColor','None');colorbar;
     
xlabel('Across shelf distance','FontSize',fs);
ylabel('Along shelf distance','FontSize',fs);
title('v ageostrophic','FontSize',fs);


figure(2);
clf;

%subplot(1,2,1);
hold on
%title('ageostrophic vel','FontSize',fs);

%subplot(1,2,2);
%contourf(d_plume_trial_90.x/1000.0,d_plume_trial_90.y/1000.0, ...
%         d_plume_trial_90.bpos(:,:,800)',25,'EdgeColor','None');colorbar;
%scale = 2.0
%quiver(d_plume_trial_90.x/1000.0,d_plume_trial_90.y/1000.0, ...
%       uag(:,:,800)',vag(:,:,800)',scale,'w');  
%quiver(d_plume_trial_90.x/1000.0,d_plume_trial_90.y/1000.0, ...
%       ug(:,:,800)',vg(:,:,800)',scale,'k');   

su = d_plume_trial_90.su;
sv = d_plume_trial_90.sv;

innerprod = ((su .* ug) + (sv .* vg))./(sqrt(ug.*ug+ vg.*vg).*sqrt(su.*su+sv.*sv));
proj = ((su .* ug) + (sv .* vg))./(sqrt(ug.*ug+ vg.*vg).*sqrt(su.*su+sv.*sv));

%contourf(d_plume_trial_90.x/1000.0,d_plume_trial_90.y/1000.0, ...
%         180/pi*acos(innerprod(:,:,800))',25,'EdgeColor','None'); colorbar;
subplot(2,2,2);
contourf(d_plume_trial_90.x/1000.0,d_plume_trial_90.y/1000.0, ...
         sqrt(uag(:,:,800).^2+vag(:,:,800).^2)',25,'EdgeColor','None'); colorbar;
      
xlabel('Across shelf distance','FontSize',fs);
ylabel('Along shelf distance','FontSize',fs);
title('Ageostrophic speed','FontSize',fs);

subplot(2,2,1);
contourf(d_plume_trial_90.x/1000.0,d_plume_trial_90.y/1000.0, ...
         sqrt(ug(:,:,800).^2+vg(:,:,800).^2)',25,'EdgeColor','None'); colorbar;
      
xlabel('Across shelf distance','FontSize',fs);
ylabel('Along shelf distance','FontSize',fs);
title('Geostrophic speed','FontSize',fs);

subplot(2,2,3);
angle =  acos((su.*ug+sv.*vg)./(sqrt(ug.*ug+vg.*vg).*sqrt(su.*su+sv.*sv)));

contourf(d_plume_trial_90.x/1000.0,d_plume_trial_90.y/1000.0, ...
        angle(:,:,800)',25,'EdgeColor','None'); colorbar;

xlabel('Across shelf distance','FontSize',fs);
ylabel('Along shelf distance','FontSize',fs);
title('Angle between u and u_g','FontSize',fs);


subplot(2,2,4);
angle =  acos((su.*ug+sv.*vg)./(sqrt(ug.*ug+vg.*vg).*sqrt(su.*su+sv.*sv)));

contourf(d_plume_trial_90.x/1000.0,d_plume_trial_90.y/1000.0, ...
        d_plume_trial_90.bpos(:,:,800)',25,'EdgeColor','None'); colorbar;

xlabel('Across shelf distance','FontSize',fs);
ylabel('Along shelf distance','FontSize',fs);
title('Angle between u and u_g','FontSize',fs);


%xlabel('Across shelf distance','FontSize',fs);
%ylabel('Along shelf distance','FontSize',fs);
%title('geostrophic vel','FontSize',fs);




figure(3);
clf;
[m,n,k] = size(su);
flatten = @(M) reshape(squeeze(M(2:(end-1),2:(end-1),:)), ...
                              (size(M,1)-2)*(size(M,2)-2)*size(M,3),1);

subplot(2,2,1);
hold on
plot(flatten(ug)',flatten(su)','b.');
plot(flatten(ug(:,:,800))',flatten(su(:,:,800))','r.');
plot([min(min(su(:,:,800))):0.001:max(max(su(:,:,800)))], ...
     [min(min(su(:,:,800))):0.001:max(max(su(:,:,800)))],'k.');
ylabel('u','FontSize',fs);
xlabel('u geo','FontSize',fs);
set(gca,'FontSize',fs);
hold off

subplot(2,2,2);
hold on
plot(flatten(vg)',flatten(sv)','b.');
plot(flatten(vg(:,:,800))',flatten(sv(:,:,800))','r.');
plot([min(min(sv(:,:,800))):0.001:max(max(sv(:,:,800)))], ...
     [min(min(sv(:,:,800))):0.001:max(max(sv(:,:,800)))],'k.');
ylabel('v','FontSize',fs);
xlabel('v geo','FontSize',fs);
set(gca,'FontSize',fs);
hold off

subplot(2,2,3);
hold on
plot(flatten(uag)',flatten(su)','b.');
plot(flatten(uag(:,:,800))',flatten(su(:,:,800))','r.');
plot([min(min(su(:,:,800))):0.001:max(max(su(:,:,800)))], ...
     [min(min(su(:,:,800))):0.001:max(max(su(:,:,800)))],'k.');
ylabel('u','FontSize',fs);
set(gca,'FontSize',fs);
xlabel('u ageo','FontSize',fs);
hold off

subplot(2,2,4);
hold on
plot(flatten(vag)',flatten(sv)','b.');
plot(flatten(vag(:,:,800))',flatten(sv(:,:,800))','r.');
plot([min(min(sv(:,:,800))):0.001:max(max(sv(:,:,800)))], ...
     [min(min(sv(:,:,800))):0.001:max(max(sv(:,:,800)))],'k.');
ylabel('v','FontSize',fs);
xlabel('v ageo ','FontSize',fs);
set(gca,'FontSize',fs);
hold off



