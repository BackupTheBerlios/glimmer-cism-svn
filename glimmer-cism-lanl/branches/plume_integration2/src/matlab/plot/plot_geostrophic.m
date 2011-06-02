function plot_geostrophic(dplume, tslice)

[ug,vg,uag,vag,agfrac] = geostrophic_vel(dplume);

figure(1);
clf;
fs = 16;

subplot(2,2,1);
contourf(dplume.x/1000.0,dplume.y/1000.0, ...
         ug(:,:,tslice)',25,'EdgeColor','None');colorbar;
xlabel('Across shelf distance','FontSize',fs);
ylabel('Along shelf distance','FontSize',fs);
title('u geostrophic','FontSize',fs);

subplot(2,2,2);
contourf(dplume.x/1000.0,dplume.y/1000.0, ...
         vg(:,:,tslice)',25,'EdgeColor','None');colorbar;
xlabel('Across shelf distance','FontSize',fs);
ylabel('Along shelf distance','FontSize',fs);
title('v geostrophic','FontSize',fs);

subplot(2,2,3);
contourf(dplume.x/1000.0,dplume.y/1000.0, ...
         uag(:,:,tslice)',25,'EdgeColor','None');colorbar;
xlabel('Across shelf distance','FontSize',fs);
ylabel('Along shelf distance','FontSize',fs);
title('u ageostrophic','FontSize',fs);

subplot(2,2,4);
contourf(dplume.x/1000.0,dplume.y/1000.0, ...
         vag(:,:,tslice)',25,'EdgeColor','None');colorbar;
     
xlabel('Across shelf distance','FontSize',fs);
ylabel('Along shelf distance','FontSize',fs);
title('v ageostrophic','FontSize',fs);


figure(2);
clf;

%subplot(1,2,1);
hold on
%title('ageostrophic vel','FontSize',fs);

%subplot(1,2,2);
%contourf(dplume.x/1000.0,dplume.y/1000.0, ...
%         dplume.bpos(:,:,tslice)',25,'EdgeColor','None');colorbar;
%scale = 2.0
%quiver(dplume.x/1000.0,dplume.y/1000.0, ...
%       uag(:,:,tslice)',vag(:,:,tslice)',scale,'w');  
%quiver(dplume.x/1000.0,dplume.y/1000.0, ...
%       ug(:,:,tslice)',vg(:,:,tslice)',scale,'k');   

su = dplume.su;
sv = dplume.sv;

innerprod = ((su .* ug) + (sv .* vg))./(sqrt(ug.*ug+ vg.*vg).*sqrt(su.*su+sv.*sv));
proj = ((su .* ug) + (sv .* vg))./(sqrt(ug.*ug+ vg.*vg).*sqrt(su.*su+sv.*sv));

%contourf(dplume.x/1000.0,dplume.y/1000.0, ...
%         180/pi*acos(innerprod(:,:,tslice))',25,'EdgeColor','None'); colorbar;
subplot(2,2,2);
contourf(dplume.x/1000.0,dplume.y/1000.0, ...
         sqrt(uag(:,:,tslice).^2+vag(:,:,tslice).^2)',25,'EdgeColor','None'); colorbar;
      
xlabel('Across shelf distance','FontSize',fs);
ylabel('Along shelf distance','FontSize',fs);
title('Ageostrophic speed','FontSize',fs);

subplot(2,2,1);
contourf(dplume.x/1000.0,dplume.y/1000.0, ...
         sqrt(ug(:,:,tslice).^2+vg(:,:,tslice).^2)',25,'EdgeColor','None'); colorbar;
      
xlabel('Across shelf distance','FontSize',fs);
ylabel('Along shelf distance','FontSize',fs);
title('Geostrophic speed','FontSize',fs);

subplot(2,2,3);
angle =  acos((su.*ug+sv.*vg)./(sqrt(ug.*ug+vg.*vg).*sqrt(su.*su+sv.*sv)));

contourf(dplume.x/1000.0,dplume.y/1000.0, ...
        angle(:,:,tslice)',25,'EdgeColor','None'); colorbar;

xlabel('Across shelf distance','FontSize',fs);
ylabel('Along shelf distance','FontSize',fs);
title('Angle between u and u_g','FontSize',fs);


subplot(2,2,4);
angle =  acos((su.*ug+sv.*vg)./(sqrt(ug.*ug+vg.*vg).*sqrt(su.*su+sv.*sv)));

contourf(dplume.x/1000.0,dplume.y/1000.0, ...
        dplume.bpos(:,:,tslice)',25,'EdgeColor','None'); colorbar;

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
plot(flatten(ug(:,:,tslice))',flatten(su(:,:,tslice))','r.');
plot([min(min(su(:,:,tslice))):0.001:max(max(su(:,:,tslice)))], ...
     [min(min(su(:,:,tslice))):0.001:max(max(su(:,:,tslice)))],'k.');
ylabel('u','FontSize',fs);
xlabel('u geo','FontSize',fs);
set(gca,'FontSize',fs);
hold off

subplot(2,2,2);
hold on
plot(flatten(vg)',flatten(sv)','b.');
plot(flatten(vg(:,:,tslice))',flatten(sv(:,:,tslice))','r.');
plot([min(min(sv(:,:,tslice))):0.001:max(max(sv(:,:,tslice)))], ...
     [min(min(sv(:,:,tslice))):0.001:max(max(sv(:,:,tslice)))],'k.');
ylabel('v','FontSize',fs);
xlabel('v geo','FontSize',fs);
set(gca,'FontSize',fs);
hold off

subplot(2,2,3);
hold on
plot(flatten(uag)',flatten(su)','b.');
plot(flatten(uag(:,:,tslice))',flatten(su(:,:,tslice))','r.');
plot([min(min(su(:,:,tslice))):0.001:max(max(su(:,:,tslice)))], ...
     [min(min(su(:,:,tslice))):0.001:max(max(su(:,:,tslice)))],'k.');
ylabel('u','FontSize',fs);
set(gca,'FontSize',fs);
xlabel('u ageo','FontSize',fs);
hold off

subplot(2,2,4);
hold on
plot(flatten(vag)',flatten(sv)','b.');
plot(flatten(vag(:,:,tslice))',flatten(sv(:,:,tslice))','r.');
plot([min(min(sv(:,:,tslice))):0.001:max(max(sv(:,:,tslice)))], ...
     [min(min(sv(:,:,tslice))):0.001:max(max(sv(:,:,tslice)))],'k.');
ylabel('v','FontSize',fs);
xlabel('v ageo ','FontSize',fs);
set(gca,'FontSize',fs);
hold off

end

