function plot_geostrophic(dplume, tslice, amb_temp,amb_salt)

[ug,vg,uag,vag,agfrac] = geostrophic_vel(dplume,amb_temp,amb_salt);

su = dplume.su;
sv = dplume.sv;
fs = 16;
m1 = 'b.';
m2 = 'b.';

if (false)

figure(1);
clf;

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

hold on
innerprod = ((su .* ug) + (sv .* vg))./(sqrt(ug.*ug+ vg.*vg).*sqrt(su.*su+sv.*sv));
proj = ((su .* ug) + (sv .* vg))./(sqrt(ug.*ug+ vg.*vg).*sqrt(su.*su+sv.*sv));

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

end


figure(3);
clf;
m1 = 'k.';
[m,n,k] = size(su);
flatten = @(M) reshape(squeeze(M(2:(end-1),2:(end-1),:)), ...
                              (size(M,1)-2)*(size(M,2)-2)*size(M,3),1);

minsu = min([min(min(su(:,:,tslice))), min(min(ug(:,:,tslice))),min(min(uag(:,:,tslice)))]);
maxsu = max([max(max(su(:,:,tslice))), max(max(ug(:,:,tslice))),max(max(uag(:,:,tslice)))]);
minsv = min([min(min(sv(:,:,tslice))), min(min(vg(:,:,tslice))),min(min(vag(:,:,tslice)))]);
maxsv = max([max(max(sv(:,:,tslice))), max(max(vg(:,:,tslice))),max(max(vag(:,:,tslice)))]);

subplot(2,2,1);
hold on
xlim([minsu maxsu]);
ylim([minsu maxsu]);
%plot(flatten(ug)',flatten(su)',m1);
plot(flatten(ug(:,:,tslice))',flatten(su(:,:,tslice))',m2);
plot(minsu:0.001:maxsu,minsu:0.001:maxsu,m1);
ylabel('x-vel.','FontSize',fs);
xlabel('geostrophic x-vel.','FontSize',fs);
set(gca,'FontSize',fs);
hold off

subplot(2,2,2);
hold on
xlim([minsv maxsv]);
ylim([minsv maxsv]);
%plot(flatten(vg)',flatten(sv)',m1);
plot(flatten(vg(:,:,tslice))',flatten(sv(:,:,tslice))',m2);
plot(minsv:0.001:maxsv,minsv:0.001:maxsv,m1);
ylabel('y-vel.','FontSize',fs);
xlabel('geostrophic y-vel.','FontSize',fs);
set(gca,'FontSize',fs);
hold off

subplot(2,2,3);
hold on
xlim([minsu maxsu]);
ylim([minsu maxsu]);
%plot(flatten(uag)',flatten(su)',m1);
plot(flatten(uag(:,:,tslice))',flatten(su(:,:,tslice))',m2);
plot(minsu:0.001:maxsu, minsu:0.001:maxsu,m1);
ylabel('x-vel.','FontSize',fs);
set(gca,'FontSize',fs);
xlabel('ageostrophic x-vel.','FontSize',fs);
hold off

subplot(2,2,4);
hold on
xlim([minsv maxsv]);
ylim([minsv maxsv]);
%plot(flatten(vag)',flatten(sv)',m1);
plot(flatten(vag(:,:,tslice))',flatten(sv(:,:,tslice))',m2);
plot([minsv:0.001:maxsv],[minsv:0.001:maxsv],m1);
ylabel('y-vel.','FontSize',fs);
xlabel('ageostrophic y-vel. ','FontSize',fs);
set(gca,'FontSize',fs);
hold off

end

