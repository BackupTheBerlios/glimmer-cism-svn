fs = 14;
clf;

%subplot(2,1,1);

hold on
plot(amb_t(-flat_ocean.draft),flat_ocean.temp,'b.');
plot(amb_t(-flat_ocean.draft),amb_t(-flat_ocean.draft),'r.');
ylim([-1 0.1]);
xlim([-1 0.1]);
xlabel('Ambient water temperature(C)','FontSize',fs);
ylabel('Plume water temperature (C)','FontSize',fs);      
title('Plume temperatures against ambient water column temperature at ice base','FontSize',fs);
hold off

%subplot(2,1,2);
%plot(-flat_ocean.draft,flat_ocean.temp,'b.');
%ylim([-1 0.1]);
%xlabel('Ice draft (m)','FontSize',fs);
%ylabel('Plume water temperature (C)','FontSize',fs);      
%title('Plume temperatures against ambient water column temperature at ice base','FontSize',fs);
