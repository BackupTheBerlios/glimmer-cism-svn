clf;
hold on
plot(d.x/1000.0,d.draft(:,1,808),'k')
text(12,-575,'0km');
plot(d.x/1000.0,d.draft(:,6,808),'r')
text(11.5,-525,'2.5km');
plot(d.x/1000.0,d.draft(:,11,808),'b')
text(11.5,-475,'5.0km')
plot(d.x/1000.0,d.draft(:,16,808),'g')
text(11.0, -425,'7.5km')
%plot(d.x/1000.0,d.draft(:,21,808),'k')
plot(d.x/1000.0,d.draft(:,40,808),'k')
text(12,-350,'20.0km')
%plot(d.x/1000.0,d.draft(:,60,808),'k')
plot(d.x/1000.0,d.draft(:,80,808),'k')
text(12.5,-250,'40km');
%plot(d.x/1000.0,d.draft(:,100,808),'k')
plot(d.x/1000.0,d.draft(:,120,808),'k')
text(12.0,-125,'60km');
fs = 14;
xlabel('Cross-shelf distance (km)','FontSize',fs);
ylabel('Ice Draft (m)','FontSize',fs);
title('Cross-shelf ice draft sections','FontSize',fs);
