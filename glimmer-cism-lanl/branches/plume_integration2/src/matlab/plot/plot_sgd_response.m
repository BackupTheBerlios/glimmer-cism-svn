%plot sgd response
lw = 1.0;
fs = 16;
ms = 10.0;

x = [0 0.5 1.0 1.5 2.0];
y = [0 0.1 0.19 0.24 0.98];
%x = [0 0.5 1.0 1.5];
%y = [0 0.1 0.19 0.24];

figure(1);
clf;
plot(x,y,'k*','MarkerSize',ms);

p = polyfit(x(1:4),y(1:4),1);

x_line = -0.125:0.1:2.25;
hold on

plot(x_line, p(1)*x_line +p(2),'b--','linewidth',lw);
xlim([-0.125 2.25]);

xlabel('Subglacial discharge rate (km^3 a^{-1})','FontSize',fs);
ylabel('Increase in total melting (km^3 a^{-1})','FontSize',fs);
set(gca,'FontSize',fs);

text(1.125,0.15,'slope = 0.16','FOntSize',fs);
text(1.25,1,'not steady','FontSize',fs);
