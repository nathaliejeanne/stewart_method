%
% Beispiel aus Kapitel 3
%


clc;
clear all;

% Initialisierung
n = 1; % an Dimension von x anpassen
t = sym('t', [1 1]);
x = sym('x',[n 1]);

T = 2;
x0 = [-1];
f = [2 1];
h = [-x; x];
nablah = [-1; 1];

ffun = matlabFunction(f, 'Vars', [t; x]);
hfun = matlabFunction(h, 'Vars', [x]);
nablahfun = matlabFunction(nablah, 'Vars', [x]);

% Ausführung
[tsol, xsol, tswitch, tvector, psivector, zvector] = stewart(T, x0, ffun, hfun, nablahfun,length(h),@ode45,1);

% Plot
Schaltzeit = tswitch

% plot(tsol,xsol(1,:),'LineWidth',2); hold on
% plot(tswitch,zeros(1,length(tswitch)),'ko','MarkerSize',8); hold off
% xlabel('$t$','Interpreter','latex');
% % ylabel(['$x(t)$'],'Interpreter','latex');
% legend({'$x$','Schaltzeit'},'interpreter','latex','Location','northwest');
% grid on


plot(tvector, zvector(1,:),"--",'LineWidth',2); hold on
plot(tvector, zvector(2,:),"--",'LineWidth',2);
plot(tvector,psivector,'LineWidth',2); hold off
xlabel('$t$','Interpreter','latex');
% ylabel(['$x(t)$'],'Interpreter','latex');
legend({'$z_1$','$z_2$', '$\psi$'},'interpreter','latex');
grid on

   
  
%     legend('psi','x(t)',"x'(t)",'Schaltzeiten');


    set(gca,'fontsize', 14)
    set(groot,'defaultAxesTickLabelInterpreter','latex'); 
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
    set(gca,'fontsize', 20)
    
%     figure_name = ['bsp_stewart_psi.pdf'];
%     exportgraphics(gcf, figure_name, 'ContentType', 'vector');

