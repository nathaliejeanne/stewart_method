%
% Beispiel P2 aus Kapitel 5
%


clc;
clear all;

% Initialisierung
n = 2; % an Dimension von x anpassen
t = sym('t', [1 1]);
x = sym('x',[n 1]);




T = 10;

x0 = [3;4];
f = [x(2), x(2);-0.2*x(2)-x(1)-4+2*cos(pi*t),-0.2*x(2)-x(1)+4+2*cos(pi*t)];
h = [-x(2); x(2)];
nablah = [0, -1; 0, 1];



ffun = matlabFunction(f, 'Vars', [t; x]);
hfun = matlabFunction(h, 'Vars', [x]);
nablahfun = matlabFunction(nablah, 'Vars', [x]);


% Ausführung

count_start_stewart = cputime;
[tsol, xsol, tswitch, tvector, psivector, zvector] = stewart(T, x0, ffun, hfun, nablahfun,length(h),@ode45,0);
Zeit_Stewart = cputime - count_start_stewart


count_start_disode = cputime;
 [tout,yout,tdis,ydis,idis,stats]=disode45(@fun_p2, @gfun_p2,[0,T], x0);
Zeit_DISODE = cputime - count_start_disode


% Plot

Schaltzeiten_Stewart = tswitch
Schaltzeiten_DISODE = tdis



%     plot(tswitch,zeros(1,length(tswitch)),'ko','MarkerSize',8); 
        xline(tswitch,'-');
    hold on


    yyaxis left;
    plot(tsol,xsol(1,:),"-.",'LineWidth',2);
    plot(tout,yout(:,1),"-",'LineWidth',2);
    ylim([2.5, 4.5]);

    yyaxis right;
    plot(tsol,xsol(2,:),"-.",'LineWidth',2);
    plot(tout,yout(:,2),"-",'LineWidth',2);
    ylim([-1.5, 1.5]);

    xlim([0, 3.5]);



    hold off
    

%     axis([0 10 -2 6])
%     axis([1.48 2.52 -1.3 4.1])
%     axis([0 4 -1.5 6.2])


%     legend({'Schaltzeit','Stewart',"DISODE45"},'location','northeast');
    legend({'Schaltzeiten','','','','','','','','','','$x(t)$ Stewart',"$x(t)$ DISODE45","$x'(t)$ Stewart","$x'(t)$ DISODE45"},'location','northeast');

    xlabel('$t$','Interpreter','latex');
%     ylabel(['x(t)'],'Interpreter','latex');
%     ylabel(['Zeitschritte'],'Interpreter','latex');

    grid on
    set(gca,'fontsize', 14)
    set(groot,'defaultAxesTickLabelInterpreter','latex'); 
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
    set(gca,'fontsize', 20)
    
%     figure_name = ['p2_small_2.pdf'];
%     exportgraphics(gcf, figure_name, 'ContentType', 'vector');


% Funktionen, die für DISODE45 benötigt werden :

% P2 : 

%
%  Definition of the vector field
%
 function f=fun_p2(t,y)
      f=[y(2);-0.2*y(2)-y(1)-4*sign(y(2))+2*cos(pi*t)] ;
 end
%
%  Definition of the switching surface
%
 function [g,isterminal,direction]=gfun_p2(t,y)
      g=y(2);
      isterminal=0;
      direction=0;
 end
