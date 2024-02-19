%
% Beispiel P1 aus Kapitel 5
%


clc;
clear all;

% Initialisierung
n = 2; % an Dimension von x anpassen
t = sym('t', [1 1]);
x = sym('x',[n 1]);


T = 20;

x0 = [1;1];
f = [x(2), x(2);-4*x(2)-x(1)-1,-4*x(2)-x(1)+1];
h = [-x(1)-x(2); x(1)+x(2)];
nablah = [-1, -1; 1, 1];



ffun = matlabFunction(f, 'Vars', [t; x]);
hfun = matlabFunction(h, 'Vars', [x]);
nablahfun = matlabFunction(nablah, 'Vars', [x]);


% Ausführung

count_start_stewart = cputime;
[tsol, xsol, tswitch, tvector, psivector, zvector] = stewart(T, x0, ffun, hfun, nablahfun,length(h),@ode45,0);
Zeit_Stewart = cputime - count_start_stewart


count_start_disode = cputime;
 [tout,yout,tdis,ydis,idis,stats]=disode45(@fun_p1, @gfun_p1,[0,T], x0);
Zeit_DISODE = cputime - count_start_disode


% Plot

Schaltzeiten_Stewart = tswitch
Schaltzeiten_DISODE = tdis



%     plot(tswitch,zeros(1,length(tswitch)),'ko','MarkerSize',8); 
        xline(tswitch,'-');
    hold on



    plot(tsol,xsol(1,:),'LineWidth',2);
    plot(tout,yout(:,1),'LineWidth',2);

    plot(tsol,xsol(2,:),'LineWidth',2);
    plot(tout,yout(:,2),'LineWidth',2);



    hold off
    
    xlim([0 T])
%     axis([2.5 4.2 -0.23 0.32])



    legend({'Schaltzeit','$x(t)$ Stewart',"$x(t)$ DISODE45","$x'(t)$ Stewart","$x'(t)$ DISODE45"},'location','northeast');

    xlabel('$t$','Interpreter','latex');
%     ylabel(['x(t)'],'Interpreter','latex');


    grid on
    set(gca,'fontsize', 14)
    set(groot,'defaultAxesTickLabelInterpreter','latex'); 
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
    set(gca,'fontsize', 20)
    
%     figure_name = ['p1_plot.pdf'];
%     exportgraphics(gcf, figure_name, 'ContentType', 'vector');



% Funktionen die DISODE45 benötigt : 

%
%  Definition of the vector field
%
 function f=fun_p1(t,y)
      f=[y(2);-y(1)-4*y(2)-sign(y(1)+y(2))] ;
%       f=[y(2);y(2)-4*y(1)-sign(y(1)+y(2))] ;
 end
%
%  Definition of the switching surface
%
 function [g,isterminal,direction]=gfun_p1(t,y)
      g=y(1)+y(2);
      isterminal=0;
      direction=0;
 end
