%
% Beispiel steife ODE
%


clc;
clear all;

% Initialisierung
n = 1; % an Dimension von x anpassen
t = sym('t', [1 1]);
x = sym('x',[n 1]);



x0 = 1;
f = [-1000*(x-cos(t)) 1000*(x-cos(t))];
h = [-x; x];
nablah = [-1; 1];



ffun = matlabFunction(f, 'Vars', [t; x]);
hfun = matlabFunction(h, 'Vars', [x]);
nablahfun = matlabFunction(nablah, 'Vars', [x]);



% Ausführung

T = 3;
count_start_stewart = cputime;
[tsol, xsol, tswitch, tvector, psivector, zvector] = stewart(T, x0, ffun, hfun, nablahfun,length(h),@ode15s,0);
Zeit_Stewart = cputime - count_start_stewart



count_start_disode = cputime;
 [tout,yout,tdis,ydis,idis,stats]=disode45(@fun_stiff, @gfun_stiff,[0,T], x0);
Zeit_DISODE = cputime - count_start_disode

% Plot

t1 = tout;
t2 = diff(tout);
h_stewart = diff(tsol);
h_disode = diff(tout);


Schaltzeiten_Stewart = tswitch
Schaltzeiten_DISODE = tdis



    plot(tswitch,zeros(1,length(tswitch)),'ko','MarkerSize',8); 
    hold on


%     plot(tsol,[nan;h_stewart],'LineWidth',2,'Color',[0.8500 0.3250 0.0980]);
%     plot(tout,[nan, h_disode],'LineWidth',2,'Color',[0.9290 0.6940 0.1250]);


    plot(tsol,xsol(1,:),'LineWidth',2);
    plot(tout,yout(:,1),'LineWidth',2);



    hold off
    
    
    axis([1.5688 1.5743 -0.01*10^(-3) 3*10^(-3)])


    legend({'Schaltzeit','Stewart',"DISODE45"},'location','northeast');
%     legend({'Stewart',"DISODE45"},'location','northwest');

    xlabel('$t$','Interpreter','latex');
%     ylabel(['diff(t)'],'Interpreter','latex');
    ylabel(['x(t)'],'Interpreter','latex');

    grid on
    set(gca,'fontsize', 14)
    set(groot,'defaultAxesTickLabelInterpreter','latex'); 
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
    set(gca,'fontsize', 20)
    
%     figure_name = ['stiff_ode.pdf'];
%     exportgraphics(gcf, figure_name, 'ContentType', 'vector');



% Funktionen benötigt für DISODE45 :


%
%  Definition of the vector field
%
 function f=fun_stiff(t,y)
      f = -1000 * (y - cos(t)) * sign(y);
 end
%
%  Definition of the switching surface
%
 function [g,isterminal,direction]=gfun_stiff(t,y)
      g = y;
      isterminal=0;
      direction=0;
 end