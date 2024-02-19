%
% Implementierung der Stewart-Methode
%


% Input:

% T : Intervalllänge
% x0 : Anfangswert
% f : Funktion f (matlabFunction)
% h : Funktion h, beinhaltet alle Diskriminanzfunktionen (matlabFunction)
% nabla_h : Funktion nabla_h, Gradient von h (matlabFunction)
% length_of_h : Anzahl der Diskriminanzfunktionen 
% odesolver : Wähle einen odesolver (ode45, ode15s, ...)
% psi_z_plot : 1 falls Vektoren t_vector, psi_vector und z_vector berechnet
% werden sollen um z und psi zu plotten


% Output:

% t_end : Zeitvektor der Lösung
% x_end : Lösungsvektor
% t_switch : Schaltzeiten
% t_vector : benötigt für Plots von psi und z
% psi_vector : Vektor, der Werte für psi enthält
% z_vector : Vektor, der Werte für z enthält



function [t_end, x_end, t_switch, t_vector, psi_vector, z_vector] = stewart(T,x_0, f, h, nabla_h,length_of_h,odesolver,psi_z_plot)
    % Dimension von x
    n = length(x_0);

    % Initialisierung
    eps = 10^(-5);
    t_star = 10^2;
    tau = 0.5; % Schrittweite
    alpha = 10^5;

    t_step = 0;
    s = 0;
    x_step = x_0;
    y = x_0;
    I = I_epsilon(x_step(:,1)); % Startmenge
    psi_vec(1) = psi(I, x_step(:,1), t_step(1), t_star); % psi 0 initialisieren

    t_switch = []; % Schaltzeiten

    % Symbolische Variablen
    t = sym('t', [1 1]);
    x = sym('x',[n 1]);
    
    % Werte von phi und z um plotten zu können
    z_vector = [];
    psi_vector = [];
    t_vector = [];
    % Array der am Ende alle Zeitschritte t und deren Lösung x enthalt
    t_end = [];
    x_end = [];

    

    % Start

    k = 1;

    while t_step(k)<T

       s = cat(2,s,t_step(k) + tau);


       % function handle rhsI mit neuem I erstellen

       M = sym(zeros(length(I)));
       nabla_h_fun = subs(nabla_h);
       f_fun = subs(f);
       for i = 1:length(I)
           for j = 1:length(I)
               m_value = nabla_h_fun(I(i),:)*f_fun(:,I(j));
               M(i,j) = m_value;
           end
       end
       M_alpha1 = M + alpha.*ones(size(I));

       e = ones(length(I),1);
       z_hat = inv(M_alpha1)*e;
       z1 = z_hat./(transpose(e)*z_hat);
       
       z_fun = matlabFunction(z1, 'Vars', {t; x}); % z_fun wird später für z_vector benötigt

       sum = 0;
       for i = 1:length(I)
            sum = sum + z1(i)*f_fun(:,I(i));
       end


       % function handle der rhs
       rhsI = matlabFunction(sum, 'Vars', {t; x});

       % ODE Solver gibt psi-Nullstelle in letztem Eintrag te und ye von t_vec, y_vec zurück.
       options = odeset('Events', @psiEvent);
       [t_vec,y_vec, te, ye, ~] = odesolver(rhsI, [t_step(k), s(k+1)], x_step(:,k), options);
    
       y = cat(2,y,transpose(y_vec(end,:)));
       

       psi_vec = cat(2,psi_vec,psi(I,y(:,k+1),s(k+1), t_star));
       

        
       % Existiert Schaltzeit? : 
        
       if ~isempty(te)
            if (psi_vec(k+1) <= 0) && (psi_vec(k+1) < psi_vec(k))

            t_step = cat(2,t_step,te);
            t_switch = cat(2,t_switch,t_step(k+1));
            x_step = cat(2,x_step,transpose(ye));
            I_eps = I_epsilon(x_step(:,k+1)); % neues I_eps bestimmen

           % Bestimmung neuer aktiver Menge :

            Ma = M_alpha(I_eps, x_step(:,k+1),t_step(k+1));
            zp = pathlcp(Ma, -ones(1, size(Ma, 1)));
            Ip = [];
            for i = 1:length(zp)
                if zp(i) > 0
                   Ip = cat(2, Ip, i); 
                end
            end
            I = Ip;

            psi_vec(k+1) = psi(I,x_step(:,k+1),s(k+1), t_star); % psi_vec(k+1) an neues I anpassen
            
            end
       else 
            % Schritt akzeptiert falls keine Schaltzeit

            x_step = cat(2,x_step,y(:,k+1));
            t_step = cat(2,t_step,s(k+1));
            
            t_end = cat(1, t_end, t_vec); 
            x_end = cat(2, x_end, transpose(y_vec)); 

            % Sollen Werte für psi und z zum plotten bestimmt werden?

            if psi_z_plot == 1

            t_vector = cat(1,t_vector,t_vec);

            % for-Schleife zur Bestimmung von psi_vector und z_vector
            
            for v = 1:length(t_vec)
               psi_vector = cat(2,psi_vector,psi(I,transpose(y_vec(v,:)),t_vec(v), t_star));
                
               z_vec = z_fun(t_vec(v), transpose(y_vec(v,:)));

               % z_var (erzeugt durch for-Schleife) wird benötigt, da z_i
               % für i nicht in I nicht bestimmt wird, z_i muss daher Null
               % gesetzt werden.
               z_var = [];
               l = 1;
               for m = 1:length_of_h
                       if ~ismember(m, I)
                           z_var = cat(1,z_var,0); 
                       else
                           z_var = cat(1,z_var,z_vec(l));
                           l = l+1;
                       end
               end
               z_vector = cat(2,z_vector,z_var);
            end  
            end
            
       end

       k = k+1;

    end


   
    % STOP bei psi = 0

    function [value, isterminal, direction] = psiEvent(s, y)
        value = psi(I, y, s, t_star); % Stop bei psi = 0
        isterminal = 1;
        direction = -1;
    end
    
    function z = z(I, y,s)
        e = ones(length(I),1);
        M_alphax = M_alpha(I,y,s);
        z_hat = M_alphax\e;
        z = z_hat./(transpose(e)*z_hat);
    end

    function M_alpha = M_alpha(I,y,s)
        M = zeros(length(I));
        y_cell = num2cell(y);
        nabla_hy = nabla_h(y_cell{:});
        fy = f(s,y_cell{:});
        for i = 1:length(I)
            for j = 1:length(I)
                M(i,j) = nabla_hy(I(i),:)*fy(:,I(j));
            end
        end
        alpha = abs(min(M(:))) + 1;
        M_alpha = M + alpha.*ones(size(I));
    end
    
    function psi = psi(I, y, s, t_star) 
        h_without_I = [];
        h_with_I = [];
        y_cell = num2cell(y);
        hy = h(y_cell{:});
        for i = 1:length(hy)
            if ~ismember(i, I)
                h_without_I = cat(2, h_without_I, hy(i));
            else
                h_with_I = cat(2, h_with_I, hy(i));
            end
        end
        psi = min([min(z(I, y,s)), min(h_without_I) - min(h_with_I), t_star - s]);
    end

    function I_eps = I_epsilon(y)
        y_cell = num2cell(y);
        hy = h(y_cell{:});
        I_eps = [];
        for i = 1:length(hy)
            if hy(i) < min(hy) + eps
                I_eps = cat(2, I_eps, i);
            end
        end
    end
end