%Часть без анимации reachset
alpha = 3;%2,3,4,5
t = 0.85;
[X,Y,x_l,y_l] = reachset(alpha, t);
hold on
plot(X,Y,'Color','b','LineWidth',3);
plot(x_l,y_l,'Color','r','LineWidth',3);
xlabel('x_1');
ylabel('x_2');
hold off
%%
%Часть с анимацией reachsetdyn

alpha = 3;
t_1 = 0.3;
t_2 = 0.5;
filename = 'name.avi';
N = 25;
reachsetdyn(alpha,t_1,t_2,N,filename);

%%
%Функции-----------------------------------------

function [out1, out2, out3, out4] = reachset(alpha, t_1)
    t_0 = 0;
    tspan = [t_0 t_1];
    options = odeset('Events',@eventfcn);
    options_psi = odeset('Events',@eventfcn_psi);
    [t_plus, x_plus] = ode45(@(t_plus,x_plus) odefcn_s_plus(t_plus,x_plus,alpha), tspan, [0, 0],options);
    hold on
   % plot(x_plus(:,1),x_plus(:,2),'Color','black');
    hold off
    n = numel(x_plus(:,2));
    XM = [];
    YM = [];
    if (t_plus(n) <= t_1)
       tau = 0; 
    else
       tau = t_plus(n); 
    end    
    alph = alpha;
    x_T = [];
    y_T = [];
    for i=1:n-1
        b = 1;
        alpha = alph;
        x_1_b = x_plus(i,1);
        x_2_b = x_plus(i,2);
        tspan = [tau,t_1];
        if(tau == t_1)
            b = 0;
        end 
        W_plus_x = [];
        W_plus_y = [];
        W_minus_x = [];
        W_minus_y = [];
        while b
            alpha = alpha*(-1);
            j = -sign(alpha);
           
            [t_psi,x_psi] = ode45(@(t_psi,x_psi) odefcn_psi(t_psi,x_psi,alpha),tspan,[x_1_b x_2_b j 0], options_psi);
            n_t = numel(t_psi);
            t_psi_e = t_psi(n_t);
            hold on 
            XM = [XM, x_psi(:,1)'];
            YM = [YM, x_psi(:,2)'];
            
            hold off
            
            if t_psi_e == t_1
                b = 0;
                
                x_T = [x_T, x_psi(n_t,1)];
                y_T = [y_T, x_psi(n_t,2)];
            else
                if sign(alpha) < 0
                    
                    W_minus_x = [W_minus_x, x_psi(n_t,1)];
                    W_minus_y = [W_minus_y, x_psi(n_t,2)]; 
                else
                    
                    W_plus_x = [W_plus_x, x_psi(n_t,1)];
                    W_plus_y = [W_plus_y, x_psi(n_t,2)];
 
                end
                
                tspan = [t_psi_e, t_1];
                x_1_b = x_psi(n_t,1);
                x_2_b = x_psi(n_t,2);
            end
        end    
    end
    
    XM = [XM, x_T];
    YM = [YM, y_T];
%---------------------------------------------------------------------------------------------------------------------------------------------------------
    options_psi = odeset('Events',@eventfcn_psi);
    tspan = [t_0 t_1];
    [t_minus, x_minus] = ode45(@(t_minus,x_minus) odefcn_s_plus(t_minus,x_minus,alpha), tspan, [0, 0],options);
    hold on
    %plot(x_minus(:,1),x_minus(:,2),'Color','red');
    hold off
    n = numel(x_minus(:,2));
    if (t_minus(n) <= t_1)
       tau = 0; 
    else
       tau = t_minus(n); 
    end    
    x_T = [];
    y_T = [];
    for i=1:n-1
        b = 1;
        alpha = -alph;
        x_1_b = x_minus(i,1);
        x_2_b = x_minus(i,2);
        tspan = [tau,t_1];
        if(tau == t_1)
            b = 0;
        end
        W_plus_x = [];
        W_plus_y = [];
        W_minus_x = [];
        W_minus_y = [];
        while b
            alpha = alpha*(-1);
            j = -sign(alpha);
           
            [t_psi,x_psi] = ode45(@(t_psi,x_psi) odefcn_psi(t_psi,x_psi,alpha),tspan,[x_1_b x_2_b j 0], options_psi);
            n_t = numel(t_psi);
            t_psi_e = t_psi(n_t);
            XM = [XM, x_psi(:,1)'];
            YM = [YM, x_psi(:,2)'];
            
            if t_psi_e == t_1
                b = 0;
                
                x_T = [x_T, x_psi(n_t,1)];
                y_T = [y_T, x_psi(n_t,2)];
            else
                if sign(alpha) < 0
                    W_minus_x = [W_minus_x, x_psi(n_t,1)];
                    W_minus_y = [W_minus_y, x_psi(n_t,2)];
                    
                else
                    W_plus_x = [W_plus_x, x_psi(n_t,1)];
                    W_plus_y = [W_plus_y, x_psi(n_t,2)];
 
                end
                
                tspan = [t_psi_e, t_1];
                x_1_b = x_psi(n_t,1);
                x_2_b = x_psi(n_t,2);
            end
            
        end      
    end
    XM = [XM, x_T];
    YM = [YM, y_T];
    XM = XM';
    YM = YM';
    k = boundary(XM,YM,0.3);
    out1 = XM(k);
    out2 = YM(k);
    out3 = [fliplr(x_minus(:,1)),x_plus(:,1)];
    out4 = [fliplr(x_minus(:,2)),x_plus(:,2)];
end

function reachsetdyn(alpha, t1, t2, N, filename)
    mov(1:N) = struct('cdata',[],'colormap',[]);
    for i=1:N
        tau = t1 + ((t2 - t1)*(i))./N;
        [x, y, x1, y1] = reachset(alpha, tau);
        
        plot(x,y,'Color','b','LineWidth',3);
        axis([-2 2 -2 2]);%фиксация
        mov(i) = getframe();
    end
    des = VideoWriter('a.avi');
    open(des);
    writeVideo(des,mov);
    close(des);
end

function dydt = odefcn_s_plus(t,y,alpha)
    dydt = zeros(2,1);
    dydt(1) = y(2);
    dydt(2) = alpha - y(1).^2 + 2.*sin(3.*y(1).^3) - y(1).*y(2);
end

function dydt = odefcn_psi(t,y,alpha)
    dydt = zeros(4,1);
    dydt(1) = y(2);
    dydt(2) = alpha - y(1).^2 + 2.*sin(3.*y(1).^3) - y(1).*y(2);
    dydt(3) = 2*y(4)*y(1) - 18*y(4)*((y(1))^(2))*cos(3*y(1)^3) + y(4)*y(2);
    dydt(4) = -y(3) + y(4)*y(1);
end

function [value, isterminal, direction] = eventfcn(t,x)
    value = x(2);
    isterminal = 1;
    direction = 0;
end

function [value, isterminal, direction] = eventfcn_psi(t,x)
    value = x(4);
    isterminal = 1;
    direction = 0;
end