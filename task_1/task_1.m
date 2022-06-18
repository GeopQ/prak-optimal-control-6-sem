%%
%Входящие параметры.
 
A = [0,1;-1, 0];
B = [5,0 ;0, 5];
f = [2, -2];

t0 = 0;
t_max = 1.5;

p1 = 0;
p2 = 0;

alpha = 1;
beta = 1;
gamma = 1;

q1 = [-10, 3];
q2 = [-8,3];
q3 = [-5.6, 2.6];

x1 = [3,0];

r = 1/2;
if ~X1_in_X0(x1,q1,q2,q3,r)        
%-------------------------------------
%Произведение замены для решения задачи в обратном времени.
C = -A;
D = -B;
g = -f;

if abs(det(D)) < 0.1
    D = D + 0.1*eye(2);
end    
C_T = C';
%Основной цикл программы.
n = 100;
arg_sc = linspace(0,2*pi - 2*pi./n,n);
T_n = zeros(1,n);
figure('Name','График оптимальной траектории X2(X1)','NumberTitle','off');
hold on
draw_X0(q1,q2,q3,r);
hold off
x_opt_mas = [];
y_opt_mas = [];
T_mas = [];
T = t_max;
t_un = linspace(t0, t_max, 100);
for i = 1:n
    w_0 = [sin(arg_sc(i)),cos(arg_sc(i))];
    tspan = [t0 t_max];
    [t,y] = ode45(@(t,y) odefcn_w(t,y,-C_T),tspan,w_0);
    t_leng = length(t);
    u_opt = zeros(2,t_leng);
    u_value = zeros(1,t_leng);
    
    for j = 1:t_leng
        M_supp = D';
        v_supp = M_supp*y(j,:)';
        inf_m = surf_P(v_supp,alpha,beta,gamma,p1,p2);
        value = inf_m(1);
        cords = [inf_m(2),inf_m(3)];
        u_opt(1,j) = cords(1);
        u_opt(2,j) = cords(2);
        u_value(j) = value;
    end
    tspan_1 = [t0 t_max];
    options = odeset('Events',@(t,y) EventsFcn(t,y,q1,q2,q3,r),'RelTol',1e-5,'AbsTol',1e-6);
    [t_opt,psi_opt,te,ye,ie] = ode45(@(t_opt, psi_opt) odefcn_OC(t_opt,psi_opt,C,D,g,t,u_opt), tspan_1, x1,options);    
    hold on
    %plot(psi_opt(:,1),psi_opt(:,2),'r');
    hold off
    
    U_x = interp1(t_opt,psi_opt(:,1),t_un);
    U_y = interp1(t_opt,psi_opt(:,2),t_un);
    x_opt_mas = [x_opt_mas, U_x'];
    y_opt_mas = [y_opt_mas, U_y'];
    T_mas = [T_mas, max(t_opt)];
    if max(t_opt) <= T
        x_help = psi_opt(:,1);
        y_help = psi_opt(:,2);
        u_opt_gr_x = interp1(t,u_opt(1,:),t_un);
        u_opt_gr_y = interp1(t,u_opt(2,:),t_un);
        psi_opt_gr_x = interp1(t,y(:,1),t_un);
        psi_opt_gr_y = interp1(t,y(:,2),t_un);
        psi_help_x = y(:,1);
        psi_help_y = y(:,2);
    end
end
%Вывод графиков из подзадач.
[m,i] = min(T_mas);
if abs(m - t_max) < 0.1
    disp('Задача не разрешима за время t_max.');
    b = 0;
    flag = 0;
else
    disp('Задача разрешима.');
    disp('Время T:');
    disp(m);
    b = 1;
    flag = 1;
end
for j = 1 : n
   if (abs(j - i) <= 0.01 && b)
     hold on
     plot(x_opt_mas(:,j),y_opt_mas(:,j),'r');
     hold off   
   else
     hold on
     plot(x_opt_mas(:,j),y_opt_mas(:,j),'b');
     hold off   
   end
end
k = 1;
b = 1;
while ((x_opt_mas(k,i) > 0) || (x_opt_mas(k,i) <= 0)) && b
    k = k + 1; 
    if  k >= size(x_opt_mas(:,i))
        b = 0; 
    end    
end
xlabel('x');
ylabel('y');
%figure('Name','График компоненты оптимальной траектории X_x(t)','NumberTitle','off');
%plot(t_un,x_opt_mas(:,i));
%xlabel('t');
%ylabel('x(t)');
%figure('Name','График компоненты оптимальной траектории X_y(t)','NumberTitle','off');
%plot(t_un,y_opt_mas(:,i));
%xlabel('t');
%ylabel('y(t)');
%figure('Name','График компоненты оптимального управления U_x(t)','NumberTitle','off');
%plot(t_un,u_opt_gr_x);
%xlabel('t');
%ylabel('u_x(t)');
%figure('Name','График компоненты оптимального управления U_y(t)','NumberTitle','off');
%plot(t_un,u_opt_gr_y);
%xlabel('t');
%ylabel('u_y(t)');
%figure('Name','График компоненты решения сопряженной системы psi_x(t)','NumberTitle','off');
%plot(t_un,psi_opt_gr_x);
%xlabel('t');
%ylabel('psi_x(t)');
%figure('Name','График компоненты решение сопряженной системы psi_y(t)','NumberTitle','off');
%plot(t_un,psi_opt_gr_y);
%xlabel('t');
%ylabel('psi_y(t)');
%Расчет погрешности условий трансверсальности.
if flag
    usl_t_l_1 = dot([psi_opt_gr_x(1),psi_opt_gr_y(1)],[x_opt_mas(1,i),y_opt_mas(1,i)]);
    usl_t_l_2 = supf_X1([psi_opt_gr_x(1),psi_opt_gr_y(1)],x1);

    disp('Погрешность условия трансверсальности на левом конце:');
    disp(abs(usl_t_l_1 - usl_t_l_2(1)));
    
    usl_t_r_1 = dot([-psi_opt_gr_x(end),-psi_opt_gr_y(end)],[x_opt_mas(k-1,i),y_opt_mas(k-1,i)]);
    usl_t_r_2 = supf_X0([-psi_opt_gr_x(end),-psi_opt_gr_y(end)],q1,q2,q3,r);
    disp('Погрешность условия трансверсальности на правом конце:');
    disp(abs(usl_t_r_1 - usl_t_r_2(1)));   
end
else
    disp('x1 расположен в целевом множестве, поэтому ничего не считаем.');
end
%%

function dydt = odefcn_w(t,y,C)
    dydt = zeros(2,1);
    dydt(1) = C(1,1).*y(1) + C(1,2).*y(2);
    dydt(2) = C(2,1).*y(1) + C(2,2).*y(2);
end

function res = surf_P(l,a,b,g,p1,p2)
    res = zeros(3,1);
    if l(1) >= 0
        res(1) = l(1).*(p1 + l(1)*sqrt((a)./(g*(l(1)^2 + l(2)^2)))) + l(2).*(p2 + l(2)*sqrt((a)./(g*(l(1)^2 + l(2)^2)))); 
        res(2) = p1 + l(1)*sqrt((a)./(g*(l(1)^2 + l(2)^2))); 
        res(3) = p2 + l(2)*sqrt((a)./(g*(l(1)^2 + l(2)^2)));
    else
        res(1) = l(1).*(p1 + sqrt((b.*l(1)^2)./(g*(l(1)^2 + l(2)^2)))) + l(2).*(p2 + l(2)*sqrt((b)./(g*(l(1)^2 + l(2)^2))));
        res(2) = p1 + l(1)*sqrt((b)./(g*(l(1)^2 + l(2)^2)));
        res(3) = p2 + l(2)*sqrt((b)./(g*(l(1)^2 + l(2)^2)));
    end       
    %res(1) = norm(l);
    %res(2) = l(1)/res(1);
    %res(3) = l(2)/res(1);
end

function res = supf_X0(l,q1,q2,q3,r)
    res = zeros(3,1);
    res(1) = max(max(dot(l,q1),dot(l,q2)), dot(l,q3)) + abs(r);
    if (dot(l,q1) > dot(l,q2))|| (dot(l,q1) > dot(l,q3))
        res(2) = q1(1) + r*l(1); 
        res(3) = q2(2) + r*l(2);  
    else
        if (dot(l,q2) > dot(l,q3))
            res(2) = q2(1) + r*l(1);
            res(3) = q3(2) + r*l(2);  
        else    
        res(2) = q3(1) + r*l(1);
        res(3) = q3(2) + r*l(2);
        end
    end    
end

function res = supf_X1(l,x)
    res(1) = dot(l,x);
    res(2) = x(1);
    res(3) = x(2);
end

function dphidt = odefcn_OC(t,phi,A,B,f,setka,u)
    dphidt = zeros(2,1);
    
    u_1 = zeros(1,length(setka));
    u_2 = zeros(1,length(setka));

    u_1(1,1:length(setka))= B(1,1)*u(1,1:length(setka)) + B(1,2)*u(1:length(setka));
    u_2(1,1:length(setka))= B(2,1)*u(1,1:length(setka)) + B(2,2)*u(1:length(setka));
    
    U_1 = interp1(setka,u_1,t);
    U_2 = interp1(setka,u_2,t);
    dphidt(1) = A(1,1)*phi(1) + A(1,2)*phi(2) + f(1) + U_1;
    dphidt(2) = A(2,1)*phi(1) + A(2,2)*phi(2) + f(2) + U_2;
end

function [value, isterminal, direction] = EventsFcn(t,y,q1,q2,q3,r)
    y = y';
    
    q1_q2 = supp_line(q1,q2);
    q1_q3 = supp_line(q1,q3);
    q2_q3 = supp_line(q2,q3);
    r_1 = supp_dist(y,q1_q2);
    r_2 = supp_dist(y,q1_q3);
    r_3 = supp_dist(y,q2_q3);
    
    if r == 0
        r = 0.01;
    end    
    r_1_b = r_1 < r;
    r_2_b = r_2 < r;
    r_3_b = r_3 < r;
    
    m_b = max(max(max(r_1_b),max(r_2_b)),max(r_3_b));
    b = 0;
    num_1 = (q1(1) - y(1))*(q2(2) - q1(2)) - (q2(1) - q1(1))*(q1(2) - y(2));
    num_2 = (q2(1) - y(1))*(q3(2) - q2(2)) - (q3(1) - q2(1))*(q2(2) - y(2));
    num_3 = (q3(1) - y(1))*(q1(2) - q3(2)) - (q1(1) - q3(1))*(q3(2) - y(2));
    if (num_1 > 0) && (num_2 > 0) && (num_3 > 0)
        b = 1; 
    end
     if (num_1 < 0) && (num_2 < 0) && (num_3 < 0)
        b = 1; 
     end
     
     if  b || m_b
        value = 0; 
    else
        value = 1;
    end      
    isterminal = 1;
    direction = 0;  
end

function res_m = supp_line(q1,q2)
    n = 100;
    if q1(1) == q2(1)
       if q1(2) == q2(2)
            res_m = [q1(1);q1(2)];
       else
           x = q1(1).*ones(1,n);
           y = linspace(q1(2),q2(2),n);
           res_m = [x;y];
       end    
    else
       if q1(2) == q2(2)
           y = q1(2).*ones(1,n);
           x = linspace(q1(1),q2(1),n);
           res_m = [x;y];
       else
           setka = linspace(q1(1),q2(1),n);
           x = setka;
           y = ((q2(2) - q1(2)).*x +q2(1).*q1(2) - q1(1).*q2(2))./(q2(1) - q1(1));
           res_m = [setka;y];
       end    
    end    
end

function res_m = supp_dist(y,q1_q2)
    l = length(q1_q2(1,:));
   
    mas_s = [y(1).*ones(1,l); y(2).*ones(1,l)];
    
    q1_q2 = q1_q2 - mas_s;
   
    q1_q2 = q1_q2.^2;
   
    sum_pr = q1_q2(1,:) + q1_q2(2,:);
   
    res_m = sum_pr.^(1/2);
    
end

function res = draw_X0(q1,q2,q3,r)
    Q = [q1;q2;q3];
    if (q1(1) == q2(1))&&(q1(1) == q3(1))||(q1(2) == q2(2))&&(q1(2) == q3(2))
        
        rank_m = 1;
    else    
        rank_m = rank(Q);
    end
    switch rank_m
        case 1
             if r == 0
                 if ((max(q1 - q2) == 0) && (max(q1 - q3)==0))
                     hold on
                     plot(q1(1),q1(2),'g.','MarkerSize',10);
                     hold off
                 else    
                 hold on
                 plot([q1(1),q2(1)],[q1(2),q2(2)],'Color','green');
                 plot([q1(1),q3(1)],[q1(2),q3(2)],'Color','green');
                 hold off
                 end
             else
                 
                 
             q1_q2 = supp_line(q1,q2);
             q1_q3 = supp_line(q1,q3);
             q2_q3 = supp_line(q2,q3);
             
             n = 3;
             arg = linspace(0,2*pi - 2*pi./n);
             x = [];
             y = [];
             for i=1:n
                 x = [x, r.*sin(arg(i))];
                 y = [y, r.*cos(arg(i))]; 
             end
             x = [x,x(1)];
             y = [y,y(1)];
             n_l = length(q1_q2(1,:));
             hold on
             supp_draw(x,y,q1_q2,n_l);
             supp_draw(x,y,q1_q3,n_l);
             supp_draw(x,y,q2_q3,n_l);
             hold off
             end   
         case 2
            if r == 0
                P = [q1;q2;q3]; 
                [k, av] = convhull(P);
                hold on
                fill(P(k,1),P(k,2),'green');
                plot(P(k,1),P(k,2),'Color','green');
                hold off
            else
               q1_q2 = supp_line(q1,q2);
               q1_q3 = supp_line(q1,q3);
               q2_q3 = supp_line(q2,q3);
               n = 100;
               arg = linspace(0,2*pi - 2*pi./n);
               x = [];
               y = [];
               for i=1:n
                   x = [x, r.*sin(arg(i))];
                   y = [y, r.*cos(arg(i))]; 
               end
               x = [x,x(1)];
               y = [y,y(1)];
               n_l = length(q1_q2(1,:));
               hold on
               P = [q1;q2;q3];
               [k, av] = convhull(P);
               fill(P(k,1),P(k,2),'green');
               plot(P(k,1),P(k,2),'Color','green');
               supp_draw(x,y,q1_q2,n_l);
               supp_draw(x,y,q1_q3,n_l);
               supp_draw(x,y,q2_q3,n_l);
               hold off
            end
    end        
    
end

function res = supp_draw(x,y, q1_q2,l)
        for i=1:l
           x_1 = x + q1_q2(1,i);
           y_1 = y + q1_q2(2,i);
            
           hold on
           fill(x_1,y_1,'green');
           plot(x_1,y_1,'Color','green');
           hold off 
        end   
end

function res = X1_in_X0(x,q1,q2,q3,r)
    q1_q2 = supp_line(q1,q2);
    q1_q3 = supp_line(q1,q3);
    q2_q3 = supp_line(q2,q3);
    r_1 = supp_dist(x,q1_q2);
    r_2 = supp_dist(x,q1_q3);
    r_3 = supp_dist(x,q2_q3);
    r_1_b = r_1 <= r;
    r_2_b = r_2 <= r;
    r_3_b = r_3 <= r;
    m_b = max(max(max(r_1_b),max(r_2_b)),max(r_3_b));
   
    b = 0;
    num_1 = (q1(1) - x(1))*(q2(2) - q1(2)) - (q2(1) - q1(1))*(q1(2) - x(2));
    num_2 = (q2(1) - x(1))*(q3(2) - q2(2)) - (q3(1) - q2(1))*(q2(2) - x(2));
    num_3 = (q3(1) - x(1))*(q1(2) - q3(2)) - (q1(1) - q3(1))*(q3(2) - x(2));
    if (num_1 > 0) && (num_2 > 0) && (num_3 > 0)
        b = 1; 
    end
     if (num_1 < 0) && (num_2 < 0) && (num_3 < 0)
        b = 1; 
     end    
     
    if  b || m_b
        res = 1; 
    else
        res = 0;
    end    
      
end
