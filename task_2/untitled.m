%%
%ПЕРВЫЙ
%СЛУЧАЙ--------------------------------------------Ищем начальные условия
%для psi через параметризацию 


alpha = 2;
T = 3;
k_1 = 0;%k_1^2 + 2*k_2 - 4g + 1 > 0
k_2 = 3;
L = 0.5;
eps = 1;
g = 0.4;
p = k_1^2 + 2*k_2 - 4*g + 1;


psi_0 = -1./2;
n = 100;
psi_0_alpha = linspace(0.1,alpha,n);
t_0_T = linspace(0,T,n);

u_1m = [];
u_2m = [];
psi_1m = [];
psi_2m = [];
t_m = [];
x_1m = [];
x_2m = [];
sw_t = [];
time_mas = [];

u_1_opt = 1;
u_2_opt = 1;
min_u = 1000;
for i = 1:n
    for j = 1:n
    
x_1 = L;
x_2 = 0;
psi_2 = psi_0_alpha(i);
tq = t_0_T(j);
psi_1 = (1/(2*sinh(tq*p/2)))*(((1 + k_1)*sinh(tq*p/2) + p*cosh(tq*p/2))*psi_2 - alpha*p*exp(-tq*(1 + k_1)/2));

%определение оптимальной компоненты u_1
if -psi_2/(2*psi_0) < 0
    u_1_opt = 0;
end    

if (-psi_2/(2*psi_0)>= 0) && (-psi_2/(2*psi_0) <= alpha)
    u_1_opt = -psi_2./(2*psi_0);
end

if -psi_2/(2*psi_0)> alpha
    u_1_opt = alpha;
end

%определение оптимальной компоненты u_2
if psi_2 * x_2 > 0
    u_2_opt = k_1;
end    

if  psi_2 * x_2 < 0
    u_2_opt = k_2;
end



b = 1;
options = odeset('Events',@eventfcn,'Events',@eventfcn_1);
if (u_1_opt < alpha) && (u_1_opt > 0)
   b_f_dif = 1; 
else
   b_f_dif = 0; 
end    
t1 = 0;
k = 1;
u_mas = [];
t_mas = [];
psi_1_mas = [];
psi_2_mas = [];
x_1_mas = [];
x_2_mas = [];
ur = [];
t_mas_sw = [];
while b
    tspan = [t1 T];
    if b_f_dif
        [t_f, x_f] = ode45(@(t_f,x_f)odefcn_main_1(t_f,x_f,g,u_2_opt),tspan,[x_1 x_2 psi_1 psi_2],options);
        
        t_mas = [t_mas, t_f'];
        u_mas = [u_mas, x_f(:,4)'];
        ur = [ur, u_2_opt*ones(1,numel(t_f))];
        x_1_mas = [x_1_mas, x_f(:,1)'];
        x_2_mas = [x_2_mas, x_f(:,2)'];
        psi_1_mas = [psi_1_mas, x_f(:,3)'];
        psi_2_mas = [psi_2_mas, x_f(:,4)'];
        t_mas_sw = [t_mas_sw, t_f(end)];
    else
        [t_f, x_f] = ode45(@(t_f,x_f)odefcn_main(t_f,x_f,g,u_1_opt,u_2_opt),tspan,[x_1 x_2 psi_1 psi_2],options);
        
        t_mas = [t_mas, t_f'];
        u_1_opt_mas = u_1_opt.*ones(1,numel(t_f));
        ur = [ur, u_2_opt*ones(1,numel(t_f))];
        u_mas = [u_mas, u_1_opt_mas];
        x_1_mas = [x_1_mas, x_f(:,1)'];
        x_2_mas = [x_2_mas, x_f(:,2)'];
        psi_1_mas = [psi_1_mas, x_f(:,3)'];
        psi_2_mas = [psi_2_mas, x_f(:,4)'];
        t_mas_sw = [t_mas_sw, t_f(end)];
    end
   
    k = k + 1;
    if t_f(end) == T
        b = 0;
       
    else    
        
        x_1 = x_f(end,1);
        x_2 = x_f(end,2);
        psi_1 = x_f(end,3);
        psi_2 = x_f(end,4);
        t1 = t_f(end);
        [out1, out2, out3] = check_pos(x_f(:,2)',x_f(:,4)',alpha);
        if out3 == 1
            u_2_opt = k_2;
        else
            u_2_opt = k_1;
        end    
        if out1 == 1
            b_f_dif = 1;
        else
            b_f_dif = 0;
            u_1_opt = out2;
        end    
    end
    
    if k == 100
        b = 0;
       
    end    
end
   
    if abs(x_1_mas(end)) + abs(x_2_mas(end)) <= eps
        if trapz(t_mas,u_mas.^2) < min_u
         
            min_u = trapz(t_mas,u_mas.^2);
            time_mas = t_mas;
            u_1m = u_mas;
            u_2m = ur;
            
            psi_1m = psi_1_mas;
            psi_2m = psi_2_mas;
            
            x_1m = x_1_mas;
            x_2m = x_2_mas;
            
            sw_t = t_mas_sw;
    
        end

    end    
    end    
end 
disp('J(u)');
disp(min_u);
figure;
plot(time_mas,x_1m);
xlabel('t');
ylabel('x_1');

figure;
plot(time_mas,x_2m);
xlabel('t');
ylabel('x_2');

figure;
plot(time_mas,psi_1m);
xlabel('t');
ylabel('psi_1');

figure;
plot(time_mas,psi_2m);
xlabel('t');
ylabel('psi_2');

figure;
plot(time_mas,u_1m);
xlabel('t');
ylabel('u_1');

figure;
plot(time_mas,u_2m);
xlabel('t');
ylabel('u_2');

figure;
plot(x_1m,x_2m);
xlabel('x_1');
ylabel('x_2');

%%
%ВТОРОЙ
%СЛУЧАЙ-----------------------------------------------------------------------------------------(Отличие только в том, что мы ищем по x_2 \in [0,eps^'])
alpha = 2;
T = 3;
k_1 = 0;%k_1^2 + 2*k_2 - 4g + 1 > 0
k_2 = 3;
L = 0.5;
eps = 1;
g = 0.4;
p = k_1^2 + 2*k_2 - 4*g + 1;

psi_0 = -1./2;

n = 100;
min_u = 1000;
x_2_end_mas = linspace(1,eps,n);%массив конечных значенией по x2

for i=1 : n
   a = cosh((T*p)/2);
   b = sinh((T*p)/2); 
   c = exp(T*(1 + k_1));
   psi_2 = x_2_end_mas(i)*2*(1 + k_1)*p/(exp(-T*(1+k_1)/2)*((c-1)*p*a - (c - 1)*(p*a)+(1+k_1)*b - (1 + k_1)*(c+1)*b));
   psi_1 = (p*a + (1 +k_1)*b)*psi_2/(2*b);
   
   x_1 = L;
   x_2 = 0;
   t1 = 0;
   %определение оптимальной компоненты u_2
    if psi_2 * x_2 >= 0
        u_2_opt = k_1;
    end    

    if  psi_2 * x_2 < 0
        u_2_opt = k_2;        
    end
    
    options = odeset('Events',@eventfcn);
    u_mas = [];
    t_mas = [];
    psi_1_mas = [];
    psi_2_mas = [];
    x_1_mas = [];
    x_2_mas = [];
    ur = [];
    t_mas_sw = [];
    while b
        tspan = [t1 T];
    
        [t_f, x_f] = ode45(@(t_f,x_f)odefcn_main_1(t_f,x_f,g,u_2_opt),tspan,[x_1 x_2 psi_1 psi_2],options);
        t_mas = [t_mas, t_f'];
        u_mas = [u_mas, x_f(:,4)'];
        ur = [ur, u_2_opt*ones(1,numel(t_f))];
        x_1_mas = [x_1_mas, x_f(:,1)'];
        x_2_mas = [x_2_mas, x_f(:,2)'];
        psi_1_mas = [psi_1_mas, x_f(:,3)'];
        psi_2_mas = [psi_2_mas, x_f(:,4)'];
        t_mas_sw = [t_mas_sw, t_f(end)];
   
    k = k + 1;
    if t_f(end) == T
        b = 0;
       
    else    
        
        x_1 = x_f(end,1);
        x_2 = x_f(end,2);
        psi_1 = x_f(end,3);
        psi_2 = x_f(end,4);
        t1 = t_f(end);
        [out1, out2, out3] = check_pos(x_f(:,2),x_f(:,4),alpha);
        if out3 == 1
            u_2_opt = k_2;
        else
            u_2_opt = k_1;
        end    
        if out1 == 1
            b_f_dif = 1;
        else
            b_f_dif = 0;
            u_1_opt = out2;
        end    
    end
    
    if k == 100
        b = 0;
       
    end    
    end
    if (abs(x_1_mas(end)) + abs(x_2_mas(end)) <= eps) %Это условие так же проверяет принадлежит же наше предположение x_2(end) \in [0,\upvarepsilon^'] где \upvarepsilon^' = \upvarepsilon - abs(x_1(T))
        if trapz(t_mas,u_mas.^2) < min_u
         
            min_u = trapz(t_mas,u_mas.^2);
            time_mas = t_mas;
            u_1m = u_mas;
            u_2m = ur;
            
            psi_1m = psi_1_mas;
            psi_2m = psi_2_mas;
            
            x_1m = x_1_mas;
            x_2m = x_2_mas;
            
            sw_t = t_mas_sw;
    
        end

    end    
end    
%%
disp('J(u)');
disp(min_u);
figure;
plot(time_mas,x_1m);
xlabel('t');
ylabel('x_1');

figure;
plot(time_mas,x_2m);
xlabel('t');
ylabel('x_2');

figure;
plot(time_mas,psi_1m);
xlabel('t');
ylabel('psi_1');

figure;
plot(time_mas,psi_2m);
xlabel('t');
ylabel('psi_2');

figure;
plot(time_mas,u_1m);
xlabel('t');
ylabel('u_1');

figure;
plot(time_mas,u_2m);
xlabel('t');
ylabel('u_2');

figure;
plot(x_1m,x_2m);
xlabel('x_1');
ylabel('x_2');
%%
%Сопряженная система.

function dydt = odefcn_main(t,y,g,u_1,u_2)
    dydt = zeros(4,1);
    dydt(1) = y(2); %x
    dydt(2) = u_1 - y(2)*(1 + u_2) - g*y(1);%x
    
    dydt(3) = g*y(4);%psi
    dydt(4) = -y(3) + y(4)*(1 + u_2);%psi
end

function dydt = odefcn_main_1(t,y,g,u_2)
    dydt = zeros(4,1);
    dydt(1) = y(2); %x
    dydt(2) = y(4) - y(2)*(1 + u_2) - g*y(1);%x
    
    dydt(3) = g*y(4);%psi
    dydt(4) = -y(3) + y(4)*(1 + u_2);%psi
end

function [value, isterminal, direction] = eventfcn(t,x)
    value = x(2);
    isterminal = 1;
    direction = 0;
end

function [value, isterminal, direction] = eventfcn_1(t,x)
    value = (x(4) - 2)*x(4);
    isterminal = 1;
    direction = 0;
end

function [out_1,out_2,out_3] = check_pos(x_2,psi_2,alpha)
        eps1 = 0.01;
        out_1 = 0;
        out_2 = 0;
        out_3 = 0;
        if psi_2(end - 3) < 0
            out_1 = 1;
            out_2 = 0;
        else
            
            if psi_2(end - 3) > alpha
                out_1 = 1;
                out_2 = 0;
            else 
                if abs(psi_2(end) - alpha) < eps1
                out_1 = 0;
                out_2 = alpha;
                else
                    
                    if abs(psi_2(end)) < eps1
                    out_1 = 0;
                    out_2 = 0;
                    end   
                end    
            end    
        end    
        
        if psi_2(end) * x_2(end) < 0
            out_3 = 1;
        else
            out_3 = 0;
        end 
            
end