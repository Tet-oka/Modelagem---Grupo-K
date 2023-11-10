% Definição dos parâmetros iniciais da integração
t = 10; 
T_sim = 1/100;
tempo = 0:T_sim:t;
%passo máximo ode
max_step = odeset('MaxStep', T_sim);
% condições iniciais de integração:
q1_0 = 0; 
q2_0 = 0;
theta_0 = 0;
q1p_0 = -300/3.6;
q2p_0 = -300/3.6; %Condição de velocidade inicial [m/s];
thetap_0 = 0;


x0 = [q1_0 q2_0 theta_0 q1p_0 q2p_0 thetap_0];
[t, y] = ode45(@f, tempo, x0, max_step);

figure(1)
plot(tempo, y(:,1))
title("q1")

figure(2)
plot(tempo, y(:,2))
title("q2")

figure(3)
plot(tempo, y(:,3))
title("theta")

figure(4)
plot(tempo, y(:,4))
title("q1 ponto")

figure(5)
plot(tempo, y(:,5))
title("q2 ponto")

figure(6)
plot(tempo, y(:,6))
title("theta ponto")
  

function dydt = f(t, y_0)
g = 9.8;
phi = 0; 
rho = 1.2923; 
S = (40^2)/1.7; 
C_pav = 1; 
M = 88000; 
m = 2*3000;
C_L = C_pav*(-0.00165*(180*y_0(3)/pi)^2+ 0.07378*180*y_0(3)/pi + 0.21999);
C_D = C_pav*(0.00017*(180*y_0(3)/pi)^2 + 0.01111*180*y_0(3)/pi + 0.15714);
u_rol = (0.0041+0.000041*y_0(2))*C_pav;
D_po = 5;
D_go = 2.2;
u_long = 60;
u_v = 0;
y_ext = 0;
yponto_ext = 0;
J_oz = 16864415;
k_r = 4*2913000;
c_r = 4*17358;
k_t = 2*4418000;
c_t = 2*25549;
dydt1 = y_0(4);
dydt2 = y_0(5);
dydt3 = y_0(6);
dydt4 = -(((c_r + c_t)*y_0(4))/m) + (c_t*y_0(5))/m + (-(g*m) + k_t*(-y_0(1) + y_0(2)) + k_r*(-y_0(1) + y_ext) + c_r*yponto_ext)/m;
dydt5 = -0.5*(-2*g*M*J_oz + S*rho*C_L*J_oz*((u_long + u_v)*(u_long + u_v)) + M*cos(phi + y_0(3))*D_go*(-(S*rho*(sin(phi + y_0(3))*C_D + cos(phi + y_0(3))*C_L)*D_po*((u_long + u_v)*(u_long + u_v))) + D_go*(2*g*M*cos(phi + y_0(3)) + sin(phi + y_0(3))*(-2*g*M + S*rho*C_L*((u_long + u_v)*(u_long + u_v)))*u_rol)) + 2*J_oz*k_t*(y_0(1) - y_0(2)))/(M*(M*(cos(phi + y_0(3))*cos(phi + y_0(3)))*(D_go*D_go) - J_oz)) + (c_t*J_oz*y_0(4))/(M*(-(M*(cos(phi + y_0(3))*cos(phi + y_0(3)))*(D_go*D_go)) + J_oz)) + (c_t*J_oz*y_0(5))/(M*M*(cos(phi + y_0(3))*cos(phi + y_0(3)))*(D_go*D_go) - M*J_oz) - (sin(phi + y_0(3))*D_go*J_oz*(y_0(6)*y_0(6)))/(M*(cos(phi + y_0(3))*cos(phi + y_0(3)))*(D_go*D_go) - J_oz);
dydt6 = (sec(phi + y_0(3))*(-(S*rho*C_D*D_po*((u_long + u_v)*(u_long + u_v))*tan(phi + y_0(3))) + S*rho*C_L*((u_long + u_v)*(u_long + u_v))*(-D_po + D_go*(1 + u_rol*tan(phi + y_0(3)))) - 2*D_go*(g*M*u_rol*tan(phi + y_0(3)) + k_t*(-y_0(1) + y_0(2)))))/(2*M*(D_go*D_go) - 2*(sec(phi + y_0(3))*sec(phi + y_0(3)))*J_oz) + (sec(phi + y_0(3))*c_t*D_go*y_0(4))/(M*(D_go*D_go) - sec(phi + y_0(3))*sec(phi + y_0(3))*J_oz) + (sec(phi + y_0(3))*c_t*D_go*y_0(5))/(-(M*(D_go*D_go)) + sec(phi + y_0(3))*sec(phi + y_0(3))*J_oz) + (M*(D_go*D_go)*tan(phi + y_0(3))*(y_0(6)*y_0(6)))/(M*(D_go*D_go) - sec(phi + y_0(3))*sec(phi + y_0(3))*J_oz);
dydt = [dydt1; dydt2; dydt3; dydt4; dydt5; dydt6];
end
