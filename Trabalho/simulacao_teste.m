% Definição dos parâmetros iniciais da integração
t = 2.2; 
T_sim = 1/100;
tempo = 0:T_sim:t;
%passo máximo ode
max_step = odeset('MaxStep', T_sim);
% condições iniciais de integração:
q1_0 = 0; 
q2_0 = 0;
theta_0 = 0;
q1p_0 = -3;
q2p_0 = -3; %Condição de velocidade inicial [m/s];
thetap_0 = 0;
%variáveis definidas também fora da ode
phi = 13*pi/180; 

x0 = [q1_0 q2_0 theta_0 q1p_0 q2p_0 thetap_0 0 300/3.6];
[t, y_L] = ode45(@fL, tempo, x0, max_step);
% [t, yl] = ode45(@fl, tempo, x0, max_step);

figure(1)
plot(tempo, y_L(:,1), "b")
% hold on
% plot(tempo, yl(:,1), "r")
% legend("Não-linear", "Linear")
title('Variação de q1')
xlabel('Tempo (s)')
ylabel('Posição (m)')

figure(2)
plot(tempo, y_L(:,2), "b")
% hold on
% plot(tempo, yl(:,2), "r")
% legend("Não-linear", "Linear")
title('Variação de q2')
xlabel('Tempo (s)')
ylabel('Posição (m)')

figure(3)
plot(tempo, y_L(:,3), "b")
% hold on
% plot(tempo, yl(:,3), "r")
% legend("Não-linear", "Linear")
title('Variação de theta em rad')
xlabel('Tempo (s)')
ylabel('Ângulo (rad)')

figure(4)
plot(tempo, 180*y_L(:,3)/pi, "b")
% hold on
% plot(tempo, 180*yl(:,1)/pi, "r")
% legend("Não-linear", "Linear")
title('Variação de theta em graus')
xlabel('Tempo (s)')
ylabel('Ângulo (rad)')

figure(5)
plot(tempo, y_L(:,4), "b")
% hold on
% plot(tempo, yl(:,4), "r")
% legend("Não-linear", "Linear")
title('Velocidade de q1')
xlabel('Tempo (s)')
ylabel('Velocidade (m/s)')

figure(6)
plot(tempo, y_L(:,5), "b")
% hold on
% plot(tempo, yl(:,5), "r")
% legend("Não-linear", "Linear")
title('Velocidade de q2')
xlabel('Tempo (s)')
ylabel('Velocidade (m/s)')

figure(7)
plot(tempo, y_L(:,6), "b")
% hold on
% plot(tempo, yl(:,6), "r")
% legend("Não-linear", "Linear")
title('Velocidade angular de theta')
xlabel('Tempo (s)')
ylabel('Velocidade (rad/s)')

figure(8)
plot(tempo, y_L(:,7), "b")
% hold on
% plot(tempo, yl(:,7), "r")
% legend("Não-linear", "Linear")
title("Deslocamento longitudinal")
xlabel('Tempo (s)')
ylabel('Posição (m)')

figure(9)
plot(tempo, y_L(:,8), "b")
% hold on
% plot(tempo, yl(:,8), "r")
% legend("Não-linear", "Linear")
title("Velocidade longitudinal")
xlabel('Tempo (s)')
ylabel('Velocidade (m/s)')

function dyLdt = fL(t, yL_0)
g = 9.8;
phi = 11*pi/180; 
rho = 1.2923; 
S = (27.66^2)/1.7; 
C_pav = 1; 
M = 88000; 
m = 2*3000;
C_L = C_pav*(-0.00165*(180*yL_0(3)/pi)^2+ 0.07378*180*yL_0(3)/pi + 0.21999);
C_D = C_pav*(0.00017*(180*yL_0(3)/pi)^2 + 0.01111*180*yL_0(3)/pi + 0.15714);
u_rol = (0.0041+0.000041*yL_0(2))*C_pav;
D_po = 5;
D_go = 2.2;
u_v = 0;
y_ext = 0;
yponto_ext = 0;
J_oz = 16864415;
k_r = 4*2913000;
c_r = 4*173580;
k_t = 2*4418000;
c_t = 2*255490;
dyLdt1 = yL_0(4);
dyLdt2 = yL_0(5);
dyLdt3 = yL_0(6);
dyLdt4 = -(((c_r + c_t)*yL_0(4))/m) + (c_t*yL_0(5))/m + (-(g*m) + k_t*(-yL_0(1) + yL_0(2)) + k_r*(-yL_0(1) + y_ext + c_r*yponto_ext))/m;
dyLdt5 = -0.5*(-2*g*M*J_oz + S*rho*C_L*J_oz*((yL_0(8) + u_v)*(yL_0(8) + u_v)) + 2*J_oz*k_t*(yL_0(1) - yL_0(2)) + M*cos(phi)*D_go*(D_go*(-((2*g*M - S*rho*C_L*((yL_0(8) + u_v)*(yL_0(8) + u_v)))*u_rol*(sin(phi) + cos(phi)*yL_0(3))) + 2*g*M*(cos(phi) - sin(phi)*yL_0(3))) - S*rho*D_po*((yL_0(8) + u_v)*(yL_0(8) + u_v))*(cos(phi)*(C_L + C_D*yL_0(3)) + sin(phi)*(C_D - C_L*yL_0(3)))))/(M*(M*(cos(phi)*cos(phi))*(D_go*D_go) - J_oz)) + (c_t*J_oz*yL_0(4))/(M*(-(M*(cos(phi)*cos(phi))*(D_go*D_go)) + J_oz)) + (c_t*J_oz*yL_0(5))/(M*M*(cos(phi)*cos(phi))*(D_go*D_go) - M*J_oz);
dyLdt6 = (sec(phi)*(-(S*rho*C_D*D_po*((yL_0(8) + u_v)*(yL_0(8) + u_v))*(tan(phi) + yL_0(3))) - 2*D_go*(k_t*(-yL_0(1) + yL_0(2)) + g*M*tan(phi)*yL_0(3) + g*M*u_rol*(tan(phi) + yL_0(3))) + S*rho*C_L*((yL_0(8) + u_v)*(yL_0(8) + u_v))*(D_po*(-1 + tan(phi)*yL_0(3)) + D_go*(1 + u_rol*(tan(phi) + yL_0(3))))))/(2*M*(D_go*D_go) - 2*(sec(phi)*sec(phi))*J_oz) + (sec(phi)*c_t*D_go*yL_0(4))/(M*(D_go*D_go) - sec(phi)*sec(phi)*J_oz) + (sec(phi)*c_t*D_go*yL_0(5))/(-(M*(D_go*D_go)) + sec(phi)*sec(phi)*J_oz);
dyLdt7 = yL_0(8);
dyLdt8 = (-(M*g - C_L*S*rho*yL_0(8)^2/2)*(u_rol)*C_pav - C_D*S*rho*yL_0(7)^2/2)/M;
dyLdt = [dyLdt1; dyLdt2; dyLdt3; dyLdt4; dyLdt5; dyLdt6; dyLdt7; dyLdt8];
end
    
