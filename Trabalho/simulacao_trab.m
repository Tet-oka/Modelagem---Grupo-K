
%% Definição dos parâmetros iniciais da integração
t = 5; 
T_sim = 1/100;
tempo = 0:T_sim:t;
%Passo máximo ODE
max_step = odeset('MaxStep', T_sim);

phi = 13*pi/180; 
g = 9.81;
rho = 1.2923; 
S = 96.37; 
C_pav = 1; 
M = 88000; 
m = 3000;
D_po = 5;
D_go = 2.2;
D_fo = 20.3;
J_oz = 16864415;
k_r = 13600000;
c_r = 9700;
k_t = 11486800;
c_t = 1021960;
m_f = 1000; 
k_tf = 5743400;
c_tf = 51098;
k_rf = 3400000;
c_rf = 2*2425;

%q1 e q2 definidos a partir da posição de equilíbrio estático do sistema em
%repouso

%Condições iniciais de integração:
q1_0 = (389300)/k_r; %q1 de equilíbrio dinâmico a 80m/s
q2_0 = (389300)/k_t + (389300)/k_r; %q2 de equilíbrio dinâmico a 80m/s
theta_0 = 0;
q3_0 = 0;
q1p_0 = -3;
q2p_0 = -3; %Condição de velocidade inicial [m/s];
thetap_0 = 0;
q3p_0 = -3;
pos_ini = 0;
vel_ini = 300/3.6;

x0 = [q1_0 q2_0 theta_0 q3_0 q1p_0 q2p_0 thetap_0 q3p_0 pos_ini vel_ini];

[t, y] = ode45(@f, tempo, x0, max_step);
[t, yL] = ode45(@fL, tempo, x0, max_step);


C_L = (-0.00165*((180*(y(:,3)+ phi)/pi).^2)) + (0.07378*180*(y(:,3)+ phi)/pi) + 0.21999;
C_D = (0.00017*((180*(y(:,3)+ phi)/pi).^2)) + (0.01111*180*(y(:,3)+ phi)/pi) + 0.15714;
L = C_L.*0.5.*S.*rho.*(y(:,10)).^2;



% figure(20)
% plot(tempo, C_L, "b")
% title("CL")
% xlabel('Tempo (s)')
% 
figure(21)
plot(tempo, L, "b")
title("L")
xlabel('Tempo (s)')

% figure(22)
% plot(tempo, C_D, "b")
% title("CD")
% xlabel('Tempo (s)')
% a_q1 = -(((c_r + c_t).*y(:,4))/m) + (c_t.*y(:,5))/m + (-(g*m) + k_t.*(-y(:,1) + y(:,2)) + k_r.*(-y(:,1) + y_ext + c_r*yponto_ext))/m;
% a_q2 = -0.5*(-2*g*M*J_oz + S*rho*C_L*J_oz.*((y(:,8) + u_v).*(y(:,8) + u_v)) + 2*J_oz*k_t.*(y(:,1) - y(:,2)) + M*cos(phi)*D_go*(D_go*(-((2*g*M - S*rho*C_L.*((y(:,8) + u_v).*(y(:,8) + u_v)))*u_rol*(sin(phi) + cos(phi).*y(:,3))) + 2*g*M*(cos(phi) - sin(phi).*y(:,3))) - S*rho*D_po.*((y(:,8) + u_v).*(y(:,8) + u_v))*(cos(phi)*(C_L + C_D.*y(:,3)) + sin(phi)*(C_D - C_L.*y(:,3)))))/(M*(M*(cos(phi)*cos(phi))*(D_go*D_go) - J_oz)) + (c_t*J_oz.*y(:,4))/(M*(-(M*(cos(phi)*cos(phi))*(D_go*D_go)) + J_oz)) + (c_t*J_oz.*y(:,5))/(M*M*(cos(phi)*cos(phi))*(D_go*D_go) - M*J_oz);

%% Plot dos gráficos

figure(1)
plot(tempo, y(:,1), "b")
hold on
plot(tempo, yL(:,1), "r")
legend("Não-linear", "Linear")
title('Variação de q1')
xlabel('Tempo (s)')
ylabel('Posição (m)')

figure(2)
plot(tempo, y(:,2), "b")
hold on
plot(tempo, yL(:,2), "r")
legend("Não-linear", "Linear")
title('Variação de q2')
xlabel('Tempo (s)')
ylabel('Posição (m)')

figure(3)
plot(tempo, y(:,3), "b")
hold on
plot(tempo, yL(:,3), "r")
legend("Não-linear", "Linear")
title('Variação de theta em rad')
xlabel('Tempo (s)')
ylabel('Ângulo (rad)')

figure(4)
plot(tempo, 180*y(:,3)/pi, "b")
hold on
plot(tempo, 180*yL(:,3)/pi, "r")
legend("Não-linear", "Linear")
title('Variação de theta em graus para Phi=13')
xlabel('Tempo (s)')
ylabel('Ângulo (graus)')

figure(5)
plot(tempo, y(:,5), "b")
hold on
plot(tempo, yL(:,5), "r")
legend("Não-linear", "Linear")
title('Velocidade de q1')
xlabel('Tempo (s)')
ylabel('Velocidade (m/s)')

figure(6)
plot(tempo, y(:,6), "b")
hold on
plot(tempo, yL(:,6), "r")
legend("Não-linear", "Linear")
title('Velocidade de q2')
xlabel('Tempo (s)')
ylabel('Velocidade (m/s)')

figure(7)
plot(tempo, y(:,7), "b")
hold on
plot(tempo, yL(:,7), "r")
legend("Não-linear", "Linear")
title('Velocidade angular de theta')
xlabel('Tempo (s)')
ylabel('Velocidade (rad/s)')

figure(8)
plot(tempo, y(:,9), "b")
hold on
plot(tempo, yL(:,9), "r")
legend("Não-linear", "Linear")
title("Deslocamento longitudinal")
xlabel('Tempo (s)')
ylabel('Posição (m)')

figure(9)
plot(tempo, y(:,10), "b")
hold on
plot(tempo, yL(:,10), "r")
legend("Não-linear", "Linear")
title("Velocidade longitudinal")
xlabel('Tempo (s)')
ylabel('Velocidade (m/s)')

figure(10)
plot(tempo, (180/pi)*abs(y(:,3)-yL(:,3)), "b")
title("Diferença entre modelo linear e não linear para theta")
xlabel('Tempo (s)')
ylabel('Diferença (graus)')

figure(11)
plot(tempo, abs(y(:,2)-yL(:,2)), "b")
title("Diferença entre modelo linear e não linear para q2")
xlabel('Tempo (s)')
ylabel('Diferença (m)')

figure(12)
plot(tempo, abs(y(:,1)-yL(:,1)), "b")
title("Diferença entre modelo linear e não linear para q1")
xlabel('Tempo (s)')
ylabel('Diferença (m)')

% figure(11)
% plot(tempo, y(:,1), "r")
% hold on
% plot(tempo, y(:,2), "g")
% xlabel('Tempo (s)')
% ylabel('Posição (m)')
% legend('Variação de q1','Variação de q2')

%% Funções ODE

function dydt = f(t, y_0)
phi = 13*pi/180; 
g = 9.81;
rho = 1.2923; 
S = 96.37; 
C_pav = 1; 
M = 88000; 
m = 3000;
C_L = ((-0.00165*((180*(y_0(3)+ phi)/pi).^2)) + (0.07378*180*(y_0(3)+ phi)/pi) + 0.21999);
C_D = ((0.00017*((180*(y_0(3)+ phi)/pi).^2)) + (0.01111*180*(y_0(3)+ phi)/pi) + 0.15714);
u_rol = (0.0041+0.000041*y_0(10))*C_pav;
D_po = 5;
D_go = 2.2;
D_fo = 20.3;
u_v = 0;
y_ext = 0;
yponto_ext = 0;
J_oz = 16864415;
k_r = 13600000;
c_r = 9700;
k_t = 11486800;
c_t = 1021960;
m_f = 1000; 
k_tf = 5743400;
c_tf = 51098;
k_rf = 3400000;
c_rf = 2*2425;

%Os dados que mudam com o chaveamento devem ser inseridos dentro do if/else
if -y_0(3) < phi
    dydt1 = y_0(5); %q1' SEM MG a partir de 0
    dydt2 = y_0(6); %q2'
    dydt3 = y_0(7); %theta'
    dydt4 = y_0(8); %q3'
    dydt5 = -(((c_r + c_t)*y_0(5))/m) + (c_t*y_0(6))/m + (k_t*(-y_0(1) + y_0(2)) + k_r*(-y_0(1) + y_ext) + c_r*yponto_ext)/m;
    dydt6 = -0.5*(S*rho*C_L*J_oz*((y_0(10) + u_v)*(y_0(10) + u_v)) + M*cos(phi + y_0(3))*D_go*(-(S*rho*(sin(phi + y_0(3))*C_D + cos(phi + y_0(3))*C_L)*D_po*((y_0(10) + u_v)*(y_0(10) + u_v))) + D_go*(2*g*M*cos(phi + y_0(3)) + sin(phi + y_0(3))*(-2*g*M + S*rho*C_L*((y_0(10) + u_v)*(y_0(10) + u_v)))*u_rol)) + 2*J_oz*k_t*(y_0(1) - y_0(2)))/(M*(M*(cos(phi + y_0(3))*cos(phi + y_0(3)))*(D_go*D_go) - J_oz)) + (c_t*J_oz*y_0(5))/(M*(-(M*(cos(phi + y_0(3))*cos(phi + y_0(3)))*(D_go*D_go)) + J_oz)) + (c_t*J_oz*y_0(6))/(M*M*(cos(phi + y_0(3))*cos(phi + y_0(3)))*(D_go*D_go) - M*J_oz) - (sin(phi + y_0(3))*D_go*J_oz*(y_0(7)*y_0(7)))/(M*(cos(phi + y_0(3))*cos(phi + y_0(3)))*(D_go*D_go) - J_oz);
    dydt7 = -((sec(phi + y_0(3))*(S*rho*D_po*((y_0(10) + u_v)*(y_0(10) + u_v))*(C_L + C_D*tan(phi + y_0(3))) + D_go*(-2*g*M + 2*g*M*u_rol*tan(phi + y_0(3)) - S*rho*C_L*((y_0(10) + u_v)*(y_0(10) + u_v))*(1 + u_rol*tan(phi + y_0(3))) + 2*k_t*(-y_0(1) + y_0(2)))))/(2*M*(D_go*D_go) - 2*(sec(phi + y_0(3))*sec(phi + y_0(3)))*J_oz)) + (sec(phi + y_0(3))*c_t*D_go*y_0(5))/(M*(D_go*D_go) - sec(phi + y_0(3))*sec(phi + y_0(3))*J_oz) + (sec(phi + y_0(3))*c_t*D_go*y_0(6))/(-(M*(D_go*D_go)) + sec(phi + y_0(3))*sec(phi + y_0(3))*J_oz) + (M*(D_go*D_go)*tan(phi + y_0(3))*(y_0(7)*y_0(7)))/(M*(D_go*D_go) - sec(phi + y_0(3))*sec(phi + y_0(3))*J_oz);
    dydt8 = 0;
    dydt9 = y_0(10);
    dydt10 = (-((M+m)*g - C_L*S*rho*((y_0(10)+u_v)^2)/2)*(u_rol) - C_D*S*rho*((y_0(10)+u_v)^2)/2)/(M+m);
    dydt = [dydt1; dydt2; dydt3; dydt4; dydt5; dydt6; dydt7; dydt8; dydt9; dydt10];
else
    dydt1 = y_0(5);
    dydt2 = y_0(6);
    dydt3 = y_0(7);
    dydt4 = y_0(8);
    dydt5 = -(((c_r + c_t)*y_0(5))/m) + (c_t*y_0(6))/m + (k_t*(-y_0(1) + y_0(2)) + k_r*(-y_0(1) + y_ext) + c_r*yponto_ext)/m;
    dydt6 = (c_t*J_oz*y_0(5))/(M*(-(M*(cos(y_0(3))*cos(y_0(3)))*(D_go*D_go)) + J_oz)) + ((-(M*(cos(y_0(3))*cos(y_0(3)))*c_tf*D_fo*D_go) + (c_t + c_tf)*J_oz)*y_0(6))/(M*(M*(cos(y_0(3))*cos(y_0(3)))*(D_go*D_go) - J_oz)) + (cos(y_0(3))*c_tf*D_fo*(-(M*(cos(y_0(3))*cos(y_0(3)))*D_fo*D_go) + J_oz)*y_0(7))/(M*(M*(cos(y_0(3))*cos(y_0(3)))*(D_go*D_go) - J_oz)) + (D_go*J_oz*(y_0(7)*y_0(7)))/(-(M*cos(y_0(3))*cot(y_0(3))*(D_go*D_go)) + csc(y_0(3))*J_oz) - (M*cos(y_0(3))*(D_go*D_go)*(2*g*M*cos(y_0(3)) + sin(y_0(3))*(-2*g*M + S*rho*C_L*((y_0(10) + u_v)*(y_0(10) + u_v)))*u_rol) + J_oz*(S*rho*C_L*((y_0(10) + u_v)*(y_0(10) + u_v)) + 2*k_t*(y_0(1) - y_0(2)) - 2*k_tf*(sin(y_0(3))*D_fo + y_0(2) - y_0(4))) + M*cos(y_0(3))*D_go*(-(S*rho*(sin(y_0(3))*C_D + cos(y_0(3))*C_L)*D_po*((y_0(10) + u_v)*(y_0(10) + u_v))) + 2*cos(y_0(3))*D_fo*k_tf*(sin(y_0(3))*D_fo + y_0(2) - y_0(4))))/(2.*M*(M*(cos(y_0(3))*cos(y_0(3)))*(D_go*D_go) - J_oz)) + (c_tf*(M*(cos(y_0(3))*cos(y_0(3)))*D_fo*D_go - J_oz)*y_0(8))/(M*(M*(cos(y_0(3))*cos(y_0(3)))*(D_go*D_go) - J_oz));
    dydt7 = (sec(y_0(3))*c_t*D_go*y_0(5))/(M*(D_go*D_go) - sec(y_0(3))*sec(y_0(3))*J_oz) + (sec(y_0(3))*(-(c_tf*D_fo) + (c_t + c_tf)*D_go)*y_0(6))/(-(M*(D_go*D_go)) + sec(y_0(3))*sec(y_0(3))*J_oz) + (c_tf*D_fo*(D_fo - D_go)*y_0(7))/(M*(D_go*D_go) - sec(y_0(3))*sec(y_0(3))*J_oz) + (y_0(7)*y_0(7))/(cot(y_0(3)) - (2*csc(2*y_0(3))*J_oz)/(M*(D_go*D_go))) - (sec(y_0(3))*(S*rho*D_po*((y_0(10) + u_v)*(y_0(10) + u_v))*(C_L + C_D*tan(y_0(3))) + D_go*(-2*g*M + 2*g*M*u_rol*tan(y_0(3)) - S*rho*C_L*((y_0(10) + u_v)*(y_0(10) + u_v))*(1 + u_rol*tan(y_0(3))) + 2*k_t*(-y_0(1) + y_0(2)) + 2*k_tf*(sin(y_0(3))*D_fo + y_0(2) - y_0(4))) - 2*D_fo*k_tf*(sin(y_0(3))*D_fo + y_0(2) - y_0(4))))/(2*M*(D_go*D_go) - 2*(sec(y_0(3))*sec(y_0(3)))*J_oz) + (sec(y_0(3))*c_tf*(D_fo - D_go)*y_0(8))/(-(M*(D_go*D_go)) + sec(y_0(3))*sec(y_0(3))*J_oz);
    dydt8 = (c_tf*y_0(6))/m_f + (cos(y_0(3))*c_tf*D_fo*y_0(7))/m_f - ((c_rf + c_tf)*y_0(8))/m_f + (k_tf*(sin(y_0(3))*D_fo + y_0(2) - y_0(4)) + k_rf*(-y_0(4) + y_ext) + c_rf*yponto_ext)/m_f;
    dydt9 = y_0(10);
    dydt10 = (-((M+m)*g - C_L*S*rho*((y_0(10)+u_v)^2)/2)*(u_rol) - C_D*S*rho*((y_0(10)+u_v)^2)/2)/(M+m);
    dydt = [dydt1; dydt2; dydt3; dydt4; dydt5; dydt6; dydt7; dydt8; dydt9; dydt10];
end
end
 

function dyLdt = fL(t, yL_0)
phi = 13*pi/180; 
g = 9.81;
rho = 1.2923; 
S = 96.37; 
C_pav = 1; 
M = 88000; 
m = 3000;
C_L = ((-0.00165*((180*(yL_0(3)+ phi)/pi).^2)) + (0.07378*180*(yL_0(3)+ phi)/pi) + 0.21999);
C_D = ((0.00017*((180*(yL_0(3)+ phi)/pi).^2)) + (0.01111*180*(yL_0(3)+ phi)/pi) + 0.15714);
u_rol = (0.0041+0.000041*yL_0(10))*C_pav;
D_po = 5;
D_go = 2.2;
D_fo = 20.3;
u_v = 0;
y_ext = 0;
yponto_ext = 0;
J_oz = 16864415;
k_r = 13600000;
c_r = 9700;
k_t = 11486800;
c_t = 1021960;
m_f = 1000; 
k_tf = 5743400;
c_tf = 51098;
k_rf = 3400000;
c_rf = 2*2425;
vel_ini = 300/3.6;
% Os dados que mudam com o chaveamento devem ser inseridos dentro do if/else
if -yL_0(3) < phi
    dyLdt1 = yL_0(5);
    dyLdt2 = yL_0(6);
    dyLdt3 = yL_0(7);
    dyLdt4 = yL_0(8);
    dyLdt5 = -(((c_r + c_t)*yL_0(5))/m) + (c_t*yL_0(6))/m + (-((S*rho*C_L*(k_r - k_t)*((yL_0(10) + u_v)*(yL_0(10) + u_v)))/k_r) + 2*(k_t*(-yL_0(1) + yL_0(2)) + k_r*(-yL_0(1) + y_ext) + c_r*yponto_ext))/(2.*m);
   %dyLdt5 = -(((c_r + c_t)*yL_0(5))/m) + (c_t*yL_0(6))/m + (k_t*(-yL_0(1) + yL_0(2)) + k_r*(-yL_0(1) + y_ext) + c_r*yponto_ext)/m;
    dyLdt6 = -0.5*(J_oz*((S*rho*C_L*(k_r - k_t)*((yL_0(10) + u_v)*(yL_0(10) + u_v)))/k_r + 2*k_t*(yL_0(1) - yL_0(2))) + M*cos(phi)*(D_go*D_go)*(-((2*g*M - S*rho*C_L*((yL_0(10) + u_v)*(yL_0(10) + u_v)))*u_rol*(sin(phi) + cos(phi)*yL_0(3))) + 2*g*M*(cos(phi) - sin(phi)*yL_0(3))) - M*S*rho*cos(phi)*D_go*D_po*((yL_0(10) + u_v)*(yL_0(10) + u_v))*(cos(phi)*(C_L + C_D*yL_0(3)) + sin(phi)*(C_D - C_L*yL_0(3))))/(M*(M*(cos(phi)*cos(phi))*(D_go*D_go) - J_oz)) + (c_t*J_oz*yL_0(5))/(M*(-(M*(cos(phi)*cos(phi))*(D_go*D_go)) + J_oz)) + (c_t*J_oz*yL_0(6))/(M*M*(cos(phi)*cos(phi))*(D_go*D_go) - M*J_oz);
    dyLdt7 = -0.5*(sec(phi)*(S*rho*D_po*k_r*((yL_0(10) + u_v)*(yL_0(10) + u_v))*(C_D*(tan(phi) + yL_0(3)) + C_L*(1 - tan(phi)*yL_0(3))) + D_go*(S*rho*C_L*k_t*((yL_0(10) + u_v)*(yL_0(10) + u_v)) + k_r*(-(S*rho*C_L*((yL_0(10) + u_v)*(yL_0(10) + u_v))*(1 + u_rol*(tan(phi) + yL_0(3)))) + 2*(-(g*M) + k_t*(-yL_0(1) + yL_0(2)) + g*M*tan(phi)*yL_0(3) + g*M*u_rol*(tan(phi) + yL_0(3)))))))/((M*(D_go*D_go) - sec(phi)*sec(phi)*J_oz)*k_r) + (sec(phi)*c_t*D_go*yL_0(5))/(M*(D_go*D_go) - sec(phi)*sec(phi)*J_oz) + (sec(phi)*c_t*D_go*yL_0(6))/(-(M*(D_go*D_go)) + sec(phi)*sec(phi)*J_oz);
    dyLdt8 = 0;
    dyLdt9 = yL_0(10);
    dyLdt10 = (-((M+m)*g - C_L*S*rho*((yL_0(10)+u_v)^2)/2)*(u_rol) - C_D*S*rho*((yL_0(10)+u_v)^2)/2)/(M+m);
    dyLdt = [dyLdt1; dyLdt2; dyLdt3; dyLdt4; dyLdt5; dyLdt6; dyLdt7; dyLdt8; dyLdt9; dyLdt10];
else
    dyLdt1 = 1;
    dyLdt2 = 1;
    dyLdt3 = 1; 
    dyLdt4 = 1;
    dyLdt5 = 1;
    dyLdt6 = 1;
    dyLdt7 = 1;
    dyLdt8 = 1;
    dyLdt9 = yL_0(10);
    dyLdt10 = (-((M+m)*g - C_L*S*rho*((yL_0(10)+u_v)^2)/2)*(u_rol) - C_D*S*rho*((yL_0(10)+u_v)^2)/2)/(M+m);
    dyLdt = [dyLdt1; dyLdt2; dyLdt3; dyLdt4; dyLdt5; dyLdt6; dyLdt7; dyLdt8; dyLdt9; dyLdt10];
end
end