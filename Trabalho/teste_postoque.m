t = 100; 
T_sim = 1/100;
tempo = 0:T_sim:t;
%passo máximo ode
max_step = odeset('MaxStep', T_sim);
% condições iniciais de integração:
q1_0 = 2; 
q2_0 =0;
theta_0 = 0;
q3_0 = 0;
q1p_0 = 0;
q2p_0 = 0; %Condição de velocidade inicial [m/s];
thetap_0 = 0;
q3p_0 = 0;
pos_ini = 0;
vel_ini = 80;



x0 = [q1_0 q2_0 theta_0 q3_0 q1p_0 q2p_0 thetap_0 q3p_0 pos_ini vel_ini];
[t, y] = ode45(@ft, tempo, x0, max_step);

figure(1)
plot(tempo, y(:,1), "b")
% hold on
% plot(tempo, yL(:,1), "r")
%legend("Não-linear", "Linear")
title('Variação de q1')
xlabel('Tempo (s)')
ylabel('Posição (m)')

figure(2)
plot(tempo, y(:,2), "b")
% hold on
% plot(tempo, yL(:,2), "r")
%legend("Não-linear", "Linear")
title('Variação de q2')
xlabel('Tempo (s)')
ylabel('Posição (m)')

figure(3)
plot(tempo, y(:,3), "b")
% hold on
% plot(tempo, yL(:,3), "r")
% legend("Não-linear", "Linear")
title('Variação de theta em rad')
xlabel('Tempo (s)')
ylabel('Ângulo (rad)')

figure(9)
plot(tempo, y(:,10), "b")
% hold on
% plot(tempo, yL(:,8), "r")
% legend("Não-linear", "Linear")
title("Velocidade longitudinal")
xlabel('Tempo (s)')
ylabel('Velocidade (m/s)')


function dydt = ft(t, y_0)
    phi = 13*pi/180; 
    g = 9.8;
    rho = 1.2923; 
    S = (25.6^2)/1.7; 
    C_pav = 8; 
    M = 88000; 
    m = 2*3000;
    C_L = (-0.00165*(180*y_0(3)/pi)^2+ 0.07378*180*y_0(3)/pi + 0.21999);
    C_D = (0.00017*(180*y_0(3)/pi)^2 + 0.01111*180*y_0(3)/pi + 0.15714);
    u_rol = (0.0041+0.000041*y_0(2))*C_pav;
    D_po = 5;
    D_go = 2.2;
    D_fo = 20.3;
    u_v = 0;
    y_ext = 0;
    yponto_ext = 0;
    J_oz = 16864415;
    k_r = 2*13600000;
    c_r = 2*9700;
    k_t = 11486800;
    c_t = 1021960;
    m_f = 1000; 
    k_tf = 5743400;
    c_tf = 51098;
    k_rf = 3400000;
    c_rf = 2*2425;
    dydt1 = y_0(5);
    dydt2 = y_0(6);
    dydt3 = y_0(7); 
    dydt4 = y_0(8);
    dydt5 = -(((c_r + c_t)*y_0(5))/m) + (c_t*y_0(6))/m + (-(g*m) + k_t*(-y_0(1) + y_0(2)) + k_r*(-y_0(1) + y_ext) + c_r*yponto_ext)/m;
    dydt6 = (c_t*J_oz*y_0(5))/(M*(-(M*(cos(y_0(3))*cos(y_0(3)))*(D_go*D_go)) + J_oz)) + ((-(M*(cos(y_0(3))*cos(y_0(3)))*c_tf*D_fo*D_go) + (c_t + c_tf)*J_oz)*y_0(6))/(M*(M*(cos(y_0(3))*cos(y_0(3)))*(D_go*D_go) - J_oz)) - (M*cos(y_0(3))*(D_go*D_go)*(2*g*M*cos(y_0(3)) + sin(y_0(3))*(-2*g*M + S*rho*C_L*((y_0(10) + u_v)*(y_0(10) + u_v)))*u_rol) + J_oz*(S*rho*C_L*((y_0(10) + u_v)*(y_0(10) + u_v)) - 2*(g*M + k_t*(-y_0(1) + y_0(2))) - 2*k_tf*(sin(y_0(3))*D_fo + y_0(2) - y_0(4))) + M*cos(y_0(3))*D_go*(-(S*rho*(sin(y_0(3))*C_D + cos(y_0(3))*C_L)*D_po*((y_0(10) + u_v)*(y_0(10) + u_v))) + 2*cos(y_0(3))*D_fo*k_tf*(sin(y_0(3))*D_fo + y_0(2) - y_0(4))))/(2*M*(M*(cos(y_0(3))*cos(y_0(3)))*(D_go*D_go) - J_oz)) + (c_tf*(M*(cos(y_0(3))*cos(y_0(3)))*D_fo*D_go - J_oz)*y_0(8))/(M*(M*(cos(y_0(3))*cos(y_0(3)))*(D_go*D_go) - J_oz)) + (cos(y_0(3))*c_tf*D_fo*(-(M*(cos(y_0(3))*cos(y_0(3)))*D_fo*D_go) + J_oz)*y_0(7))/(M*(M*(cos(y_0(3))*cos(y_0(3)))*(D_go*D_go) - J_oz)) - (sin(y_0(3))*D_go*J_oz*(y_0(7)*y_0(7)))/(M*(cos(y_0(3))*cos(y_0(3)))*(D_go*D_go) - J_oz);
    dydt7 = (sec(y_0(3))*c_t*D_go*y_0(5))/(M*(D_go*D_go) - sec(y_0(3))*sec(y_0(3))*J_oz) + (cos(y_0(3))*(c_tf*D_fo - (c_t + c_tf)*D_go)*y_0(6))/(M*(cos(y_0(3))*cos(y_0(3)))*(D_go*D_go) - J_oz) + (sec(y_0(3))*(2*sin(y_0(3))*(D_fo*D_fo)*k_tf - (S*rho*C_D*D_po*((y_0(10) + u_v)*(y_0(10) + u_v)) + 2*g*M*D_go*u_rol)*tan(y_0(3)) + S*rho*C_L*((y_0(10) + u_v)*(y_0(10) + u_v))*(-D_po + D_go*(1 + u_rol*tan(y_0(3)))) - 2*D_fo*k_tf*(sin(y_0(3))*D_go - y_0(2) + y_0(4)) + 2*D_go*(k_t*(y_0(1) - y_0(2)) + k_tf*(-y_0(2) + y_0(4)))))/(2*M*(D_go*D_go) - 2*(sec(y_0(3))*sec(y_0(3)))*J_oz) + (cos(y_0(3))*c_tf*(D_fo - D_go)*y_0(8))/(-(M*(cos(y_0(3))*cos(y_0(3)))*(D_go*D_go)) + J_oz) + (c_tf*D_fo*(D_fo - D_go)*y_0(7))/(M*(D_go*D_go) - sec(y_0(3))*sec(y_0(3))*J_oz) + (y_0(7)*y_0(7))/(cot(y_0(3)) - (2*csc(2*y_0(3))*J_oz)/(M*(D_go*D_go)));
    dydt8 = (c_tf*y_0(6))/m_f - ((c_rf + c_tf)*y_0(8))/m_f + (-(g*m_f) + k_tf*(sin(y_0(3))*D_fo + y_0(2) - y_0(4)) + k_rf*(-y_0(4) + y_ext) + c_rf*yponto_ext)/m_f + (cos(y_0(3))*c_tf*D_fo*y_0(7))/m_f;
    dydt9 = y_0(10);
    dydt10 = (-(M*g - C_L*S*rho*y_0(10)^2/2)*(u_rol)*C_pav - C_D*S*rho*y_0(7)^2/2)/M;
    dydt = [dydt1; dydt2; dydt3; dydt4; dydt5; dydt6; dydt7; dydt8; dydt9; dydt10];
end