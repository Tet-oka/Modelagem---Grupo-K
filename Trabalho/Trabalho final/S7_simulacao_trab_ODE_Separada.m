
%% Definição dos parâmetros iniciais da integração
t = 12; 
T_sim = 1/100;
tempo = 0:T_sim:t;
%Passo máximo ODE
max_step = odeset('MaxStep', T_sim);

phi = 13*pi/180; 
g = 9.81;
rho = 1.2923; 
M = 88000; 
WL = 285; %Wing loading para Concorde pousando com 110 ton com A = (25.6^2/1.7)
S = M/WL; 
C_pav = 1; 
m = 4000;
D_po = 5;
D_go = 3.5;
D_fo = 18.2;
D_co = 29.2;
J_oz = 16864415*M/88000;
k_r = 13600000;
c_r = 9700;
k_t = 11486800;
c_t = 102196;
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
thetap_0 = -0.1;
q3p_0 = 0;
pos_ini = 0;
vel_ini = 300/3.6;

x0 = [q1_0 q2_0 theta_0 q3_0 q1p_0 q2p_0 thetap_0 q3p_0 pos_ini vel_ini];

[t, y] = ode45(@f, tempo, x0, max_step);
[t, yL] = ode45(@fL, tempo, x0, max_step);

%% Condições para a troca da ODE entre modelo incial e meio carro

h=1;
while y(h,9)~=y(h+1,9)
    h=h+1;
end

j=1;
while yL(j,9)~=yL(j+1,9)
    j=j+1;
end

y0=y(h-1,:);
yL0=yL(j-1,:);

y=y(1:h-1,:);
yL=yL(1:j-1,:);

tempofr = tempo(1:(h-1));
tempofrL = tempo(1:(j-1));
tempot = tempo(h:end);
tempotL = tempo(j:end);

yL0(3)=yL0(3)+13*pi/180;
yL0(8)=yL0(7)*D_fo*cos(yL0(3));

xx0=y0;
xxL0=yL0;

[tt, yy] = ode45(@ff, tempot, xx0, max_step);
[tt, yyL] = ode45(@ffL, tempotL, xxL0, max_step);


for i = 1:length(yyL)
    yyL(i,3)=yyL(i,3)-13*pi/180;
end

yc = [y;yy];
ycL = [yL;yyL];

min1 = (min(yc(:,1)));
max1 = (max(yc(:,1)));
min2 = (min(yc(:,2)));
max2 = (max(yc(:,2)));
min3 = (min(yc(:,3)));
max3 = (max(yc(:,3)));
min4 = (min(yc(:,4)));
max4 = (max(yc(:,4)));
min5 = (min(yc(:,5)));
max5 = (max(yc(:,5)));
min6 = (min(yc(:,6)));
max6 = (max(yc(:,6)));
min7 = (min(yc(:,7)));
max7 = (max(yc(:,7)));
min8 = (min(yc(:,8)));
max8 = (max(yc(:,8)));
min9 = (min(yc(:,9)));
max9 = (max(yc(:,9)));
min10 = (min(yc(:,10)));
max10 = (max(yc(:,10)));


%% Cálculo das acelerações

dq2 = diff(yc(:,6))./diff(tempo); 
dq2L = diff(ycL(:,6))./diff(tempo);

dtheta = diff(yc(:,7))./diff(tempo);
dthetaL = diff(ycL(:,7))./diff(tempo);

G_2ponto = dq2 + (D_go.*dtheta.*cos(yc(2:end,3)));

Cockpit_2ponto = dq2 + (D_co.*dtheta.*cos(yc(2:end,3)));
Cockpit_2pontoL = dq2L + (D_co.*dthetaL.*cos(ycL(2:end,3)));
rmscockpit = rms(Cockpit_2ponto(1,:).')

% C_L = (-0.00165*((180*(y(:,3)+ phi)/pi).^2)) + (0.07378*180*(y(:,3)+ phi)/pi) + 0.21999;
% C_D = (0.00017*((180*(y(:,3)+ phi)/pi).^2)) + (0.01111*180*(y(:,3)+ phi)/pi) + 0.15714;
% L = C_L.*0.5.*S.*rho.*(y(:,10)).^2;



% figure(20)
% plot(tempo, C_L, "b")
% title("CL")
% xlabel('Tempo (s)')
% 
% figure(21)
% plot(tempo, L, "b")
% title("L")
% xlabel('Tempo (s)')

% figure(22)
% plot(tempo, C_D, "b")
% title("CD")
% xlabel('Tempo (s)')
% a_q1 = -(((c_r + c_t).*y(:,4))/m) + (c_t.*y(:,5))/m + (-(g*m) + k_t.*(-y(:,1) + y(:,2)) + k_r.*(-y(:,1) + y_ext + c_r*yponto_ext))/m;
% a_q2 = -0.5*(-2*g*M*J_oz + S*rho*C_L*J_oz.*((y(:,8) + u_v).*(y(:,8) + u_v)) + 2*J_oz*k_t.*(y(:,1) - y(:,2)) + M*cos(phi)*D_go*(D_go*(-((2*g*M - S*rho*C_L.*((y(:,8) + u_v).*(y(:,8) + u_v)))*u_rol*(sin(phi) + cos(phi).*y(:,3))) + 2*g*M*(cos(phi) - sin(phi).*y(:,3))) - S*rho*D_po.*((y(:,8) + u_v).*(y(:,8) + u_v))*(cos(phi)*(C_L + C_D.*y(:,3)) + sin(phi)*(C_D - C_L.*y(:,3)))))/(M*(M*(cos(phi)*cos(phi))*(D_go*D_go) - J_oz)) + (c_t*J_oz.*y(:,4))/(M*(-(M*(cos(phi)*cos(phi))*(D_go*D_go)) + J_oz)) + (c_t*J_oz.*y(:,5))/(M*M*(cos(phi)*cos(phi))*(D_go*D_go) - M*J_oz);

 %% Plot dos gráficos

% figure(1)
% plot(tempo, yc(:,1), "b")
% hold on
% plot(tempo, ycL(:,1), "r")
% grid on
% grid minor
% ylim([11/10*min1 11/10*max1]);
% p = patch([0 0 T_sim*length(tempofr) T_sim*length(tempofr)],[11/10*min1 11/10*max1 11/10*max1 11/10*min1],'');
% set(p,'FaceAlpha',0.1)
% set(p,'EdgeColor','none')
% legend("Não-linear", "Linear", "Free-roll")
% title('Variação de q1')
% xlabel('Tempo (s)')
% ylabel('Posição (m)')
% 
% figure(2)
% plot(tempo, yc(:,2), "b")
% hold on
% plot(tempo, ycL(:,2), "r")
% grid on
% grid minor
% ylim([11/10*min2 11/10*max2]);
% p = patch([0 0 T_sim*length(tempofr) T_sim*length(tempofr)],[11/10*min2 11/10*max2 11/10*max2 11/10*min2],'');
% set(p,'FaceAlpha',0.1)
% set(p,'EdgeColor','none')
% legend("Não-linear", "Linear", "Free-roll")
% title('Variação de q2')
% xlabel('Tempo (s)')
% ylabel('Posição (m)')
% 
% figure(3)
% plot(tempo, yc(:,3), "b")
% hold on
% plot(tempo, ycL(:,3), "r")
% grid on
% grid minor
% ylim([11/10*min3 11/10*max3]);
% p = patch([0 0 T_sim*length(tempofr) T_sim*length(tempofr)],[11/10*min3 11/10*max3 11/10*max3 11/10*min3],'');
% set(p,'FaceAlpha',0.1)
% set(p,'EdgeColor','none')
% legend("Não-linear", "Linear", "Free-roll")
% title('Variação de theta em rad')
% xlabel('Tempo (s)')
% ylabel('Ângulo (rad)')
% 
% figure(4)
% plot(tempo, 180*yc(:,3)/pi, "b")
% hold on
% plot(tempo, 180*ycL(:,3)/pi, "r")
% grid on
% grid minor
% ylim([11/10*min2*180/pi 11/10*max3*180/pi]);
% p = patch([0 0 T_sim*length(tempofr) T_sim*length(tempofr)],[11/10*min3*180/pi 11/10*max3*180/pi 11/10*max3*180/pi 11/10*min3*180/pi],'');
% set(p,'FaceAlpha',0.1)
% set(p,'EdgeColor','none')
% legend("Não-linear", "Linear", "Free-roll")
% title('Variação de theta em graus para Phi=13')
% xlabel('Tempo (s)')
% ylabel('Ângulo (graus)')
% 
% figure(5)
% plot(tempo, yc(:,5), "b")
% hold on
% plot(tempo, ycL(:,5), "r")
% grid on
% grid minor
% ylim([11/10*min5 11/10*max5]);
% p = patch([0 0 T_sim*length(tempofr) T_sim*length(tempofr)],[11/10*min5 11/10*max5 11/10*max5 11/10*min5],'');
% set(p,'FaceAlpha',0.1)
% set(p,'EdgeColor','none')
% legend("Não-linear", "Linear", "Free-roll")
% title('Velocidade de q1')
% xlabel('Tempo (s)')
% ylabel('Velocidade (m/s)')
% 
% figure(6)
% plot(tempo, yc(:,6), "b")
% hold on
% plot(tempo, ycL(:,6), "r")
% grid on
% grid minor
% ylim([11/10*min6 11/10*max6]);
% p = patch([0 0 T_sim*length(tempofr) T_sim*length(tempofr)],[11/10*min6 11/10*max6 11/10*max6 11/10*min6],'');
% set(p,'FaceAlpha',0.1)
% set(p,'EdgeColor','none')
% legend("Não-linear", "Linear", "Free-roll")
% title('Velocidade de q2')
% xlabel('Tempo (s)')
% ylabel('Velocidade (m/s)')
% 
% figure(7)
% plot(tempo, yc(:,7), "b")
% hold on
% plot(tempo, ycL(:,7), "r")
% grid on
% grid minor
% ylim([11/10*min7 11/10*max7]);
% p = patch([0 0 T_sim*length(tempofr) T_sim*length(tempofr)],[11/10*min7 11/10*max7 11/10*max7 11/10*min7],'');
% set(p,'FaceAlpha',0.1)
% set(p,'EdgeColor','none')
% legend("Não-linear", "Linear", "Free-roll")
% title('Velocidade angular de theta')
% xlabel('Tempo (s)')
% ylabel('Velocidade (rad/s)')
% 
% figure(8)
% plot(tempo, yc(:,9), "b")
% hold on
% plot(tempo, ycL(:,9), "r")
% grid on
% grid minor
% ylim([11/10*min9 11/10*max9]);
% p = patch([0 0 T_sim*length(tempofr) T_sim*length(tempofr)],[11/10*min9 11/10*max9 11/10*max9 11/10*min9],'');
% set(p,'FaceAlpha',0.1)
% set(p,'EdgeColor','none')
% legend("Não-linear", "Linear", "Free-roll")
% title("Deslocamento longitudinal")
% xlabel('Tempo (s)')
% ylabel('Posição (m)')
% 
% figure(9)
% plot(tempo, yc(:,10), "b")
% hold on
% plot(tempo, ycL(:,10), "r")
% grid on
% grid minor
% set(p,'FaceAlpha',0.1)
% set(p,'EdgeColor','none')
% legend("Não-linear", "Linear", "Free-roll")
% title("Velocidade longitudinal")
% xlabel('Tempo (s)')
% ylabel('Velocidade (m/s)')
% 
% figure(10)
% plot(tempo, (180/pi)*abs(yc(:,3)-ycL(:,3)), "b")
% grid on
% grid minor
% title("Diferença entre modelo linear e não linear para theta")
% p = patch([0 0 T_sim*length(tempofr) T_sim*length(tempofr)],[-0 0.5 0.5 0],'');
% set(p,'FaceAlpha',0.1)
% set(p,'EdgeColor','none')
% legend("", "Free-roll")
% xlabel('Tempo (s)')
% ylabel('Diferença (graus)')
% 
% figure(11)
% plot(tempo, abs(yc(:,2)-ycL(:,2)), "b")
% hold on
% plot(tempo, abs(yc(:,1)-ycL(:,1)), "b")
% grid on
% grid minor
% legend("q2", "q1")
% title("Diferença entre modelo linear e não linear para q1")
% xlabel('Tempo (s)')
% ylabel('Diferença (m)')
% 
% figure(12)
% plot(tempo, yc(:,1), "r")
% hold on
% plot(tempo, yc(:,2), "g")
% xlabel('Tempo (s)')
% ylabel('Posição (m)')
% legend('Variação de q1','Variação de q2')
% 
% figure(13)
% plot(tempo, yc(:,4), "b")
% hold on
% plot(tempo, ycL(:,4), "r")
% grid on
% grid minor
% legend("Não-linear", "Linear")
% title('Variação de q3')
% xlabel('Tempo (s)')
% ylabel('Posição (m)')
% 
% figure(14)
% plot(tempo, yc(:,8), "b")
% hold on
% plot(tempo, ycL(:,8), "r")
% grid on
% grid minor
% legend("Não-linear", "Linear")
% title('Variação de q3ponto')
% xlabel('Tempo (s)')
% ylabel('Velcidade (m/s)')
% 
% figure(15)
% plot(tempo(2:end), (dq2(:,1)./g), "b")
% grid on
% grid minor
% hold on
% plot(tempo, ycL(:,8), "r")
% legend("Não-linear", "Linear")
% p = patch([0 0 T_sim*length(tempofr) T_sim*length(tempofr)],[-3 3 3 -3],'');
% set(p,'FaceAlpha',0.1)
% set(p,'EdgeColor','none')
% title('Carga inercial em q2',['C_{pav}=',int2str(C_pav)])
% legend("", "Free-roll")
% xlabel('Tempo (s)')
% ylabel('Aceleração (m/s^2)')


%% Plot dos gráficos antes do toque

% figure(101)
% plot(tempofr, y(:,1), "b")
% hold on
% plot(tempofrL, yL(:,1), "r")
% grid on
% grid minor
% legend("Não-linear", "Linear")
% title('Variação de q1 durante o free-roll')
% xlabel('Tempo (s)')
% ylabel('Posição (m)')
% 
% figure(102)
% plot(tempofr, y(:,2), "b")
% hold on
% plot(tempofrL, yL(:,2), "r")
% grid on
% grid minor
% legend("Não-linear", "Linear")
% title('Variação de q2 durante o free-roll')
% xlabel('Tempo (s)')
% ylabel('Posição (m)')
% 
% figure(103)
% plot(tempofr, y(:,3), "b")
% hold on
% plot(tempofrL, yL(:,3), "r")
% grid on
% grid minor
% legend("Não-linear", "Linear")
% title('Variação de theta em rad durante o free-roll')
% xlabel('Tempo (s)')
% ylabel('Ângulo (rad)')
% 
% figure(104)
% plot(tempofr, 180*y(:,3)/pi, "b")
% hold on
% plot(tempofrL, 180*yL(:,3)/pi, "r")
% grid on
% grid minor
% legend("Não-linear", "Linear")
% title('Variação de theta em graus para Phi=13 durante o free-roll')
% xlabel('Tempo (s)')
% ylabel('Ângulo (graus)')
% 
% figure(105)
% plot(tempofr, y(:,5), "b")
% hold on
% plot(tempofrL, yL(:,5), "r")
% grid on
% grid minor
% legend("Não-linear", "Linear")
% title('Velocidade de q1 durante o free-roll')
% xlabel('Tempo (s)')
% ylabel('Velocidade (m/s)')
% 
% figure(106)
% plot(tempofr, y(:,6), "b")
% hold on
% plot(tempofrL, yL(:,6), "r")
% grid on
% grid minor
% legend("Não-linear", "Linear")
% title('Velocidade de q2 durante o free-roll')
% xlabel('Tempo (s)')
% ylabel('Velocidade (m/s)')
% 
% figure(107)
% plot(tempofr, y(:,7), "b")
% hold on
% plot(tempofrL, yL(:,7), "r")
% grid on
% grid minor
% legend("Não-linear", "Linear")
% title('Velocidade angular de theta durante o free-roll')
% xlabel('Tempo (s)')
% ylabel('Velocidade (rad/s)')
% 
% figure(108)
% plot(tempofr, y(:,9), "b")
% hold on
% plot(tempofrL, yL(:,9), "r")
% grid on
% grid minor
% legend("Não-linear", "Linear")
% title("Deslocamento longitudinal durante o free-roll")
% xlabel('Tempo (s)')
% ylabel('Posição (m)')
% 
% figure(109)
% plot(tempofr, y(:,10), "b")
% hold on
% plot(tempofrL, yL(:,10), "r")
% grid on
% grid minor
% legend("Não-linear", "Linear")
% title("Velocidade longitudinal durante o free-roll")
% xlabel('Tempo (s)')
% ylabel('Velocidade (m/s)')
% 
% 
% figure(110)
% plot(tempofrL, (180/pi)*abs(y(1:(j-1),3)-yL(1:(j-1),3)), "b")
% grid on
% grid minor
% title("Diferença entre modelo linear e não linear para theta durante o free-roll")
% xlabel('Tempo (s)')
% ylabel('Diferença (graus)')
% 
% figure(111)
% plot(tempofrL, abs(y(1:(j-1),1)-yL(1:(j-1),1)), "m")
% hold on
% plot(tempofrL, abs(y(1:(j-1),2)-yL(1:(j-1),2)), "g")
% grid on
% grid minor
% legend("q1", "q2")
% title("Diferença entre modelo linear e não linear para q1 e q2 durante o free-roll")
% xlabel('Tempo (s)')
% ylabel('Diferença (m)')
% 
% figure(112)
% plot(tempofr, y(:,1), "m")
% hold on
% plot(tempofr, y(:,2), "g")
% title("Variação de q1 e q2 durante o free-roll")
% xlabel('Tempo (s)')
% ylabel('Posição (m) durante o free-roll')
% legend('Variação de q1','Variação de q2')
% 
% Não fazem nada durante o free-roll
% figure(113)
% plot(tempofr, y(:,4), "b")
% hold on
% plot(tempofrL, yL(:,4), "r")
% grid on
% grid minor
% legend("Não-linear", "Linear")
% title('Variação de q3 durante o free-roll')
% xlabel('Tempo (s)')
% ylabel('Posição (m)')
% 
% figure(114)
% plot(tempofr, y(:,8), "b")
% hold on
% plot(tempofrL, yL(:,8), "r")
% grid on
% grid minor
% legend("Não-linear", "Linear")
% title('Variação de q3ponto durante o free-roll')
% xlabel('Tempo (s)')
% ylabel('Velcidade (m/s)')
% 
% figure(115)
% plot(tempofr(2:end), (dq2(1:(h-2),1)./g), "b")
% hold on
% plot(tempofrL(2:end), (dq2L(1:(j-2),1)./g), "r")
% grid on
% grid minor
% title('Carga inercial em q2 durante o free-roll',['C_{pav}=',int2str(C_pav)])
% legend("Não-linear", "Linear")
% xlabel('Tempo (s)')
% ylabel('Aceleração (m/s^2)')
% 
% figure(201)
% plot(tempofr(2:end), (dq2(1:(h-2),1)./g),   "m")
% hold on
% grid on
% grid minor
% title('Carga inercial em q2 durante o free-roll')
% legend("c_t=c_{t0}", "c_t=c_{t0}*0.6", "c_t=c_{t0}*1.4", "c_t=c_{t0}*1.8")
% xlabel('Tempo (s)')
% ylabel('Carga inercial')
% 
% 
% figure(202)
% plot(tempofr(2:end), (Cockpit_2ponto(1:(h-2),1)),   "m")
% hold on
% grid on
% grid minor
% title('Aceleração vertical no cockpit durante o free-roll')
% legend("c_t=c_{t0}", "c_t=c_{t0}*0.6", "c_t=c_{t0}*1.4", "c_t=c_{t0}*1.8")
% xlabel('Tempo (s)')
% ylabel('Aceleração (m/s^2)')
% 
% figure(203)
% plot(tempofr(2:end), (dtheta(1:(h-2),1).*(180/pi)),   "m")
% hold on
% grid on
% grid minor
% title('Aceleração angular de arfagem do avião')
% legend("c_t=c_{t0}", "c_t=c_{t0}*0.6", "c_t=c_{t0}*1.4", "c_t=c_{t0}*1.8")
% xlabel('Tempo (s)')
% ylabel('Aceleração angular (graus/s^2)')
% 
% figure(204)
% plot(tempofr, (phi*180/pi)+180*y(:,3)/pi,   "m")
% hold on
% grid on
% grid minor
% legend("c_t=c_{t0}", "c_t=c_{t0}*0.6", "c_t=c_{t0}*1.4", "c_t=c_{t0}*1.8")
% title("Variação de theta em graus para Phi= "+ phi*180/pi +" durante o free-roll")
% xlabel('Tempo (s)')
% ylabel('Ângulo (graus)')
% 
% figure(205)
% plot(tempofr,  y(:,2), "r")
% hold on
% grid on
% grid minor
% legend("phi_0=13°", "phi_0=11°", "phi_0=9°", "phi_0=7°")
% title("Variação de q_2 durante o free-roll para diferentes ângulos de arfagem iniciais")
% xlabel('Tempo (s)')
% ylabel('Posição (m)')

figure(206)
plot(tempofrL,  yL(:,2), "m")
hold on
grid on
grid minor
legend("theta ponto=0", "theta ponto=-0.05", "theta ponto=-0.1", "theta ponto=-0.15")
title("Variação de q_2 durante o free-roll para diferentes velocidades angulares iniciais")
xlabel('Tempo (s)')
ylabel('Posição (m)')




%% Funções ODE

function dydt = f(t, y_0)
phi = 13*pi/180; 
g = 9.81;
rho = 1.2923; 
C_pav = 1; 
M = 88000; 
m = 4000;
WL = 285;
S = M/WL; 
%C_L = 0.0398*180*(y_0(3)+ phi)/pi -0.00633;
C_L = (0.04264*((180/pi)*(y_0(3) + phi)))-0.0158;
% C_D = 0.0088*180*(y_0(3)+ phi)/pi -0.00767;
C_D = (0.01394*((180/pi)*(y_0(3) + phi))) - 0.0248;
u_rol = (0.0041+0.000041*y_0(10))*C_pav;
D_po = 5;
D_go = 3.5;
D_fo = 18.2;
u_v = 0; %estoura em 3.9
y_ext = 0;
yponto_ext = 0;
J_oz = 16864415*M/88000;
k_r = 13600000;
c_r = 9700;
k_t = 11486800;
c_t = 102196;
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
    dydt1 = 0;
    dydt2 = 0;
    dydt3 = 0;
    dydt4 = 0;
    dydt5 = 0;
    dydt6 = 0;
    dydt7 = 0;
    dydt8 = 0;
    dydt9 = 0;
    dydt10 = 0;
    dydt = [dydt1; dydt2; dydt3; dydt4; dydt5; dydt6; dydt7; dydt8; dydt9; dydt10];
end
end

function dyLdt = fL(t, yL_0)
phi = 13*pi/180; 
g = 9.81;
rho = 1.2923; 
C_pav = 1; 
M = 88000; 
m = 4000;
WL = 285;
S = M/WL; 
%C_L = 0.0398*180*(yL_0(3)+ phi)/pi -0.00633;
C_L = (0.04264*((180/pi)*(yL_0(3) + phi)))-0.0158;
% C_D = 0.0088*180*(yL_0(3)+ phi)/pi -0.00767;
C_D = (0.01394*((180/pi)*(yL_0(3) + phi))) - 0.0248;
u_rol = (0.0041+0.000041*yL_0(10))*C_pav;
D_po = 5;
D_go = 3.5;
D_fo = 18.2;
u_v = 0; %estoura em 3.9
y_ext = 0;
yponto_ext = 0;
J_oz = 16864415*M/88000;
k_r = 13600000;
c_r = 9700;
k_t = 11486800;
c_t = 102196;
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
    dyLdt5 = -(((c_r+c_t)*yL_0(5))/m)+(c_t*yL_0(6))/m+(k_t*(-yL_0(1)+yL_0(2))+k_r*(-yL_0(1)+y_ext)+c_r*yponto_ext)/m;
    dyLdt6 = -0.5*(S*rho*C_L*J_oz*((yL_0(10)+u_v)*(yL_0(10)+u_v))+2*J_oz*k_t*(yL_0(1)-yL_0(2))+M*cos(phi)*D_go*(D_go*(-((2*g*M-S*rho*C_L*((yL_0(10)+u_v)*(yL_0(10)+u_v)))*u_rol*(sin(phi)+cos(phi)*yL_0(3)))+2*g*M*(cos(phi)-sin(phi)*yL_0(3)))-S*rho*D_po*((yL_0(10)+u_v)*(yL_0(10)+u_v))*(cos(phi)*(C_L+C_D*yL_0(3))+sin(phi)*(C_D-C_L*yL_0(3)))))/(M*(M*(cos(phi)*cos(phi))*(D_go*D_go)-J_oz))+(c_t*J_oz*yL_0(5))/(M*(-(M*(cos(phi)*cos(phi))*(D_go*D_go))+J_oz))+(c_t*J_oz*yL_0(6))/(M*M*(cos(phi)*cos(phi))*(D_go*D_go)-M*J_oz);
    dyLdt7 = -((sec(phi)*(S*rho*D_po*((yL_0(10)+u_v)*(yL_0(10)+u_v))*(C_D*(tan(phi)+yL_0(3))+C_L*(1-tan(phi)*yL_0(3)))+D_go*(-(S*rho*C_L*((yL_0(10)+u_v)*(yL_0(10)+u_v))*(1+u_rol*(tan(phi)+yL_0(3))))+2*(-(g*M)+k_t*(-yL_0(1)+yL_0(2))+g*M*tan(phi)*yL_0(3)+g*M*u_rol*(tan(phi)+yL_0(3))))))/(2*M*(D_go*D_go)-2*(sec(phi)*sec(phi))*J_oz))+(sec(phi)*c_t*D_go*yL_0(5))/(M*(D_go*D_go)-sec(phi)*sec(phi)*J_oz)+(sec(phi)*c_t*D_go*yL_0(6))/(-(M*(D_go*D_go))+sec(phi)*sec(phi)*J_oz);
    dyLdt8 = 0;
    dyLdt9 = yL_0(10);
    dyLdt10 = (-((M+m)*g - C_L*S*rho*((yL_0(10)+u_v)^2)/2)*(u_rol) - C_D*S*rho*((yL_0(10)+u_v)^2)/2)/(M+m);
    dyLdt = [dyLdt1; dyLdt2; dyLdt3; dyLdt4; dyLdt5; dyLdt6; dyLdt7; dyLdt8; dyLdt9; dyLdt10];
else
    dyLdt1 = 0;
    dyLdt2 = 0;
    dyLdt3 = 0;
    dyLdt4 = 0;
    dyLdt5 = 0;
    dyLdt6 = 0;
    dyLdt7 = 0;
    dyLdt8 = 0;
    dyLdt9 = 0;
    dyLdt10 = 0;
    dyLdt = [dyLdt1; dyLdt2; dyLdt3; dyLdt4; dyLdt5; dyLdt6; dyLdt7; dyLdt8; dyLdt9; dyLdt10];
end
end

function dydt = ff(t, y_0)
phi = 13*pi/180; 
g = 9.81;
rho = 1.2923; 
M = 88000; 
m = 4000;
WL = 285;
S = M/WL; 
C_pav = 1; 
%C_L = 0.0398*180*(y_0(3)+ phi)/pi -0.00633;
C_L = (0.04264*((180/pi)*(y_0(3) + phi)))-0.0158;
% C_D = 0.0088*180*(y_0(3)+ phi)/pi -0.00767;
C_D = (0.01394*((180/pi)*(y_0(3) + phi))) - 0.0248;
u_rol = (0.0041+0.000041*y_0(10))*C_pav;
D_po = 5;
D_go = 3.5;
D_fo = 18.2;
u_v = 0; %estoura em 3.9 
y_ext = 0;
yponto_ext = 0;
J_oz = 16864415*M/88000;
k_r = 13600000;
c_r = 9700;
k_t = 11486800;
c_t = 102196;
m_f = 1000; 
k_tf = 5743400;
c_tf = 51098;
k_rf = 3400000;
c_rf = 2*2425;

    dydt1 = y_0(5);
    dydt2 = y_0(6);
    dydt3 = y_0(7);
    dydt4 = y_0(8);
    dydt5 = -(((c_r + c_t)*y_0(5))/m) + (c_t*y_0(6))/m + (k_t*(-y_0(1) + y_0(2)) + k_r*(-y_0(1) + y_ext) + c_r*yponto_ext)/m;
    dydt6 = (c_t*J_oz*y_0(5))/(M*(-(M*(cos(phi + y_0(3))*cos(phi + y_0(3)))*(D_go*D_go)) + J_oz)) + ((-(M*(cos(phi + y_0(3))*cos(phi + y_0(3)))*c_tf*D_fo*D_go) + (c_t + c_tf)*J_oz)*y_0(6))/(M*(M*(cos(phi + y_0(3))*cos(phi + y_0(3)))*(D_go*D_go) - J_oz)) + (cos(phi + y_0(3))*c_tf*D_fo*(-(M*(cos(phi + y_0(3))*cos(phi + y_0(3)))*D_fo*D_go) + J_oz)*y_0(7))/(M*(M*(cos(phi + y_0(3))*cos(phi + y_0(3)))*(D_go*D_go) - J_oz)) - (sin(phi + y_0(3))*D_go*J_oz*(y_0(7)*y_0(7)))/(M*(cos(phi + y_0(3))*cos(phi + y_0(3)))*(D_go*D_go) - J_oz) - (M*cos(phi + y_0(3))*(D_go*D_go)*(2*g*M*cos(phi + y_0(3)) + sin(phi + y_0(3))*(-2*g*M + S*rho*C_L*((y_0(10) + u_v)*(y_0(10) + u_v)))*u_rol) + J_oz*(S*rho*C_L*((y_0(10) + u_v)*(y_0(10) + u_v)) + 2*k_t*(y_0(1) - y_0(2)) - 2*k_tf*(sin(phi + y_0(3))*D_fo + y_0(2) - y_0(4))) + M*cos(phi + y_0(3))*D_go*(-(S*rho*(sin(phi + y_0(3))*C_D + cos(phi + y_0(3))*C_L)*D_po*((y_0(10) + u_v)*(y_0(10) + u_v))) + 2*cos(phi + y_0(3))*D_fo*k_tf*(sin(phi + y_0(3))*D_fo + y_0(2) - y_0(4))))/(2.*M*(M*(cos(phi + y_0(3))*cos(phi + y_0(3)))*(D_go*D_go) - J_oz)) + (c_tf*(M*(cos(phi + y_0(3))*cos(phi + y_0(3)))*D_fo*D_go - J_oz)*y_0(8))/(M*(M*(cos(phi + y_0(3))*cos(phi + y_0(3)))*(D_go*D_go) - J_oz));
    dydt7 = (sec(phi + y_0(3))*c_t*D_go*y_0(5))/(M*(D_go*D_go) - sec(phi + y_0(3))*sec(phi + y_0(3))*J_oz) + (cos(phi + y_0(3))*(c_tf*D_fo - (c_t + c_tf)*D_go)*y_0(6))/(M*(cos(phi + y_0(3))*cos(phi + y_0(3)))*(D_go*D_go) - J_oz) + (c_tf*D_fo*(D_fo - D_go)*y_0(7))/(M*(D_go*D_go) - sec(phi + y_0(3))*sec(phi + y_0(3))*J_oz) + (M*(D_go*D_go)*tan(phi + y_0(3))*(y_0(7)*y_0(7)))/(M*(D_go*D_go) - sec(phi + y_0(3))*sec(phi + y_0(3))*J_oz) - (sec(phi + y_0(3))*(S*rho*D_po*((y_0(10) + u_v)*(y_0(10) + u_v))*(C_L + C_D*tan(phi + y_0(3))) + D_go*(-2*g*M + 2*g*M*u_rol*tan(phi + y_0(3)) - S*rho*C_L*((y_0(10) + u_v)*(y_0(10) + u_v))*(1 + u_rol*tan(phi + y_0(3))) + 2*k_t*(-y_0(1) + y_0(2)) + 2*k_tf*(sin(phi + y_0(3))*D_fo + y_0(2) - y_0(4))) - 2*D_fo*k_tf*(sin(phi + y_0(3))*D_fo + y_0(2) - y_0(4))))/(2*M*(D_go*D_go) - 2*(sec(phi + y_0(3))*sec(phi + y_0(3)))*J_oz) + (cos(phi + y_0(3))*c_tf*(D_fo - D_go)*y_0(8))/(-(M*(cos(phi + y_0(3))*cos(phi + y_0(3)))*(D_go*D_go)) + J_oz);
    dydt8 = (c_tf*y_0(6))/m_f + (cos(phi + y_0(3))*c_tf*D_fo*y_0(7))/m_f - ((c_rf + c_tf)*y_0(8))/m_f + (k_tf*(sin(phi + y_0(3))*D_fo + y_0(2) - y_0(4)) + k_rf*(-y_0(4) + y_ext) + c_rf*yponto_ext)/m_f;
    dydt9 = y_0(10);
    dydt10 = (-((M+m)*g - C_L*S*rho*((y_0(10)+u_v)^2)/2)*(u_rol) - C_D*S*rho*((y_0(10)+u_v)^2)/2)/(M+m);
    dydt = [dydt1; dydt2; dydt3; dydt4; dydt5; dydt6; dydt7; dydt8; dydt9; dydt10];

end

function dyLdt = ffL(t, yL_0)
phi = 13*pi/180; 
g = 9.81;
rho = 1.2923; 
M = 88000; 
m = 4000;
WL = 285;
S = M/WL; 
C_pav = 1; 
%C_L = 0.0398*180*(yL_0(3)+ phi)/pi -0.00633;
C_L = (0.04264*((180/pi)*(yL_0(3) + phi)))-0.0158;
% C_D = 0.0088*180*(yL_0(3)+ phi)/pi -0.00767;
C_D = (0.01394*((180/pi)*(yL_0(3) + phi))) - 0.0248;
u_rol = (0.0041+0.000041*yL_0(10))*C_pav;
D_po = 5;
D_go = 3.5;
D_fo = 18.2;
u_v = 0; %estoura em 3.9
y_ext = 0;
yponto_ext = 0;
J_oz = 16864415*M/88000;
k_r = 13600000;
c_r = 9700;
k_t = 11486800;
c_t = 102196;
m_f = 1000; 
k_tf = 5743400;
c_tf = 51098;
k_rf = 3400000;
c_rf = 2*2425;
vel_ini = 300/3.6;

    dyLdt1 = yL_0(5);
    dyLdt2 = yL_0(6);
    dyLdt3 = yL_0(7); 
    dyLdt4 = yL_0(8);
    dyLdt5 = -(((c_r+c_t)*yL_0(5))/m)+(c_t*yL_0(6))/m+(k_t*(-yL_0(1)+yL_0(2))+k_r*(-yL_0(1)+y_ext)+c_r*yponto_ext)/m;
    dyLdt6 = (c_t*J_oz*yL_0(5))/(M*(-(M*(D_go*D_go))+J_oz))+((-(M*c_tf*D_fo*D_go)+(c_t+c_tf)*J_oz)*yL_0(6))/(M*(M*(D_go*D_go)-J_oz))+(c_tf*D_fo*(-(M*D_fo*D_go)+J_oz)*yL_0(7))/(M*(M*(D_go*D_go)-J_oz))+(M*(D_go*D_go)*(-2*g*M+(2*g*M-S*rho*C_L*((yL_0(10)+u_v)*(yL_0(10)+u_v)))*u_rol*yL_0(3))-J_oz*(S*rho*C_L*((yL_0(10)+u_v)*(yL_0(10)+u_v))-2*(k_t*(-yL_0(1)+yL_0(2))+k_tf*(yL_0(2)+D_fo*yL_0(3)-yL_0(4))))+M*D_go*(S*rho*D_po*((yL_0(10)+u_v)*(yL_0(10)+u_v))*(C_L+C_D*yL_0(3))-2*D_fo*k_tf*(yL_0(2)+D_fo*yL_0(3)-yL_0(4))))/(2.*M*(M*(D_go*D_go)-J_oz))+(c_tf*(M*D_fo*D_go-J_oz)*yL_0(8))/(M*(M*(D_go*D_go)-J_oz));
    dyLdt7 = (c_t*D_go*yL_0(5))/(M*(D_go*D_go)-J_oz)+((c_tf*D_fo-(c_t+c_tf)*D_go)*yL_0(6))/(M*(D_go*D_go)-J_oz)+(c_tf*D_fo*(D_fo-D_go)*yL_0(7))/(M*(D_go*D_go)-J_oz)+(-(S*rho*D_po*((yL_0(10)+u_v)*(yL_0(10)+u_v))*(C_L+C_D*yL_0(3)))+D_go*(S*rho*C_L*((yL_0(10)+u_v)*(yL_0(10)+u_v))*(1+u_rol*yL_0(3))-2*(-(g*M)+k_t*(-yL_0(1)+yL_0(2))+D_fo*k_tf*yL_0(3)+g*M*u_rol*yL_0(3)+k_tf*(yL_0(2)-yL_0(4))))+2*D_fo*k_tf*(yL_0(2)+D_fo*yL_0(3)-yL_0(4)))/(2*M*(D_go*D_go)-2*J_oz)+(c_tf*(D_fo-D_go)*yL_0(8))/(-(M*(D_go*D_go))+J_oz);
    dyLdt8 = (c_tf*yL_0(6))/m_f+(c_tf*D_fo*yL_0(7))/m_f-((c_rf+c_tf)*yL_0(8))/m_f+(k_tf*(yL_0(2)+D_fo*yL_0(3)-yL_0(4))+k_rf*(-yL_0(4)+y_ext)+c_rf*yponto_ext)/m_f;
    dyLdt9 = yL_0(10);
    dyLdt10 = (-((M+m)*g - C_L*S*rho*((yL_0(10)+u_v)^2)/2)*(u_rol) - C_D*S*rho*((yL_0(10)+u_v)^2)/2)/(M+m);
    dyLdt = [dyLdt1; dyLdt2; dyLdt3; dyLdt4; dyLdt5; dyLdt6; dyLdt7; dyLdt8; dyLdt9; dyLdt10];

end

