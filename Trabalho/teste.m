%% Parâmetros do problema

g = 9.8;
alpha = 0; % 0 a 25°
rho = 0; % A definir
S = 0; % A definir
C_pav = 1; % 1 para concreto liso, 1.2 para concreto gasto e 1.5 para asfalto quente
M = 88*1000; % Massa máxima de aterrissagem [Kg]

% Onde tem "vel" precisa ajeitar pra resposta da edo

% Coeficientes de Sustentação e Arrasto
C_l = 0.00165*alpha^2+0.07378*alpha+0.21999;
C_d = 0.00017*alpha^2+0.01111*alpha+0.15714;
%Mi_rol = (0.0041+0.000041*vel)*C_pav; 

% Cálculo das forças
% L = C_l*S*rho*vel^2/2;
% D = C_d*S*rho*vel^2/2;
% N = M*g-L;
% F_rol = Mi_rol*N;


%% Definição dos parâmetros iniciais da integração

X_0 = [0 0 11*pi/180 -3 -3 0] %[q1 q2 \theta q1_ponto q2_ponto \theta_ponto]

t = 2; % Período referente a etapa de Free-Roll 
T_sim = 1/100; % Passo de integração
tempo = 0:T_sim:t; % Vetor tempo com passo de integração igual a T_sim

% número pontos de integração para o Euler Explícito
q = size(tempo(1,:));
q = q(2);

%passo máximo ode
max_step = odeset('MaxStep', T_sim);

x_0 = 0; %Condição de posição inicial [m]
x_ponto_0 = 300/3.6; %Condição de velocidade inicial [m/s]

%% Aplicação de Runge-Kutta (4-5)

[instante_tempo , X_simulado] = ode45(@f_n_lin, tempo, X_0, max_step);
plot(instante_tempo, X_simulado_1)

%% Defininido os espaços de estados

function dxdt = f_n_lin(t, x)
J_oz = 16864415.983;
g = 9.8;
u_long = 300/3.6; % Velocidade longitudinal
u_v = 0; % Velocidade do vento
u_rol = u_long;
rho = 1.293; % Densidade do ar
D_go = 2.2; % Distência GO
D_po = 5; % Distância PO
S = 27;
alpha = 11; % 0 a 25°
phi = alpha;
C_pav = 1; % 1 para concreto liso, 1.2 para concreto gasto e 1.5 para asfalto quente
M = 88*1000; % Massa máxima de aterrissagem [Kg]
C_L = 0.00165*alpha^2+0.07378*alpha+0.21999;
C_d = 0.00017*alpha^2+0.01111*alpha+0.15714;
Mi_rol = (0.0041+0.000041*u_long)*C_pav; 
m = 4000;
k_t = 2913000;
k_r = 4418000;
c_t = 17400;
c_r = 8800000; % Constante amortecimento trem pouso
y_ext=0;
yponto_ext = 0;
dxdt_1 = x(4);
dxdt_2 = x(5);
dxdt_3 = x(6);
dxdt_4 = -((((c_r + c_r)*x(4))/m) + (c_r*x(5))/m + (-(g*m) + k_t*(-x(1) + x(2)) + k_r*(-x(1) + y_ext) + c_r*yponto_ext))/m;
dxdt_5 = -(((g*M - (S*rho*C_L*((u_long + u_v)*(u_long + u_v)))/2. + k_t*(-x(1) + x(2)))/M) + (c_r*x(4))/M - (c_r*x(5)))/M;
dxdt_6 = ((S*rho*sin(phi + x(3))*C_d*D_po*((u_long + u_v)*(u_long + u_v)) - S*rho*C_L*((u_long + u_v)*(u_long + u_v))*(cos(phi + x(3))*(D_go - D_po) + sin(phi + x(3))*D_go*u_rol) + 2*D_go*(g*M*sin(phi + x(3))*u_rol + cos(phi + x(3))*k_t*(-x(1) + x(2))))/(2*J_oz) - (cos(phi + x(3))*c_r*D_go*x(4))/J_oz + (cos(phi + x(3))*c_r*D_go*x(5)))/J_oz;
dxdt = [dxdt_1; dxdt_2; dxdt_3; dxdt_4; dxdt_5; dxdt_6]
end