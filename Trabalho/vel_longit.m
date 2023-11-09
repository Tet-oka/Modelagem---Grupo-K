

%% Definição dos parâmetros iniciais da integração

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
x0 = [x_0 x_ponto_0];

%% Aplicação de Runge-Kutta (4-5)

[t, y] = ode45(@f_n_lin, tempo, x0, max_step);
figure(2)
plot(t, y(:,2),"b")
%% Defininido os espaços de estados
%Não linearizado
function dxdt = f_n_lin(t, y_0)
g = 9.8;
alpha = 0; 
rho = 1.2923; 
S = (25^2)/1.7; 
C_pav = 1; 
M = 91000; 
C_l = -0.00165*alpha^2+0.07378*alpha+0.21999;
C_d = 0.00017*alpha^2+0.01111*alpha+0.15714;
dxdt_1 = y_0(2);
dxdt_2 = (-(M*g - C_l*S*rho*y_0(2)^2/2)*(0.0041+0.000041*y_0(2))*C_pav - C_d*S*rho*y_0(2)^2/2)/M;
dxdt =  [dxdt_1; dxdt_2];
end