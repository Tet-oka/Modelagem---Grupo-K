
%Definição dos parâmetros iniciais da integração
t = 5; % Período referente a etapa de Free-Roll 
T_sim = 1/100; % Passo de integração
tempo = 0:T_sim:t; % Vetor tempo com passo de integração igual a T_sim
%passo máximo ode
max_step = odeset('MaxStep', T_sim);
x_0 = 0; %Condição de posição inicial [m]
x_ponto_0 = 300/3.6; %Condição de velocidade inicial [m/s]
x0 = [x_0 x_ponto_0];

%Aplicação de Runge-Kutta (4-5)
[t, y_v] = ode45(@f, tempo, x0, max_step);
plot(t, y_v(:,2))

%Defininido os espaços de estados e parâmetros
function dudt = f(t, u_0)
g = 9.8;
alpha = 0; 
rho = 1.2923; 
S = (25^2)/1.7; 
C_pav = 1; 
M = 91000; 
C_l = -0.00165*alpha^2+0.07378*alpha+0.21999;
C_d = 0.00017*alpha^2+0.01111*alpha+0.15714;
u_rol = 0.0041+0.000041*u_0(2);
dudt1 = u_0(2);
dudt2 = (-(M*g - C_l*S*rho*u_0(2)^2/2)*(u_rol)*C_pav - C_d*S*rho*u_0(2)^2/2)/M;
dudt =  [dudt1; dudt2];
end