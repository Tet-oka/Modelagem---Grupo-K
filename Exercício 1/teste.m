%% Parâmetros do problema

% Dados do enunciado
wp = 1;
g = 9.8;
L1 = g/(wp^2);
L2 = L1*20/21;
L = L1/L2;
mi = 2;                             % definido arbitrariamente
m1 = mi*L1;
m2 = mi*L2;


%% Definição dos parâmetros iniciais da integração

% Ao substituir os valor de \lambda = (21/20) e de \omega_p = 1, teremos que
%a freqência natural no primeiro modo é de 0,4Hz. 

% O período de simulação utilizado será de 4 períodos completos de
%oscilação, que deve ser suficiente para verificar diferenças entre os
%métodos de integração e entre as equações linearizadas e não linearizadas
F_1modo = 0.4;
t = 1/F_1modo;
t = 20*t;

% Se usarmos 100 passos de integração por segundo:
T_sim = 1/100;

% criando o vetor tempo com 100 passos por segundo e 4 períodos de oscilação
tempo = 0:T_sim:t;

% número pontos de integração
q = size(tempo(1,:));
q = q(2);

% condições iniciais de integração:

% Cenário 1:
y_0_1 = [10*pi/180 0*pi/180 0 0];

% Cenário 2
y_0_2 = [75*pi/180 0*pi/180 0 0];


%% Aplicação de Runge-Kutta (4-5) em C1 e C2

y_0 = y_0_1;    % aplicando para o cenário 1
[t_runge_n_lin_C1, y_runge_n_lin_C1] = ode45(@f_n_lin, tempo, y_0);     % não linearizado
[t_runge_lin_C1, y_runge_lin_C1] = ode45(@f_lin, tempo, y_0);           % linearizado

y_0 = y_0_2;     % aplicando para o cenário 2
[t_runge_n_lin_C2, y_runge_n_lin_C2] = ode45(@f_n_lin, tempo, y_0);     % não linearizado
[t_runge_lin_C2, y_runge_lin_C2] = ode45(@f_lin, tempo, y_0);           % linearizado

%% Aplicação do Euler Explícito em C1 e C2

% criar as matrizes que receberão os valores calculados
y_euler_n_lin_C1 = zeros(q,4);       %C1 não linearizado
y_euler_lin_C1 = zeros(q,4);         %C1 linearizado
y_euler_n_lin_C2 = zeros(q,4);       %C2 não linearizado
y_euler_lin_C2 = zeros(q,4);         %C2 linearizado

% para o cenário C1
y_0 = y_0_1;

% aplicando o Euler Explícito para não linearizado
for i = 0:q-1

    % valores de f
    dydt_n_lin_1 = y_0(3);
    dydt_n_lin_2 = y_0(4);
    dydt_n_lin_3 = (-3*g*((4*m1 + 5*m2)*sin(y_0(1)) + 3*m2*sin(y_0(1) - 2*y_0(2))))/((8*m1 + 15*m2 - 9*m2*cos(2*(y_0(1) - y_0(2))))*L1) + (9*m2*sin(2*(y_0(1) - y_0(2)))*(y_0(3)*y_0(3)))/(-8*m1 - 15*m2 + 9*m2*cos(2*(y_0(1) - y_0(2)))) + (6*m2*sin(y_0(1) - y_0(2))*L2*(y_0(4)*y_0(4)))/((-4*(m1 + 3*m2) + 9*m2*(cos(y_0(1) - y_0(2))*cos(y_0(1) - y_0(2))))*L1);
    dydt_n_lin_4 = (9*g*(m1 + 2*m2)*sin(2*y_0(1) - y_0(2)) - 3*g*(m1 + 6*m2)*sin(y_0(2)))/((8*m1 + 15*m2 - 9*m2*cos(2*(y_0(1) - y_0(2))))*L2) + (6*(m1 + 3*m2)*sin(y_0(1) - y_0(2))*L1*(y_0(3)*y_0(3)))/((4*(m1 + 3*m2) - 9*m2*(cos(y_0(1) - y_0(2))*cos(y_0(1) - y_0(2))))*L2) + (9*m2*sin(2*(y_0(1) - y_0(2)))*(y_0(4)*y_0(4)))/(8*m1 + 15*m2 - 9*m2*cos(2*(y_0(1) - y_0(2))));

    % y(i+1) = y(i) + T_sim * f
    y_euler_n_lin_C1(i+1,1) = y_0(1,1) + T_sim*dydt_n_lin_1;
    y_euler_n_lin_C1(i+1,2) = y_0(1,2) + T_sim*dydt_n_lin_2;
    y_euler_n_lin_C1(i+1,3) = y_0(1,3) + T_sim*dydt_n_lin_3;
    y_euler_n_lin_C1(i+1,4) = y_0(1,4) + T_sim*dydt_n_lin_4;

    y_0 = [ y_euler_n_lin_C1(i+1,1)  y_euler_n_lin_C1(i+1,2)  y_euler_n_lin_C1(i+1,3)  y_euler_n_lin_C1(i+1,4)];
end

y_0 = y_0_1;
% aplicando o Euler Explícito para linearizado
for i = 0:q-1

    % valores de f
    dydt_lin_1 = y_0(3);
    dydt_lin_2 = y_0(4);
    dydt_lin_3 = (3*(wp^2)*((-2*(1+2*L)*y_0(1))+(L*(3+2*L)*y_0(2))))/(4+3*L);
    dydt_lin_4 = (3*(wp^2)*(((3+4*L*(2+L))*y_0(1))-(2*((1+L)^3)*y_0(2))))/(L*(4+3*L));

    % y(i+1) = y(i) + T_sim * f
    y_euler_lin_C1(i+1,1) = y_0(1,1) + T_sim*dydt_lin_1;
    y_euler_lin_C1(i+1,2) = y_0(1,2) + T_sim*dydt_lin_2;
    y_euler_lin_C1(i+1,3) = y_0(1,3) + T_sim*dydt_lin_3;
    y_euler_lin_C1(i+1,4) = y_0(1,4) + T_sim*dydt_lin_4;

    y_0 = [ y_euler_lin_C1(i+1,1)  y_euler_lin_C1(i+1,2)  y_euler_lin_C1(i+1,3)  y_euler_lin_C1(i+1,4)];
end

% para o cenário C2
y_0 = y_0_2;

% aplicando o Euler Explícito para não linearizado
for i = 0:q-1

    % valores de f
    dydt_n_lin_1 = y_0(3);
    dydt_n_lin_2 = y_0(4);
    dydt_n_lin_3 = (-3*g*((4*m1 + 5*m2)*sin(y_0(1)) + 3*m2*sin(y_0(1) - 2*y_0(2))))/((8*m1 + 15*m2 - 9*m2*cos(2*(y_0(1) - y_0(2))))*L1) + (9*m2*sin(2*(y_0(1) - y_0(2)))*(y_0(3)*y_0(3)))/(-8*m1 - 15*m2 + 9*m2*cos(2*(y_0(1) - y_0(2)))) + (6*m2*sin(y_0(1) - y_0(2))*L2*(y_0(4)*y_0(4)))/((-4*(m1 + 3*m2) + 9*m2*(cos(y_0(1) - y_0(2))*cos(y_0(1) - y_0(2))))*L1);
    dydt_n_lin_4 = (9*g*(m1 + 2*m2)*sin(2*y_0(1) - y_0(2)) - 3*g*(m1 + 6*m2)*sin(y_0(2)))/((8*m1 + 15*m2 - 9*m2*cos(2*(y_0(1) - y_0(2))))*L2) + (6*(m1 + 3*m2)*sin(y_0(1) - y_0(2))*L1*(y_0(3)*y_0(3)))/((4*(m1 + 3*m2) - 9*m2*(cos(y_0(1) - y_0(2))*cos(y_0(1) - y_0(2))))*L2) + (9*m2*sin(2*(y_0(1) - y_0(2)))*(y_0(4)*y_0(4)))/(8*m1 + 15*m2 - 9*m2*cos(2*(y_0(1) - y_0(2))));

    % y(i+1) = y(i) + T_sim * f
    y_euler_n_lin_C2(i+1,1) = y_0(1,1) + T_sim*dydt_n_lin_1;
    y_euler_n_lin_C2(i+1,2) = y_0(1,2) + T_sim*dydt_n_lin_2;
    y_euler_n_lin_C2(i+1,3) = y_0(1,3) + T_sim*dydt_n_lin_3;
    y_euler_n_lin_C2(i+1,4) = y_0(1,4) + T_sim*dydt_n_lin_4;

    y_0 = [ y_euler_n_lin_C2(i+1,1)  y_euler_n_lin_C2(i+1,2)  y_euler_n_lin_C2(i+1,3)  y_euler_n_lin_C2(i+1,4)];
end

y_0 = y_0_2;
% aplicando o Euler Explícito para linearizado
for i = 0:q-1

    % valores de f
    dydt_lin_1 = y_0(3);
    dydt_lin_2 = y_0(4);
    dydt_lin_3 = (3*(wp^2)*((-2*(1+2*L)*y_0(1))+(L*(3+2*L)*y_0(2))))/(4+3*L);
    dydt_lin_4 = (3*(wp^2)*(((3+4*L*(2+L))*y_0(1))-(2*((1+L)^3)*y_0(2))))/(L*(4+3*L));

    % y(i+1) = y(i) + T_sim * f
    y_euler_lin_C2(i+1,1) = y_0(1,1) + T_sim*dydt_lin_1;
    y_euler_lin_C2(i+1,2) = y_0(1,2) + T_sim*dydt_lin_2;
    y_euler_lin_C2(i+1,3) = y_0(1,3) + T_sim*dydt_lin_3;
    y_euler_lin_C2(i+1,4) = y_0(1,4) + T_sim*dydt_lin_4;

    y_0 = [ y_euler_lin_C2(i+1,1)  y_euler_lin_C2(i+1,2)  y_euler_lin_C2(i+1,3)  y_euler_lin_C2(i+1,4)];
end
%% Energia Mecânica Runge-Kutta C1 não linearizado

% Energia Cinética:

% Energia cinética total
K_runge_n_lin_C1 = zeros(q,1);
p=1;
while p <= q
    K_runge_n_lin_C1(p,1) = (1/6)*((m1+3*m2)*(L1^2)*((y_runge_n_lin_C1(p,3))^2)+(3*m2*cos((y_runge_n_lin_C1(p,1)-y_runge_n_lin_C1(p,2)))*L1*L2*(y_runge_n_lin_C1(p,3))*(y_runge_n_lin_C1(p,4)))+(m2*(L2^2)*((y_runge_n_lin_C1(p,4))^2)));
    p = p+1;
end

% Energia potencial

V_runge_n_lin_C1 = zeros(q,1);
p=1;

while p <= q
    V_runge_n_lin_C1(p,1) = -(1/2)*g*(((m1+m2)*cos(y_runge_n_lin_C1(p,1))*L1)+(m2*cos(y_runge_n_lin_C1(p,2))*L2));
    p = p+1;
end

% Energia mecânica 

E_runge_n_lin_C1 = zeros(q,1);
p = 1;

while p <= q
    E_runge_n_lin_C1(p,1) = K_runge_n_lin_C1(p,1)+V_runge_n_lin_C1(p,1);
    p = p+1;
end

%% Energia Mecânica Runge-Kutta C2 não linearizado

% Energia Cinética

% Energia cinética total
K_runge_n_lin_C2 = zeros(q,1);
p=1;
while p <= q
    K_runge_n_lin_C2(p,1) = (1/6)*((m1+3*m2)*(L1^2)*((y_runge_n_lin_C2(p,3))^2)+(3*m2*cos((y_runge_n_lin_C2(p,1)-y_runge_n_lin_C2(p,2)))*L1*L2*(y_runge_n_lin_C2(p,3))*(y_runge_n_lin_C2(p,4)))+(m2*(L2^2)*((y_runge_n_lin_C2(p,4))^2)));
    p = p+1;
end

% Energia potencial

V_runge_n_lin_C2 = zeros(q,1);
p=1;

while p <= q
    V_runge_n_lin_C2(p,1) = -(1/2)*g*(((m1+m2)*cos(y_runge_n_lin_C2(p,1))*L1)+(m2*cos(y_runge_n_lin_C2(p,2))*L2));
    p = p+1;
end

% Energia mecânica 

E_runge_n_lin_C2 = zeros(q,1);
p = 1;

while p <= q
    E_runge_n_lin_C2(p,1) = K_runge_n_lin_C2(p,1)+V_runge_n_lin_C2(p,1);
    p = p+1;
end

%% Energia Mecânica Runge-Kutta C1 linearizado

% Energia Cinética

% Energia cinética total
K_runge_lin_C1 = zeros(q,1);
p=1;
while p <= q
    K_runge_lin_C1(p,1) = (1/6)*((m1+3*m2)*(L1^2)*((y_runge_lin_C1(p,3))^2)+(3*m2*cos((y_runge_lin_C1(p,1)-y_runge_lin_C1(p,2)))*L1*L2*(y_runge_lin_C1(p,3))*(y_runge_lin_C1(p,4)))+(m2*(L2^2)*((y_runge_lin_C1(p,4))^2)));
    p = p+1;
end

% Energia potencial

V_runge_lin_C1 = zeros(q,1);
p=1;

while p <= q
    V_runge_lin_C1(p,1) = -(1/2)*g*(((m1+m2)*cos(y_runge_lin_C1(p,1))*L1)+(m2*cos(y_runge_lin_C1(p,2))*L2));
    p = p+1;
end

% Energia mecânica 

E_runge_lin_C1 = zeros(q,1);
p = 1;

while p <= q
    E_runge_lin_C1(p,1) = K_runge_lin_C1(p,1)+V_runge_lin_C1(p,1);
    p = p+1;
end

%% Energia Mecânica Runge-Kutta C2 linearizado

% Energia Cinética

% Energia cinética total
K_runge_lin_C2 = zeros(q,1);
p=1;
while p <= q
    K_runge_lin_C2(p,1) = (1/6)*((m1+3*m2)*(L1^2)*((y_runge_lin_C2(p,3))^2)+(3*m2*cos((y_runge_lin_C2(p,1)-y_runge_lin_C2(p,2)))*L1*L2*(y_runge_lin_C2(p,3))*(y_runge_lin_C2(p,4)))+(m2*(L2^2)*((y_runge_lin_C2(p,4))^2)));
    p = p+1;
end

% Energia potencial

V_runge_lin_C2 = zeros(q,1);
p=1;

while p <= q
    V_runge_lin_C2(p,1) = -(1/2)*g*(((m1+m2)*cos(y_runge_lin_C2(p,1))*L1)+(m2*cos(y_runge_lin_C2(p,2))*L2));
    p = p+1;
end

% Energia mecânica 

E_runge_lin_C2 = zeros(q,1);
p = 1;

while p <= q
    E_runge_lin_C2(p,1) = K_runge_lin_C2(p,1)+V_runge_lin_C2(p,1);
    p = p+1;
end

%% Energia Mecânica Euler C1 não linearizado

% Energia Cinética:

% Energia cinética total
K_euler_n_lin_C1 = zeros(q,1);
p=1;
while p <= q
    K_euler_n_lin_C1(p,1) = (1/6)*((m1+3*m2)*(L1^2)*((y_euler_n_lin_C1(p,3))^2)+(3*m2*cos((y_euler_n_lin_C1(p,1)-y_euler_n_lin_C1(p,2)))*L1*L2*(y_euler_n_lin_C1(p,3))*(y_euler_n_lin_C1(p,4)))+(m2*(L2^2)*((y_euler_n_lin_C1(p,4))^2)));
    p = p+1;
end

% Energia potencial

V_euler_n_lin_C1 = zeros(q,1);
p=1;

while p <= q
    V_euler_n_lin_C1(p,1) = -(1/2)*g*(((m1+m2)*cos(y_euler_n_lin_C1(p,1))*L1)+(m2*cos(y_euler_n_lin_C1(p,2))*L2));
    p = p+1;
end

% Energia mecânica 

E_euler_n_lin_C1 = zeros(q,1);
p = 1;

while p <= q
    E_euler_n_lin_C1(p,1) = K_euler_n_lin_C1(p,1)+V_euler_n_lin_C1(p,1);
    p = p+1;
end

%% Energia Mecânica Euler C2 não linearizado

% Energia Cinética:

% Energia cinética total
K_euler_n_lin_C2 = zeros(q,1);
p=1;
while p <= q
    K_euler_n_lin_C2(p,1) = (1/6)*((m1+3*m2)*(L1^2)*((y_euler_n_lin_C2(p,3))^2)+(3*m2*cos((y_euler_n_lin_C2(p,1)-y_euler_n_lin_C2(p,2)))*L1*L2*(y_euler_n_lin_C2(p,3))*(y_euler_n_lin_C2(p,4)))+(m2*(L2^2)*((y_euler_n_lin_C2(p,4))^2)));
    p = p+1;
end

% Energia potencial

V_euler_n_lin_C2 = zeros(q,1);
p=1;

while p <= q
    V_euler_n_lin_C2(p,1) = -(1/2)*g*(((m1+m2)*cos(y_euler_n_lin_C2(p,1))*L1)+(m2*cos(y_euler_n_lin_C2(p,2))*L2));
    p = p+1;
end

% Energia mecânica 

E_euler_n_lin_C2 = zeros(q,1);
p = 1;

while p <= q
    E_euler_n_lin_C2(p,1) = K_euler_n_lin_C2(p,1)+V_euler_n_lin_C2(p,1);
    p = p+1;
end

%% Energia Mecânica Euler C1 linearizado

% Energia Cinética:

% Energia cinética total
K_euler_lin_C1 = zeros(q,1);
p=1;
while p <= q
    K_euler_lin_C1(p,1) = (1/6)*((m1+3*m2)*(L1^2)*((y_euler_lin_C1(p,3))^2)+(3*m2*cos((y_euler_lin_C1(p,1)-y_euler_lin_C1(p,2)))*L1*L2*(y_euler_lin_C1(p,3))*(y_euler_lin_C1(p,4)))+(m2*(L2^2)*((y_euler_lin_C1(p,4))^2)));
    p = p+1;
end

% Energia potencial

V_euler_lin_C1 = zeros(q,1);
p=1;

while p <= q
    V_euler_lin_C1(p,1) = -(1/2)*g*(((m1+m2)*cos(y_euler_lin_C1(p,1))*L1)+(m2*cos(y_euler_lin_C1(p,2))*L2));
    p = p+1;
end

% Energia mecânica 

E_euler_lin_C1 = zeros(q,1);
p = 1;

while p <= q
    E_euler_lin_C1(p,1) = K_euler_lin_C1(p,1)+V_euler_lin_C1(p,1);
    p = p+1;
end

%% Energia Mecânica Euler C2 linearizado

% Energia Cinética

% Energia cinética total
K_euler_lin_C2 = zeros(q,1);
p=1;
while p <= q
    K_euler_lin_C2(p,1) = (1/6)*((m1+3*m2)*(L1^2)*((y_euler_lin_C2(p,3))^2)+(3*m2*cos((y_euler_lin_C2(p,1)-y_euler_lin_C2(p,2)))*L1*L2*(y_euler_lin_C2(p,3))*(y_euler_lin_C2(p,4)))+(m2*(L2^2)*((y_euler_lin_C2(p,4))^2)));
    p = p+1;
end

% Energia potencial

V_euler_lin_C2 = zeros(q,1);
p=1;

while p <= q
    V_euler_lin_C2(p,1) = -(1/2)*g*(((m1+m2)*cos(y_euler_lin_C2(p,1))*L1)+(m2*cos(y_euler_lin_C2(p,2))*L2));
    p = p+1;
end

% Energia mecânica 

E_euler_lin_C2 = zeros(q,1);
p = 1;

while p <= q
    E_euler_lin_C2(p,1) = K_euler_lin_C2(p,1)+V_euler_lin_C2(p,1);
    p = p+1;
end

%% Plot dos gráficos

%gráficos buter
figure(4)
plot(tempo, K_runge_n_lin_C1(:,1),"b")
hold on
plot(tempo, K_runge_lin_C1(:,1),"m")
% hold on
% plot(tempo, K_euler_n_lin_C1(:,1),"r")
% hold on
% plot(tempo, K_euler_lin_C1(:,1),"g")
% legend("Não Linear por Runge-Kutta (4,5)", "Linear por Runge-Kutta (4,5)", "Não Linear por Euler-Explícito", "Linear por Euler-Explícito")
title("Energia Cinética para o Caso 1")

figure(5)
plot(tempo, V_runge_n_lin_C1(:,1),"b")
hold on
plot(tempo, V_runge_lin_C1(:,1),"m")
% hold on
% plot(tempo, V_euler_n_lin_C1(:,1),"r")
% hold on
% plot(tempo, V_euler_lin_C1(:,1),"g")
% legend("Não Linear por Runge-Kutta (4,5)", "Linear por Runge-Kutta (4,5)", "Não Linear por Euler-Explícito", "Linear por Euler-Explícito")
title("Energia Potencial para o Caso 1")

figure(6)
plot(tempo, E_runge_n_lin_C1(:,1),"b")
hold on
plot(tempo, E_runge_lin_C1(:,1),"m")
% hold on
% plot(tempo, E_euler_n_lin_C1(:,1),"r")
% hold on
% plot(tempo, E_euler_lin_C1(:,1),"g")
% legend("Não Linear por Runge-Kutta (4,5)", "Linear por Runge-Kutta (4,5)", "Não Linear por Euler-Explícito", "Linear por Euler-Explícito")
title("Energia Mecânica para o Caso 1")

% Plotagens para o item h
% 
% figure(7)
% plot(tempo, K_runge_n_lin_C2(:,1),"b")
% hold on
% plot(tempo, K_euler_n_lin_C2(:,1),"m")
% legend("Runge-Kutta (4,5) - M1", "Euler-Explícito - M2")
% title("Energia Cinética para o Caso 2 por M1 e M2")
% 
% figure(8)
% plot(tempo, V_runge_n_lin_C2(:,1),"b")
% hold on
% plot(tempo, V_euler_n_lin_C2(:,1),"m")
% legend("Runge-Kutta (4,5) - M1", "Euler-Explícito - M2")
% title("Energia Potencial para o Caso 2 por M1 e M2")
% 
% figure(9)
% plot(tempo, E_runge_n_lin_C2(:,1),"b")
% hold on
% plot(tempo, E_euler_n_lin_C2(:,1),"m")
% legend("Runge-Kutta (4,5) - M1", "Euler-Explícito - M2")
% title("Energia Mecânica para o Caso 2 por M1 e M2")
% 
% figure(10)
% plot(tempo, y_runge_n_lin_C1(:,1),"g")
% hold on
% plot(tempo, y_euler_n_lin_C1(:,1),"r")
% legend("Runge-Kutta", "Euler Explícito")
% title("Comparação entre os métodos de integração de Euler Explícito e Runge-Kutta para o caso não linearizado")
% 
% figure(11)
% plot(tempo, y_runge_lin_C1(:,1),"g")
% hold on
% plot(tempo, y_euler_lin_C1(:,1),"r")
% legend("Runge-Kutta", "Euler Explícito")
% title("Comparação entre os métodos de integração de Euler Explícito e Runge-Kutta para o caso linearizado")
% 
% figure(12)
% plot(tempo, E_runge_n_lin_C2(:,1),"b")
% hold on
% plot(tempo, E_runge_lin_C2(:,1),"m")
% hold on
% plot(tempo, E_euler_n_lin_C2(:,1),"r")
% hold on
% plot(tempo, E_euler_lin_C2(:,1),"g")
% legend("Não Linear por Runge-Kutta (4,5)", "Linear por Runge-Kutta (4,5)", "Não Linear por Euler-Explícito", "Linear por Euler-Explícito")
% title("Energia Mecânica para o Caso 2")
% 
% figure(13)
% plot(tempo, y_runge_n_lin_C1(:,2),"b")
% hold on
% plot(tempo, y_runge_lin_C1(:,2),"m")
% ylim([-1,1])
% legend("Não Linear por Runge-Kutta (4,5)", "Linear por Runge-Kutta (4,5)")
% title("Posição barra 2 Caso 1 Runge-Kutta")
% 
% figure(14)
% plot(tempo, y_euler_n_lin_C1(:,2),"r")
% hold on
% plot(tempo, y_euler_lin_C1(:,2),"g")
% ylim([-1,1])
% legend("Não Linear por Euler Explícito", "Linear por Euler Explícito")
% title("Posição barra 2 Caso 1 Euler")
% 
% figure(15)
% plot(tempo, y_runge_n_lin_C2(:,2),"b")
% hold on
% plot(tempo, y_runge_lin_C2(:,2),"m")
% %ylim([-1,1])
% legend("Não Linear por Runge-Kutta (4,5)", "Linear por Runge-Kutta (4,5)")
% title("Posição barra 2 Caso 2 Runge-Kutta")
% 
% figure(16)
% plot(tempo, y_euler_n_lin_C2(:,2),"r")
% hold on
% plot(tempo, y_euler_lin_C2(:,2),"g")
% % ylim([-1,1])
% legend("Não Linear por Euler Explícito", "Linear por Euler Explícito")
% title("Posição barra 2 Caso 2 Euler")

% Item g gráficos
% figure(1)
% plot(tempo, y_runge_n_lin_C1(:,2),"b")
% hold on
% plot(tempo, y_runge_lin_C1(:,2),"r")
% legend("Não Linear", "Linear")
% title("Posição barra 2 no cenário 1 por Runge-Kutta (4-5)")
% 
% figure(2)
% plot(tempo, y_euler_n_lin_C1(:,2),"g")
% hold on
% plot(tempo, y_euler_lin_C1(:,2),"m")
% legend("Não Linear", "Linear")
% title("Posição barra 2 no cenário 1 por Euler Explícito")
% 
% figure(3)
% plot(tempo, y_runge_n_lin_C1(:,2),"b")
% hold on
% plot(tempo, y_euler_n_lin_C1(:,2),"g")
% legend("Runge-Kutta (4-5)", "Euler Explícito")
% title("Comparação dos métodos de integração no cenário 1 não linearizado")
% 
% figure(4)
% plot(tempo, y_runge_lin_C1(:,2),"r")
% hold on
% plot(tempo, y_euler_lin_C1(:,2),"m")
% legend("Runge-Kutta (4-5)", "Euler Explícito")
% title("Comparação dos métodos de integração no cenário 1 linearizado")
% 
% figure(5)
% plot(tempo, y_runge_n_lin_C1(:,2),"b")
% hold on
% plot(tempo, y_runge_lin_C1(:,2),"r")
% hold on
% plot(tempo, y_euler_n_lin_C1(:,2),"g")
% hold on
% plot(tempo, y_euler_lin_C1(:,2),"m")
% xlim([0,2])
% legend("Não Linear - Runge Kutta", "Não Linear - Euler Explícito", "Linear - Runge Kutta", "Linear - Euler Explícito")
% title("Posição barra 2 no cenário 1 por Runge-Kutta (4-5) e Euler Explícito")
%% Defininido os espaçoes de estados

%Não linearizado
function dydt_n_lin = f_n_lin(t, y_0)
wp = 1;
g = 9.8;
mi = 2;
L1 = g/(wp^2);
L2 = L1*20/21;
m1 = L1*mi;
m2 = L2*mi;
dydt_n_lin_1 = y_0(3);
dydt_n_lin_2 = y_0(4);
dydt_n_lin_3 = (-3*g*((4*m1 + 5*m2)*sin(y_0(1)) + 3*m2*sin(y_0(1) - 2*y_0(2))))/((8*m1 + 15*m2 - 9*m2*cos(2*(y_0(1) - y_0(2))))*L1) + (9*m2*sin(2*(y_0(1) - y_0(2)))*(y_0(3)*y_0(3)))/(-8*m1 - 15*m2 + 9*m2*cos(2*(y_0(1) - y_0(2)))) + (6*m2*sin(y_0(1) - y_0(2))*L2*(y_0(4)*y_0(4)))/((-4*(m1 + 3*m2) + 9*m2*(cos(y_0(1) - y_0(2))*cos(y_0(1) - y_0(2))))*L1);
dydt_n_lin_4 = (9*g*(m1 + 2*m2)*sin(2*y_0(1) - y_0(2)) - 3*g*(m1 + 6*m2)*sin(y_0(2)))/((8*m1 + 15*m2 - 9*m2*cos(2*(y_0(1) - y_0(2))))*L2) + (6*(m1 + 3*m2)*sin(y_0(1) - y_0(2))*L1*(y_0(3)*y_0(3)))/((4*(m1 + 3*m2) - 9*m2*(cos(y_0(1) - y_0(2))*cos(y_0(1) - y_0(2))))*L2) + (9*m2*sin(2*(y_0(1) - y_0(2)))*(y_0(4)*y_0(4)))/(8*m1 + 15*m2 - 9*m2*cos(2*(y_0(1) - y_0(2))));
dydt_n_lin =  [dydt_n_lin_1; dydt_n_lin_2; dydt_n_lin_3; dydt_n_lin_4];
end

% Linearizado
function dydt_lin = f_lin(t, y_0)
wp = 1;
g = 9.8;
mi = 2;
L1 = g/(wp^2);
L2 = L1*20/21;
L = L1/L2;
m1 = L1*mi;
m2 = L2*mi;
dydt_lin_1 = y_0(3);
dydt_lin_2 = y_0(4);
dydt_lin_3 = (3*(wp^2)*((-2*(2+L)*y_0(1))+(3*y_0(2))))/(3+4*L);
dydt_lin_4 = (3*L*(wp^2)*(((3*(2+L))*y_0(1))-(2*(3+L)*y_0(2))))/(3+4*L);
dydt_lin =  [dydt_lin_1; dydt_lin_2; dydt_lin_3; dydt_lin_4];
end

