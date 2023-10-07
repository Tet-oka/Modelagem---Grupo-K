%% Parâmetros do problema

wp = 1;
g = 9.8;
L1 = g/(wp^2);
L2 = L1*20/21;

%% Definição dos parâmetros iniciais da integração - item e

%Ao substituir os valor de \lambda = (21/20) e de \omega_p = 1, teremos que
%a freqência natural no primeiro modo é de 0,4Hz. 

%O período de simulação utilizado será de 4 períodos completos de
%oscilação, que deve ser suficiente para verificar diferenças entre os
%métodos de integração e entre as equações linearizadas e não linearizadas
F_1modo = 0.4;
t = 1/F_1modo;
t = 4*t;

%Se usarmos 100 passos de integração por segundo:
F_sim = 1/100;

%criando o vetor tempo com 100 passos por segundo e 4 períodos de oscilação
tempo = 0:F_sim:t;

%condições iniciais de integração:

%Cenário 1:
y_0_1 = [3*pi/180 7*pi/180 0 0];

%Cenário 2
y_0_2 = [3*pi/180 7*pi/180 0.1 -0.1];


%% Aplicação dos métodos de integração

%Runge-Kutta 4-5 não linearizado
[t_runge_n_lin_C1, y_runge_n_lin_C1] = ode45(@f_n_lin_C1, tempo, y_0_1);

[t_runge_n_lin_C2, y_runge_n_lin_C2] = ode45(@f_n_lin_C2, tempo, y_0_2);

plot(t_runge_n_lin_C1, y_runge_n_lin_C1(:,1),"g")
hold on
plot(t_runge_n_lin_C2, y_runge_n_lin_C2(:,1),"r")
%% Defininido os espaçoes de estados

%Não linearizado C1
function dydt_n_lin_C1 = f_n_lin_C1(t, y_0_1)
wp = 1;
g = 9.8;
L1 = g/(wp^2);
L2 = L1*20/21;
dydt_n_lin_1 = y_0_1(3);
dydt_n_lin_2 = y_0_1(4);
dydt_n_lin_3 = (-3*g*((4*L1 + 5*L2)*sin(y_0_1(1)) + 3*L2*sin(y_0_1(1) - 2*y_0_1(2))))/((8*L1 + 15*L2 - 9*L2*cos(2*(y_0_1(1) - y_0_1(2))))*L1) + (9*L2*sin(2*(y_0_1(1) - y_0_1(2)))*(y_0_1(3)*y_0_1(3)))/(-8*L1 - 15*L2 + 9*L2*cos(2*(y_0_1(1) - y_0_1(2)))) + (6*L2*sin(y_0_1(1) - y_0_1(2))*L2*(y_0_1(4)*y_0_1(4)))/((-4*(L1 + 3*L2) + 9*L2*(cos(y_0_1(1) - y_0_1(2))*cos(y_0_1(1) - y_0_1(2))))*L1);
dydt_n_lin_4 = (9*g*(L1 + 2*L2)*sin(2*y_0_1(1) - y_0_1(2)) - 3*g*(L1 + 6*L2)*sin(y_0_1(2)))/((8*L1 + 15*L2 - 9*L2*cos(2*(y_0_1(1) - y_0_1(2))))*L2) + (6*(L1 + 3*L2)*sin(y_0_1(1) - y_0_1(2))*L1*(y_0_1(3)*y_0_1(3)))/((4*(L1 + 3*L2) - 9*L2*(cos(y_0_1(1) - y_0_1(2))*cos(y_0_1(1) - y_0_1(2))))*L2) + (9*L2*sin(2*(y_0_1(1) - y_0_1(2)))*(y_0_1(4)*y_0_1(4)))/(8*L1 + 15*L2 - 9*L2*cos(2*(y_0_1(1) - y_0_1(2))));
dydt_n_lin_C1 =  [dydt_n_lin_1; dydt_n_lin_2; dydt_n_lin_3; dydt_n_lin_4];
end

%Não linearizado C2
function dydt_n_lin_C2 = f_n_lin_C2(t, y_0_2)
wp = 1;
g = 9.8;
L1 = g/(wp^2);
L2 = L1*20/21;
dydt_n_lin_1 = y_0_2(3);
dydt_n_lin_2 = y_0_2(4);
dydt_n_lin_3 = (-3*g*((4*L1 + 5*L2)*sin(y_0_2(1)) + 3*L2*sin(y_0_2(1) - 2*y_0_2(2))))/((8*L1 + 15*L2 - 9*L2*cos(2*(y_0_2(1) - y_0_2(2))))*L1) + (9*L2*sin(2*(y_0_2(1) - y_0_2(2)))*(y_0_2(3)*y_0_2(3)))/(-8*L1 - 15*L2 + 9*L2*cos(2*(y_0_2(1) - y_0_2(2)))) + (6*L2*sin(y_0_2(1) - y_0_2(2))*L2*(y_0_2(4)*y_0_2(4)))/((-4*(L1 + 3*L2) + 9*L2*(cos(y_0_2(1) - y_0_2(2))*cos(y_0_2(1) - y_0_2(2))))*L1);
dydt_n_lin_4 = (9*g*(L1 + 2*L2)*sin(2*y_0_2(1) - y_0_2(2)) - 3*g*(L1 + 6*L2)*sin(y_0_2(2)))/((8*L1 + 15*L2 - 9*L2*cos(2*(y_0_2(1) - y_0_2(2))))*L2) + (6*(L1 + 3*L2)*sin(y_0_2(1) - y_0_2(2))*L1*(y_0_2(3)*y_0_2(3)))/((4*(L1 + 3*L2) - 9*L2*(cos(y_0_2(1) - y_0_2(2))*cos(y_0_2(1) - y_0_2(2))))*L2) + (9*L2*sin(2*(y_0_2(1) - y_0_2(2)))*(y_0_2(4)*y_0_2(4)))/(8*L1 + 15*L2 - 9*L2*cos(2*(y_0_2(1) - y_0_2(2))));
dydt_n_lin_C2 =  [dydt_n_lin_1; dydt_n_lin_2; dydt_n_lin_3; dydt_n_lin_4];
end
