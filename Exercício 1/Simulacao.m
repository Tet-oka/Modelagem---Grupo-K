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
t = 4*t;

% Se usarmos 100 passos de integração por segundo:
T_sim = 1/100;

% criando o vetor tempo com 100 passos por segundo e 4 períodos de oscilação
tempo = 0:T_sim:t;

% número pontos de integração
q = size(tempo(1,:));
q = q(2);

% condições iniciais de integração:

% Cenário 1:
y_0_1 = [3*pi/180 7*pi/180 0 0];

% Cenário 2
y_0_2 = [3*pi/180 7*pi/180 0.1 -0.1];


%% Aplicação de Runge-Kutta (4-5) em C1 e C2

y_0 = y_0_1;    % aplicando para o cenário 1
[t_runge_n_lin_C1, y_runge_n_lin_C1] = ode45(@f_n_lin, tempo, y_0);     % não linearizado
[t_runge_lin_C1, y_runge_lin_C1] = ode45(@f_lin, tempo, y_0);           % linearizado

y_0= y_0_2;     % aplicando para o cenário 2
[t_runge_n_lin_C2, y_runge_n_lin_C2] = ode45(@f_n_lin, tempo, y_0);     % não linearizado
[t_runge_lin_C2, y_runge_lin_C2] = ode45(@f_lin, tempo, y_0);           % linearizado



%% Para os cálculos de energia

I1 = (m1*L1^2)/12;                      % Momento de inércia
I2 = (m2*L2^2)/12;                      % Momento de inércia

% Para C1 não linearizado:
cos1_C1 = cos(y_runge_n_lin_C1(:,1));  % Vetor de cossenos de theta_1
cos2_C1= cos(y_runge_n_lin_C1(:,2));   % Vetor de cossenos de Theta_2
sin1_C1 = sin(y_runge_n_lin_C1(:,1));  % Vetor de senos de theta_1
sin2_C1= sin(y_runge_n_lin_C1(:,2));   % Vetor de senos de theta_2
w1quadrado_C1 = y_runge_n_lin_C1.^2;   % Vetor com os quadrados das velocidades angulares da barra 1
w2quadrado_C1 = y_runge_n_lin_C1.^2;   % Vetor com os quadrados das velocidades angulares da barra 2

% Para C2 não linearizado:
cos1_C2 = cos(y_runge_n_lin_C2(:,1));  % Vetor de cossenos de theta_1
cos2_C2= cos(y_runge_n_lin_C2(:,2));   % Vetor de cossenos de Theta_2
sin1_C2 = sin(y_runge_n_lin_C2(:,1));  % Vetor de senos de theta_1
sin2_C2= sin(y_runge_n_lin_C2(:,2));   % Vetor de senos de theta_2
w1quadrado_C2 = y_runge_n_lin_C2.^2;   % Vetor com os quadrados das velocidades angulares da barra 1
w2quadrado_C2 = y_runge_n_lin_C2.^2;   % Vetor com os quadrados das velocidades angulares da barra 2

%% Energia Mecânica C1 não linearizado

% Energia Cinética:

% Velocidade CG 1
V1_n_lin_C1 = zeros(q,1);
p = 1;
while p <= q
    v = y_runge_n_lin_C1(p,3)*L1/2;
    V1_n_lin_C1(p,1) = v^2;
    p = p+1;
end

% Velocidade CG 2 em x
V2x_n_lin_C1 = zeros(q,1);
p = 1;
while p <= q
    v = (y_runge_n_lin_C1(p,3)*L1*cos1_C1(p,1))+(y_runge_n_lin_C1(p,4)*(L2/2)*cos2_C1(p,1));
    V2x_n_lin_C1(p,1) = v^2;
    p = p+1;
end

% Velocidade CG 2 em y
V2y_n_lin_C1 = zeros(q,1);
p = 1;
while p <= q
    v = (y_runge_n_lin_C1(p,3)*L1*sin1_C1(p,1))+(y_runge_n_lin_C1(p,4)*(L2/2)*sin2_C1(p,1));
    V2y_n_lin_C1(p,1) = v^2;
    p = p+1;
end

% Energia cinética total
K_n_lin_C1 = zeros(q,1);
p=1;
while p <= q
    k = ((m1*V1_n_lin_C1(p,1))/2)+(m2*((V2x_n_lin_C1(p,1))+(V2y_n_lin_C1(p,1)))/2)+(I1*(w1quadrado_C1(p,3))/2)+(I2*(w2quadrado_C1(p,4))/2);
    K_n_lin_C1(p,1) = k;
    p = p+1;
end

% Energia potencial

V_n_lin_C1 = zeros(q,1);
p=1;

while p <= q
    V_n_lin_C1(p,1) = (m1*g*L1*(1-cos1_C1(p,1))/2)+(m2*g*((L1*(1-cos1_C1(p,1)))+(L2*(1-cos2_C1(p,1)/2))));
    p = p+1;
end

% Energia mecânica 

E_n_lin_C1 = zeros(q,1);
p = 1;

while p <= q
    E_n_lin_C1(p,1) = K_n_lin_C1(p,1)+V_n_lin_C1(p,1);
    p = p+1;
end

%% Energia Mecânica C2 não linearizado

% Energia Cinética

% Velocidade CG 1
V1_n_lin_C2 = zeros(q,1);
p = 1;
while p <= q
    v = y_runge_n_lin_C2(p,3)*L1/2;
    V1_n_lin_C2(p,1) = v^2;
    p = p+1;
end

% Velocidade CG 2 em x
V2x_n_lin_C2 = zeros(q,1);
p = 1;
while p <= q
    v = (y_runge_n_lin_C2(p,3)*L1*cos1_C2(p,1))+(y_runge_n_lin_C2(p,4)*(L2/2)*cos2_C2(p,1));
    V2x_n_lin_C2(p,1) = v^2;
    p = p+1;
end

% Velocidade CG 2 em y
V2y_n_lin_C2 = zeros(q,1);
p = 1;
while p <= q
    v = (y_runge_n_lin_C2(p,3)*L1*sin1_C2(p,1))+(y_runge_n_lin_C2(p,4)*(L2/2)*sin2_C2(p,1));
    V2y_n_lin_C2(p,1) = v^2;
    p = p+1;
end

% Energia cinética total
K_n_lin_C2 = zeros(q,1);
p=1;
while p <= q
    k = ((m1*V1_n_lin_C2(p,1))/2)+(m2*((V2x_n_lin_C2(p,1))+(V2y_n_lin_C2(p,1)))/2)+(I1*(w1quadrado_C2(p,3))/2)+(I2*(w2quadrado_C2(p,4))/2);
    K_n_lin_C2(p,1) = k;
    p = p+1;
end

% Energia potencial

V_n_lin_C2 = zeros(q,1);
p=1;

while p <= q
    V_n_lin_C2(p,1) = (m1*g*L1*(1-cos1_C2(p,1))/2)+(m2*g*((L1*(1-cos1_C2(p,1)))+(L2*(1-cos2_C2(p,1)/2))));
    p = p+1;
end

% Energia mecânica 

E_n_lin_C2 = zeros(q,1);
p = 1;

while p <= q
    E_n_lin_C2(p,1) = K_n_lin_C2(p,1)+V_n_lin_C2(p,1);
    p = p+1;
end
%% Plot dos gráficos

figure(1)
plot(tempo, y_runge_n_lin_C2(:,1),"g")
hold on
plot(tempo, y_runge_lin_C2(:,1),"r")

figure(2)
plot(tempo, K_n_lin_C1(:,1),"b")
hold on
plot(tempo, V_n_lin_C1(:,1),"m")

figure(3)
plot(tempo, E_n_lin_C1(:,1),"k")

figure(4)
plot(tempo, E_n_lin_C2(:,1),"m")

%% Testes


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
dydt_lin_3 = (3*(wp^2)*((-2*(1+2*L)*y_0(1))+(L*(3+2*L)*y_0(2))))/(4+3*L);
dydt_lin_4 = (((3*wp^2)*((3+(4*L*(2+L)*y_0(1)))-(2*((1+L)^3)*y_0(2))))/(L*(4+3*L)));
dydt_lin =  [dydt_lin_1; dydt_lin_2; dydt_lin_3; dydt_lin_4];
end

