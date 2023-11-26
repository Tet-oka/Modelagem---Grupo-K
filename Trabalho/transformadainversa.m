%%F ACIMA TOMADA COMO A DECOMPOSIÇÃO EM FRAÇÃO PARCIAL DO TERMO (2,1) DA FUNÇÃO DE TRANSFERÊNCIA, OU SEJA, A RESPOSTA À IMPULSO DE POSIÇÃO NA COORDENADA THETA
syms s
Ft=(-0.00108284/(-0.193739+s))+(0.00108314/(0.193739+s))+((-0.813703-0.00355412*s)/(72.0567+0.36599*s+s^2))+((0.948302+0.00355382*s)/(8413.84+38.1225*s+s^2));
fq=ilaplace(Ft);

t_t = 0:0.1:30;

FUNt = matlabFunction(fq);
transt = feval(FUNt, t_t);

T1 = (4995101593999441.*exp(-(1745045776414265.*t_t)/9007199254740992))/4611686018427387904;
T2 = (-4993718088193913.*exp((1745045776414265.*t_t)/9007199254740992))/4611686018427387904;
T3 = (4097621377953287.*exp(-(36599.*t_t)/200000).*(cos((2^(1/2).*24746987492607220882933^(1/2).*t_t)/26214400000) + (722735512897644203980292096.*2^(1/2).*24746987492607220882933^(1/2).*sin((2^(1/2).*24746987492607220882933^(1/2).*t_t)/26214400000))/5964928528802938777666596178386444163))/1152921504606846976; 
T4 = (4097275501501905.*exp(-(15249.*t_t)/800).*(cos((1383066871628038169^(1/2).*t_t)/13107200) + (2661334583892191449759744.*1383066871628038169^(1/2).*sin((1383066871628038169^(1/2).*t_t)/13107200))/1133361202012088190548213633242389))/1152921504606846976;

figure(1)
plot(t_t, transt, "b")
grid on
grid minor
title("Transformada inversa de Laplace para theta- Impulso")
xlabel('Tempo (s)')
ylabel('Amplitude de Theta (radianos)')


figure(2)
plot(t_t, T1, "r")
hold on
plot(t_t, T2, "g")
hold on
plot(t_t, T3, "b")
hold on
plot(t_t, T4, "y")
grid on
grid minor
legend( "Termo com exp(-0,1937) - T1", "Termo com exp(+0,1937) - T2", "Termo com exp(-0.183) - T3", "Termo com exp(-19,061) - T4")
title("Termos da transformada inversa de Laplace para theta - Impulso")
xlabel('Tempo (s)')
ylabel('Amplitude de Theta (radianos)')

figure(3)
plot(t_t, T2, "g")
hold on
plot(t_t, transt, "y")
grid on
grid minor
title("Região de atuação do modelo na transformada inversa de Laplace para theta - Impulso")
xlabel('Tempo (s)')
ylabel('Amplitude de Theta (radianos)')
p = patch([0 0 5 5],[-0.4 0.1 0.1 -0.4],'');
ylim([-0.4 0.1])
set(p,'FaceAlpha',0.1)
set(p,'EdgeColor','none')
legend("Termo com exponencial positiva", "Resposta do sistema", "Região passível de análise free-roll")

%%F ACIMA TOMADA COMO A DECOMPOSIÇÃO EM FRAÇÃO PARCIAL DO TERMO (2,1) DA FUNÇÃO DE TRANSFERÊNCIA, OU SEJA, A RESPOSTA À IMPULSO DE POSIÇÃO NA COORDENADA THETA
syms s
Fq=((72.7837+0.317715*s)/(72.0567+0.36599*s+1*s^2)+(-84.7795-0.317715*s)/(8413.84+38.1225*s+1*s^2));
fq=ilaplace(Fq);

t_q = 0:0.1:30;

FUNq = matlabFunction(fq);
transq = feval(FUNq, t_q);

figure(4)
plot(t_q, transq, "b")
grid on
grid minor
title("Transformada inversa de Laplace para q2- Impulso")
xlabel('Tempo (s)')
ylabel('Amplitude de q2 (metros)')

