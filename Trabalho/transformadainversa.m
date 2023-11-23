syms s
F=(0.00108284/(-0.193739+s))+(0.00108314/(0.193739+s))+((-0.813703-0.00355412*s)/(72.0567+0.36599*s+s^2))+((0.948302+0.00355382*s)/(8413.84+38.1225*s+s^2));
%F ACIMA TOMADA COMO A DECOMPOSIÇÃO EM FRAÇÃO PARCIAL DO TERMO (2,1) DA
%FUNÇÃO DE TRANSFERÊNCIA, OU SEJA, A RESPOSTA À IMPULSO DE POSIÇÃO NA
%COORDENADA THETA
f=ilaplace(F);

t_p = 0:0.1:30;

FUN = matlabFunction(f);
trans = feval(FUN, t_p);

T1 = (4995101593999441.*exp(-(1745045776414265.*t_p)/9007199254740992))/4611686018427387904;
T2 = (4993718088193913.*exp((1745045776414265.*t_p)/9007199254740992))/4611686018427387904;
T3 = (4097621377953287.*exp(-(36599.*t_p)/200000).*(cos((2^(1/2).*24746987492607220882933^(1/2).*t_p)/26214400000) + (722735512897644203980292096.*2^(1/2).*24746987492607220882933^(1/2).*sin((2^(1/2).*24746987492607220882933^(1/2).*t_p)/26214400000))/5964928528802938777666596178386444163))/1152921504606846976; 
T4 = (4097275501501905.*exp(-(15249.*t_p)/800).*(cos((1383066871628038169^(1/2).*t_p)/13107200) + (2661334583892191449759744.*1383066871628038169^(1/2).*sin((1383066871628038169^(1/2).*t_p)/13107200))/1133361202012088190548213633242389))/1152921504606846976;

figure(1)
plot(t_p, trans, "b")
grid on
grid minor
title("Transformada Inversa de Lapace - Impulso")
xlabel('Tempo (s)')
ylabel('Amplitude de Theta (radianos)')


figure(2)
plot(t_p, T1, "r")
hold on
plot(t_p, T2, "g")
hold on
plot(t_p, T3, "b")
hold on
plot(t_p, T4, "y")
grid on
grid minor
legend( "T1", "T2", "T3", "T4")
title("Transformada Inversa de Lapace - Impulso")
xlabel('Tempo (s)')
ylabel('Amplitude de Theta (radianos)')

figure(3)
plot(t_p, T2, "g")
hold on
plot(t_p, trans, "y")
grid on
grid minor
legend("T2", "Resposta do sistema")
title("Transformada Inversa de Lapace - Impulso")
xlabel('Tempo (s)')
ylabel('Amplitude de Theta (radianos)')