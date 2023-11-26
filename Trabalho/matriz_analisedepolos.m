%%Análise de polos para variação de velocidade  longitudinal 
resp = [];
intervalovel = 0:0.1:120;

for n = intervalovel
    u = n;
    a = 13*3.14/180;
    v = 0;
    %primeira coluna
    p11 = 0;
    p21 = 0;
    p31 = 0;
    p41 = -125434/15;
    p51 = -(2.20134*10^9/(-1.68644*10^7 + 425920*cos(a)^2));
    p61 = (5.05419*10^7*sec(a)/(851840 - 3.37288*10^7*sec(a)^2));
    
    %segunda coluna
    p12 = 0;
    p22 = 0;
    p32 = 0;
    p42 = 54434/15;
    p52 = (2.20134*(10^9))/(-1.68644*(10^7)+425920*(cos(a)^2));
    p62 = (-5.05419*(10^7)*sec(a))/(851840-3.37288*(10^7)*(sec(a)^2));
    
    %terceira coluna
    p13 = 0;
    p23 = 0;
    p33 = 0;
    p43 = 0;
    p53 = (688.225*(0.15714+0.01111*a+0.00017*a^2)*(u+v)^2*cos(a)^2)/(-1.68644*10^7+425920*cos(a)^2)+(2.42*(0.0041+0.000041*u)*(1.72656*10^6-125.132*(0.21999+0.07378*a-0.00165*a^2)*(u+v)^2)*cos(a)^2)/(-1.68644*10^7+425920*cos(a)^2)+(4.17828*10^6*cos(a)*sin(a))/(-1.68644*10^7+425920*cos(a)^2)-(688.225*(0.21999+0.07378*a-0.00165*a^2)*(u+v)^2*cos(a)*sin(a))/(-1.68644*10^7+425920*cos(a)^2);
    p63 = -((3.79843*10^6*(0.0041+0.000041*u)*sec(a))/(851840-3.37288*10^7*sec(a)^2))-(625.659*(0.15714+0.01111*a+0.00017*a^2)*(u+v)^2*sec(a))/(851840-3.37288*10^7*sec(a)^2)+(275.29*(0.21999+0.07378*a-0.00165*a^2)*(0.0041+0.000041*u)*(u+v)^2*sec(a))/(851840-3.37288*10^7*sec(a)^2)-(3.79843*10^6*sec(a)*tan(a))/(851840-3.37288*10^7*sec(a)^2)+(625.659*(0.21999+0.07378*a-0.00165*a^2)*(u+v)^2*sec(a)*tan(a))/(851840-3.37288*10^7*sec(a)^2);
    
    %quarta coluna 
    p14 = 1;
    p24 = 0;
    p34 = 0;
    p44 = -(13987/375);
    p54 = 1.9585*10^7/(1.68644*10^7-425920*cos(a)^2);
    p64 = (224831*sec(a))/(425920-1.68644*10^7*sec(a)^2);
    
    %quinta coluna
    p15 = 0;
    p25 = 1;
    p35 = 0;
    p45 = 25549/750;
    p55 = 1.72348*10^12/(-1.48407*10^12+3.7481*10^10*cos(a)^2);
    p65 = (224831*sec(a))/(-425920+1.68644*10^7*sec(a)^2);
    
    %sexta coluna
    p16 = 0;
    p26 = 0;
    p36 = 1;
    p46 = 0;
    p56 = 0;
    p66 = 0;
    
    Polos = [p11 p21 p31 p41 p51 p61; p12 p22 p32 p42 p52 p62; p13 p23 p33 p43 p53 p63; p14 p24 p34 p44 p54 p64; p15 p25 p35 p45 p55 p65; p16 p26 p36 p46 p56 p66];
    z = eig(Polos);
    resp = cat(2, resp, z);
end

resp_rel = resp(3:6,:);
resp_imag = imag(resp_rel);
resp_real = real(resp);
relf = abs(resp_rel);
for i=1:length(relf)
    relf(4,i) = -relf(4,i);
end

figure(1)
plot(intervalovel, resp_real(3:4,:), "b")
hold on
plot(intervalovel, relf(3,:), "r")
hold on
plot(intervalovel, relf(4,:), "c")
grid on
grid minor
legend( "", "Parte real do complexo conjugado", "Real positivo", "Real negativo")
title("Polos relevantes do sistema")
xlabel('Velocidade longitudinal (m/s)')


%%%%%%%%%%%%%%%%%%%%%%Análise de polos para variação de velocidade do vento
resp = [];
intervalovento = 0:0.1:50;

for n = intervalovento
    u = 83.35;
    a = 13*pi/180;
    v = n;
    %primeira coluna
    p11 = 0;
    p21 = 0;
    p31 = 0;
    p41 = -125434/15;
    p51 = -(2.20134*10^9/(-1.68644*10^7 + 425920*cos(a)^2));
    p61 = (5.05419*10^7*sec(a)/(851840 - 3.37288*10^7*sec(a)^2));
    
    %segunda coluna
    p12 = 0;
    p22 = 0;
    p32 = 0;
    p42 = 54434/15;
    p52 = (2.20134*(10^9))/(-1.68644*(10^7)+425920*(cos(a)^2));
    p62 = (-5.05419*(10^7)*sec(a))/(851840-3.37288*(10^7)*(sec(a)^2));
    
    %terceira coluna
    p13 = 0;
    p23 = 0;
    p33 = 0;
    p43 = 0;
    p53 = (688.225*(0.15714+0.01111*a+0.00017*a^2)*(u+v)^2*cos(a)^2)/(-1.68644*10^7+425920*cos(a)^2)+(2.42*(0.0041+0.000041*u)*(1.72656*10^6-125.132*(0.21999+0.07378*a-0.00165*a^2)*(u+v)^2)*cos(a)^2)/(-1.68644*10^7+425920*cos(a)^2)+(4.17828*10^6*cos(a)*sin(a))/(-1.68644*10^7+425920*cos(a)^2)-(688.225*(0.21999+0.07378*a-0.00165*a^2)*(u+v)^2*cos(a)*sin(a))/(-1.68644*10^7+425920*cos(a)^2);
    p63 = -((3.79843*10^6*(0.0041+0.000041*u)*sec(a))/(851840-3.37288*10^7*sec(a)^2))-(625.659*(0.15714+0.01111*a+0.00017*a^2)*(u+v)^2*sec(a))/(851840-3.37288*10^7*sec(a)^2)+(275.29*(0.21999+0.07378*a-0.00165*a^2)*(0.0041+0.000041*u)*(u+v)^2*sec(a))/(851840-3.37288*10^7*sec(a)^2)-(3.79843*10^6*sec(a)*tan(a))/(851840-3.37288*10^7*sec(a)^2)+(625.659*(0.21999+0.07378*a-0.00165*a^2)*(u+v)^2*sec(a)*tan(a))/(851840-3.37288*10^7*sec(a)^2);
    
    %quarta coluna 
    p14 = 1;
    p24 = 0;
    p34 = 0;
    p44 = -(13987/375);
    p54 = 1.9585*10^7/(1.68644*10^7-425920*cos(a)^2);
    p64 = (224831*sec(a))/(425920-1.68644*10^7*sec(a)^2);
    
    %quinta coluna
    p15 = 0;
    p25 = 1;
    p35 = 0;
    p45 = 25549/750;
    p55 = 1.72348*10^12/(-1.48407*10^12+3.7481*10^10*cos(a)^2);
    p65 = (224831*sec(a))/(-425920+1.68644*10^7*sec(a)^2);
    
    %sexta coluna
    p16 = 0;
    p26 = 0;
    p36 = 1;
    p46 = 0;
    p56 = 0;
    p66 = 0;
    
    Polos = [p11 p21 p31 p41 p51 p61; p12 p22 p32 p42 p52 p62; p13 p23 p33 p43 p53 p63; p14 p24 p34 p44 p54 p64; p15 p25 p35 p45 p55 p65; p16 p26 p36 p46 p56 p66];
    z = eig(Polos);
    resp = cat(2, resp, z);
end

resp_rel = resp(3:6,:);
resp_imag = imag(resp_rel);
resp_real = real(resp);

relf = abs(resp_rel);
for i=1:length(relf)
    relf(4,i) = -relf(4,i);
end

% plot(intervalo, resp_imag(1:2,:), "b--")
% hold on
figure(2)
plot(intervalovento, resp_real(3:4,:), "b")
hold on
plot(intervalovento, relf(3,:), "r")
hold on
plot(intervalovento, relf(4,:), "c")
legend( "", "Parte real do complexo conjugado", "Real positivo", "Real negativo" )
grid on
grid minor
title("Polos relevantes do sistema para Phi " + 180*a/pi + "° e velocidade longitudinal " + u + "m/s")
xlabel('Velocidade do vento (m/s)')

%%%%%%%%%%%%%%%%%%%%%%Análise de polos para variação de ângulo
resp = [];
intervaloang = 0:0.01:pi/3;

for n = intervaloang
    u = 83.35;
    a = n;
    v = 0;
    %primeira coluna
    p11 = 0;
    p21 = 0;
    p31 = 0;
    p41 = -125434/15;
    p51 = -(2.20134*10^9/(-1.68644*10^7 + 425920*cos(a)^2));
    p61 = (5.05419*10^7*sec(a)/(851840 - 3.37288*10^7*sec(a)^2));
    
    %segunda coluna
    p12 = 0;
    p22 = 0;
    p32 = 0;
    p42 = 54434/15;
    p52 = (2.20134*(10^9))/(-1.68644*(10^7)+425920*(cos(a)^2));
    p62 = (-5.05419*(10^7)*sec(a))/(851840-3.37288*(10^7)*(sec(a)^2));
    
    %terceira coluna
    p13 = 0;
    p23 = 0;
    p33 = 0;
    p43 = 0;
    p53 = (688.225*(0.15714+0.01111*a+0.00017*a^2)*(u+v)^2*cos(a)^2)/(-1.68644*10^7+425920*cos(a)^2)+(2.42*(0.0041+0.000041*u)*(1.72656*10^6-125.132*(0.21999+0.07378*a-0.00165*a^2)*(u+v)^2)*cos(a)^2)/(-1.68644*10^7+425920*cos(a)^2)+(4.17828*10^6*cos(a)*sin(a))/(-1.68644*10^7+425920*cos(a)^2)-(688.225*(0.21999+0.07378*a-0.00165*a^2)*(u+v)^2*cos(a)*sin(a))/(-1.68644*10^7+425920*cos(a)^2);
    p63 = -((3.79843*10^6*(0.0041+0.000041*u)*sec(a))/(851840-3.37288*10^7*sec(a)^2))-(625.659*(0.15714+0.01111*a+0.00017*a^2)*(u+v)^2*sec(a))/(851840-3.37288*10^7*sec(a)^2)+(275.29*(0.21999+0.07378*a-0.00165*a^2)*(0.0041+0.000041*u)*(u+v)^2*sec(a))/(851840-3.37288*10^7*sec(a)^2)-(3.79843*10^6*sec(a)*tan(a))/(851840-3.37288*10^7*sec(a)^2)+(625.659*(0.21999+0.07378*a-0.00165*a^2)*(u+v)^2*sec(a)*tan(a))/(851840-3.37288*10^7*sec(a)^2);
    
    %quarta coluna 
    p14 = 1;
    p24 = 0;
    p34 = 0;
    p44 = -(13987/375);
    p54 = 1.9585*10^7/(1.68644*10^7-425920*cos(a)^2);
    p64 = (224831*sec(a))/(425920-1.68644*10^7*sec(a)^2);
    
    %quinta coluna
    p15 = 0;
    p25 = 1;
    p35 = 0;
    p45 = 25549/750;
    p55 = 1.72348*10^12/(-1.48407*10^12+3.7481*10^10*cos(a)^2);
    p65 = (224831*sec(a))/(-425920+1.68644*10^7*sec(a)^2);
    
    %sexta coluna
    p16 = 0;
    p26 = 0;
    p36 = 1;
    p46 = 0;
    p56 = 0;
    p66 = 0;
    
    Polos = [p11 p21 p31 p41 p51 p61; p12 p22 p32 p42 p52 p62; p13 p23 p33 p43 p53 p63; p14 p24 p34 p44 p54 p64; p15 p25 p35 p45 p55 p65; p16 p26 p36 p46 p56 p66];
    z = eig(Polos);
    resp = cat(2, resp, z);
end

resp_rel = resp(3:6,:);
resp_imag = imag(resp_rel);
resp_real = real(resp);

relf = abs(resp_rel);
for i=1:length(relf)
    relf(4,i) = -relf(4,i);
end

% plot(intervalo, resp_imag(1:2,:), "b--")
% hold on
figure(3)
plot(intervaloang*180/pi, resp_real(3:4,:), "b")
hold on
plot(intervaloang*180/pi, relf(3,:), "r")
hold on
plot(intervaloang*180/pi, relf(4,:), "c")
grid on
grid minor
legend( "", "Parte real do complexo conjugado", "Real positivo", "Real negativo" )
title("Polos relevantes do sistema para vento " + v + "m/s e velocidade longitudinal " + u + "m/s")
xlabel('Ângulo (graus)')
