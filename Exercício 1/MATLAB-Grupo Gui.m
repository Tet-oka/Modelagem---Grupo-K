tspan = 0:0.01:10;
y_0 = [0.05 0 0 0];
%Linearizada
[t1,y1] = ode45(@f1, tspan, y_0); %linearizado ode45
[t2,y2] = ode23t(@f1, tspan, y_0); %linearizado ode23t
%Não linearizada
[t3,y3] = ode45(@f2, tspan, y_0); %N linearizado ode45
[t4,y4] = ode23t(@f2, tspan, y_0); %N linearizado ode23t

%Equações linearizada e não linearizada para ode45
figure(1)
plot(t1,y1(:,1),"b")
hold on 
plot(t3,y3(:,1),"r")
legend("Linearizada", "Não Linearizada")
xlabel("tempo [s]")
ylabel("theta_1 [rad]")

%Equações linearizada e não linearizada para ode23t
figure(2) 
plot(t2,y2(:,1),"b")
hold on
plot(t4,y4(:,1),"r")
legend("Linearizada", "Não Linearizada")
xlabel("tempo [s]")
ylabel("theta_1 [rad]")

%Erro diferentes ODE's linearizada
figure(3)
plot(t1, abs(y1(:,1) - y2(:,1)), "b") %erro

%Erro diferentes ODE's não linearizada
figure(4)
plot(t3, abs(y4(:,1) - y3(:,1)), "b") %erro

%Linearizado
function dy1 = f1(t,y)
    w_p = 1;
    lambda = 9/5;
    dy1_1 = y(3);
    dy1_2 = y(4);
    dy1_3 = (3*w_p^2*(lambda*(2*lambda+3)*y(2)-2*(2*lambda+1)*y(1)))/(3*lambda+4);
    dy1_4 = (3*w_p^2*((4*lambda*(lambda+2)+3)*y(1)-2*(lambda+1)^3*y(2)))/(lambda*(3*lambda+4));
    dy1 = [dy1_1;dy1_2;dy1_3;dy1_4];
end

%Não Linearizado
function dy2 = f2(t,y)
    w_p = 1;
    lamb = 9/5;
    dy2_1 = y(3);
    dy2_2 = y(4);
    dy2_3 = ((3*w_p^2*((5*lamb+4)*sin(y(1))+lamb*(3*sin(y(1)-2*y(2))-4*lamb*sin(y(2)))))/(-15*lamb+9*lamb*cos(2*(y(1)-y(2)))-8))+((6*lamb^2*y(4)^2*sin(y(1)-y(2)))/(-12*lamb+9*lamb*(cos(y(1)-y(2)))^2-4))+((24*lamb^2*y(3)*y(4)*sin(y(1)-y(2)))/(-15*lamb+9*lamb*cos(2*(y(1)-y(2)))-8))+((6*lamb*y(3)^2*sin(y(1)-y(2))*(lamb+3*cos(y(1)-y(2))))/(-12*lamb+9*lamb*(cos(y(1)-y(2)))^2-4));
    dy2_4 = ((3*(-2*lamb*(2+lamb)*sin(y(1))-6*lamb^2*sin(y(1)-2*y(2))-3*(1+2*lamb)*sin(2*y(1)-y(2))+(1+6*lamb+4*lamb^3)*sin(y(2)))*w_p^2)/(lamb*(-8-15*lamb+9*lamb*cos(2*(y(1)-y(2))))))-((6*(2*(2+6*lamb+lamb^3)+9*lamb^2*cos(y(1)-y(2)))*sin(y(1)-y(2))*y(3)^2)/(lamb*(-8-15*lamb+9*lamb*cos(2*(y(1)-y(2))))))-((6*lamb*(2*lamb+3*cos(y(1)-y(2)))*sin(y(1)-y(2))*y(3)*y(4))/(-4-12*lamb+9*(cos(y(1)-y(2)))^2))-((6*lamb*(2*lamb+3*cos(y(1)-y(2)))*sin(y(1)-y(2))*y(4)^2)/(-8-15*lamb+9*lamb*cos(2*(y(1)-y(2)))));
    dy2 = [dy2_1;dy2_2;dy2_3;dy2_4];
end

