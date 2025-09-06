clc, clearvars
T = 250+273; 
P = 10;
Tc = 407.5;
Pc = 111.3;
R = 0.08206;

a = (27*R^2)*(Tc^2)/(64*Pc);
b = R*Tc/(8*Pc);

v = linspace(1,10);

func = (P + (a./(v.^2))).*(v-b) - R*T;

plot(v,func)

x_answer = (func<0.5) & (func>(-0.005));

func(x_answer)