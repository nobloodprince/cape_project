%2nd problem of 2nd assignment
clc, clearvars
%given constants

D = 0.154;
P0 = 1601325;
P5 = 101325;
L01 =100;
L12 = 300;
L23 = 300;
L45 = 300;
L13 = 1200;
L24 = 1200;
L34 = 1200;
f = 0.005;
p = 1000;

%defining functions
%kij = @(L) (2*f*p*L)/((pi^2)*(D^5));
delta_P = @(L, qij) ((2*f*p*L)/((pi^2)*(D^5)))*(qij)^2;

%q01 = 0; q12 = 0; q13 = 0; q23 = 0; q24 = 0; q34 = 0; q45 = 0;

F = @(q01, q12, q13, q23, q24, q34, q45, P1, P2, P3, P4) [
    P1  - (P0 + delta_P(L01, q01));
    P2 - (P1 + delta_P(L12, q12));
    P3 - (P2 + delta_P(L23, q23));
    P3 - (P1 + delta_P(L13, q13));
    P4 - (P3 + delta_P(L34, q34));
    P4 - (P2 + delta_P(L24, q24));
    P5 - (P4 + delta_P(L45, q45));

    q13 + q12 - q01;
    q23 + q24 - q12;
    q23 + q13 - q34;
    q24 + q34 - q45];

F_wrapped = @(v) F(v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), v(9), v(10), v(11));
vo = ones(1, 11);
solution = fsolve(F_wrapped, vo);

disp(solution)