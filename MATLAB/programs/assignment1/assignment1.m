%assignment 1 self written code

clearvars, clc

%constants
T = 250+273;
P = 10;
Tc = 407.5;
Pc = 111.3;
R = 0.08206;
tol = 10^(-6);

%calculation of a and b
a = 27*(R^2)*(Tc^2)/(64*Pc);
b = R*Tc/(8*Pc);

%bisection method
tic

func1 = @(v) (P+a/(v^2))*(v-b)-R*T;

y1 = 1; y2=5;
i1 = 0;

while(abs(y1-y2)>tol)
    i1= i1+1;
    yo = (y1+y2)/2;
    if func1(y1)*func1(yo)<0
    y2 = yo;
    end

    if  func1(y2)*func1(yo)<0
    y1 = yo;
    end
end

t1 = toc;

disp("bisection method")
fprintf('v: %.6f\n', yo)
fprintf('Time (s): %.6f\n', t1)
fprintf('Iterations: %.0f\n', i1)
fprintf('\n')

%fixed point iteration
tic

func2 = @(v) R*T/(P + a/(v^2)) + b;
xo = 1;
diff_ = 1;
i2 = 0;

while(diff_> tol)
    i2 = i2+1;
    x1 = func2(xo);
    diff_ = abs(xo-x1);
    xo = x1;
end

t2 = toc;

disp("fixed point method")
fprintf('v: %.6f\n', xo)
fprintf('Time (s): %.6f\n', t2)
fprintf('Iterations: %.0f\n', i2)
fprintf('\n')


%Newton's Method
%func3 = (P+a/(v^2))(v-b)-R*T;
%d_func3 = (P*(v^3)-a*v+2*a*b)/(v^3);
tic;

func4 = @(v) v - ((P+a/(v^2))*(v-b)-R*T)/((P*(v^3)-a*v+2*a*b)/(v^3));

zo = 1;
deviation = 1;
i3=0;

while(deviation>tol)
    i3 = i3+1;
    z1 = func4(zo);
    deviation = abs(z1-zo);
    zo = z1;
end

t3 = toc;

disp("newton method")
fprintf('v: %.6f\n', zo)
fprintf('Time (s): %.6f\n', t3)
fprintf('Iterations: %.0f\n', i3)
fprintf('\n')

%using in built Matlab function
tic;

func5 = @(v) (P+a/(v^2))*(v-b)-R*T;
k = fzero(func5,4);
t4 = toc;

disp("MATLAB inbuilt fzero function")
fprintf('v: %.6f\n', k)
fprintf('Time (s): %.6f\n', t4)
fprintf('Iterations: %.0f\n', false)
fprintf('\n')