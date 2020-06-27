clear
clc

n = 11;
freq = 60;
E = 0.15;
a = 0.28e-2;
b = 0.72e-2;
Z = 48400;
ws = freq*2*pi;
C = 777e-9;

% x = x, x' = y, x'' = dy
dxdt = @(t,x)[
    x(2)
    ws*E*cos(ws.*t) - x(2)./(Z.*C) - (a.*x(1) + b.*x(1).^n)./C;
    ];

opts = odeset('RelTol',1e-9,'AbsTol',1e-12);
[t,s] = ode15s(dxdt,[0 1],[0,1],opts);

plot(t,s(:,2));