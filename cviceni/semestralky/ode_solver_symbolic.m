clear all;
close all;
clc;

% Symbolic ODE solver
pkg load symbolic;


a=0; b=5;
ua = 0;
qb = 2;

%syms t;
%f = exp(7*t)*(sin(t)+7*cos(t)) / (1+7*7);
%df = diff(f);
%simplify(f)


syms y(x)
%ode = -diff(y,2) + 3*y == x;
%ode = -3*diff(y,2) - 7*y == -7*11;
ode = -3*diff(y,2) - 2*diff(y,1) ==  cos(x);
%ode = -3*diff(y,2) + 2*diff(y,1) ==  exp(2/3*x);

u = dsolve(ode)
%u = dsolve(ode, y(a) == ua, diff(y,1)(b) == qb)

ff=function_handle(rhs(u));
C1 = 1;
C2 = 2;


x1=linspace(a,b,100);
y = ff(C1,C2,x1);
plot(x1,y);