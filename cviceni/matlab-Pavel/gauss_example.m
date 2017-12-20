clear all;
close all;
clc;

% interval
a = 0
b = 5

% CHANGE functions, CHANGE quadrature order

% define function and exact integral

func = @(x) x.^2 +3*x - 5;
int_func = @(x) x.^3/3 +3/2*x.^2 - 5*x;

%func = @(x) x.^6 + 2*x.^5 + 3*x;
%int_func = @(x) x.^7/7 + 2/6*x.^6 + 3/2*x.^2;

%func = @(x) 5*x + sin(x);
%int_func = @(x) 5/2*x.^2 - cos(x);

% numerical intergration
% number of quadrature points
nq = 2
% order of polynomial that will be integrated exactly
exact_order = 2*nq-1

[xx, ww] = gauss_1d(nq,a,b);
s = 0;
for q = 1:nq
  s = s + func(xx(q)) * ww(q);
end

% show results
format long;
s
int_func(b) - int_func(a)