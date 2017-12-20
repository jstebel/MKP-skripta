clear all;
clc;

phi_ref = @(x)[ 2*x.*x - 3*x + 1;
               -4*x.*x + 4*x;
                2*x.*x - x];
phidot_ref = @(x)[ 4*x-3;
                  -8*x+4;
                   4*x-1];

n = 100;
x = linspace(0,1,n);
phi = phi_ref(x);
figure;
hold on;
plot(x,phi(1,1:n),'-r', 'DisplayName', '\phi_1(x)');
plot(x,phi(2,1:n),'-b', 'DisplayName', '\phi_2(x)');
plot(x,phi(3,1:n),'-g', 'DisplayName', '\phi_3(x)');
title('P_2 bazove funkce na ref. elementu [0,1]');
xlabel('x');
legend('Location','east');
hold off;

% save to eps:
%print -deps -color basis_p2.eps
