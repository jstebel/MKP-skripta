clear all;
close all;
clc;

% delka intervalu
%L = 5
% pocet deleni intervalu
n = [10 20 40 80 160 320 640 1280];
% defaultni presnost octave/matlab je float,
% proto nepocitejte, ze chyba klesne radove pod 1e-8
% muze se stat u KP vyssiho radu

fprintf('n     |u-uh|       rad chyby   |u''-uh''|    rad chyby\n');
err = zeros(length(n),1);     % vektor chyb
err_der = zeros(length(n),1); % vektor chyb derivaci
err_H1 = zeros(length(n),1);
derivative = 0; % 0 -off, 1- on

for i = 1:length(n)
  % postupne vypocet pro ruzne jemna deleni intervalu
  
  % VYBER ULOHY
%  [x,u,e] = mkp1d_1_func(0,5,n(i));
%  [x,u,e] = mkp1d_2_func(0,5,n(i));
%  [x,u,e] = mkp1d_3_func(0,10,n(i));
%  [x,u,e] = mkp1d_4_func(1,3,n(i));
%  [x,u,e] = mkp1d_5_func(2,7,n(i));
  [x,u,e] = mkp1d_6_func(-5,5,n(i));
  
  err(i) = e(1);
  if(length(e) > 1)
    derivative = 1;
    err_der(i) = e(2);
    err_H1(i) = norm(e);
  end

  if (i > 1)
    if(derivative)
      fprintf('%4d %10.3e %8.2g %15.3e %8.2g\n', n(i), err(i), log(err(i-1)/err(i))/log(2),err_der(i), log(err_der(i-1)/err_der(i))/log(2));
    else
      fprintf('%4d %10.3e %8.2g\n', n(i), err(i), log(err(i-1)/err(i))/log(2));
    end
  else
    if(derivative)
      fprintf('%4d %10.3e        -      %10.3e        -\n', n(i), err(i), err_der(i));
    else
      fprintf('%4d %10.3e        -\n', n(i), err(i));
    end
  end
end

% plot convergence graph

%figure;
%plot(n,err,'-bo', 'DisplayName', '|u-uh|');
%hold on;
%plot(n,err_der,'-ro', 'DisplayName', '|u''-uh''|');
%legend('Location','northeast');

figure;
loglog(n,err,'-bo', 'DisplayName', '|u-uh|_{L^2}');
hold on;
if(derivative)
  loglog(n,err_der,'-ro', 'DisplayName', '|u''-uh''|_{L^2}');
end
%loglog(n,err_H1,'-gx', 'DisplayName', '|u-uh|_{H^1}');
title('Graf konvergence 1d MKP');
xlabel('Pocet elementu');
ylabel('Chyba v L^2 norme');
legend('Location','northeast');