clear all;
close all;
clc;

% delka intervalu
L = 5
% pocet deleni intervalu
N = [10 20 40 80 160 320 640 1280];

fprintf('n     chyba       rad chyby\n');
err_old = 0;  %chyba z predesleho vypoctu
err = 0;      %chyba z aktualniho vypoctu

for n = N
  % postupne vypocet pro ruzne jemna deleni intervalu
  
  % VYBER ULOHY
%  [x,u,err] = mkp1d_laplace_dirBC_3_func(L,n);
%  [x,u,err] = mkp1d_laplace_dirBC_4_func(L,n);
  [x,u,err] = mkp1d_laplace_BC_func(L,n);

  if (err_old > 0)
    fprintf('%4d %10.3e %8.2g\n', n, err, log(err_old/err)/log(2));
  else
    fprintf('%4d %10.3e\n', n, err);
  end
  
  err_old = err;
end