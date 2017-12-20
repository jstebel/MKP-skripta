function [x,u,e] = mkp1d_laplace_dirBC_4_func(L,n)

% Poissonova rovnice (1D interval (0,L))
% metodou konecnych prvku
%
% Nenulova (nehomogenni) Neumannova okr. podminka v x=0
% a nenulova (nehomogenni) Dirichletova okr. podminka v x=L
% K ... vodivost
% f ... zdrojovy clen
%
% -Ku'' = f v (0,L)
%  -(-Ku')(0) = qN                       % ma vyznam toku ven ("doleva")
%   u(L) = u1
%
% fyzikalni motivace - uloha vedeni tepla:
% u ... teplota
% qN ... dany tepelny tok na leve hranici (smerem ven)
% u1 ... dana teplota na pravem konci
%
% n  = pocet delicich uzlu intervalu

K = 0.3;
qN = -0.6;
u1 = 2;
h = L/n;

% definice zdrojoveho clenu f
% f = -1
f = - ones(n+1,1);
% f = 1 na [0,L/2], -1 na [L/2,L]
%f(1:n/2) = -f(1:n/2);

% definice presneho reseni pro porovnani
C1 = qN;
C2 = K*u1 - C1*L - L*L/2;
exact = @(x) 1/K * (0.5*x.*x + C1*x + C2);

% Vypoctene objekty:
% x = delici body intervalu
% u = vypoctena teplota
% A = matice diskretni ulohy
% b = vektor diskretni ulohy

A = zeros(n+1);
b = zeros(n+1,1);

%%%%% Vypocet %%%%%

  % Sestaveni matice
  
  % vypocet skalarnich soucinu bazovych funkci
  % (pouzivame po castech linearni funkce, tzv. Courantovu bazi)
  for i=2:n
    A(i,i)   = 2;
    A(i,i-1) = -1;
    A(i,i+1) = -1;
  end
  A = K/h * A;
  
  % Sestaveni prave strany
  
  % vypocet skalarnich soucinu s bazovymi funkcemi, krome okrajovych
  for i=2:n
    b(i) = f(i)*h;
  end
  
  % Neumann BC
  A(1,1) = K/h;
  A(1,2) =-K/h;
  b(1) = f(1)*h/2 - qN;
  
  % Dirichlet BC
  A(n+1,n+1) = 1;
  b(n+1) = u1;
    % eliminace u1
    b(n) = b(n) - u1 * (-K/h);
    A(n,n+1) = 0;
  
  % Zobraz soustavu
  %output_precision(1); %Octave
  %[A b]
  
  % Vyreseni algebraicke soustavy
  u = A\b;

  % compute error
  x = linspace(0,L,n+1);
  e = compute_error_1d(x,u,exact);
  
  % Vykresleni reseni
%  x = linspace(0,L,n+1);
%  key=sprintf('n=%d',n);
%  figure;
%  hold on;
%  % Vykresleni reseni
%  plot(x,u,'-b', 'DisplayName', key);
%  
%  xt = linspace(0,L,1e4);
%  plot(xt,exact(xt),'-r','DisplayName', 'exact');
%  legend('Location','southeast');
%  hold off;
endfunction
  