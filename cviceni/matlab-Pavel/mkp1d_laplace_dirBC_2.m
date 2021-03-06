% Poissonova rovnice (1D interval (0,1))
% metodou konecnych prvku
%
% Nenulova (nehomogenni) Dirichletova okr. podminka
%
% -u'' = -1 v (0,1)
%  u(0) = u0, u(1) = u1
%
% fyzikalni motivace - uloha vedeni tepla:
% u ... teplota
% u0 ... dana teplota na levem konci
% u1 ... dana teplota na pravem konci
%
% n  = pocet delicich uzlu intervalu

u0 = 0.1
u1 = 0.2
n = 10
h = 1/n

% definice presneho reseni pro porovnani
exact_t = @(x,u0,u1) 0.5*x.*x +(u1-u0-0.5)*x + u0;
exact = @(x) exact_t(x,u0,u1);

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
  A = 1/h * A;
  
  % Sestaveni prave strany
  
  % vypocet skalarnich soucinu s bazovymi funkcemi
  for i=2:n
    b(i,1) = -h;
  end
  
  % dirichletovy okrajove podminky
  A(1,1) = 1;
  b(1) = u0;
    % eliminace u0
    b(2) = b(2) - u0 * A(2,1);
    A(2,1) = 0;
  A(n+1,n+1) = 1;
  b(n+1) = u1;
    % eliminace u1
    b(n) = b(n) - u1 * A(n,n+1);
    A(n,n+1) = 0;
  
  % Zobraz soustavu
  output_precision(1)
  [A b]
  
  % Vyreseni algebraicke soustavy
  u = A\b;
  
  % Vykresleni reseni
  x = linspace(0,1,n+1);
  key=sprintf('n=%d',n);
  figure;
  hold on;
  % Vykresleni reseni
  plot(x,u,'-b', 'DisplayName', key);
  
  xt = linspace(0,1,1e4);
  plot(xt,exact(xt),'-r','DisplayName', 'exact');
  legend('Location','southeast');
  hold off;

