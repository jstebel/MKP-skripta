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
% n ... pocet intervalu
% n+1 ... pocet bodu
% n-1 ... pocet neznamych

u0 = 0
u1 = 0.2
n = 10
h = 1/n

% definice presneho reseni pro porovnani
exact_t = @(x,u0,u1) 0.5*x.*x +(u1-u0-0.5)*x + u0;
exact = @(x) exact_t(x,u0,u1);

% Vypoctene objekty:
% x = delici body intervalu
% u = reseni
% A = matice diskretni ulohy
% b = vektor diskretni ulohy

A = zeros(n-1);
b = zeros(n-1,1);

%%%%% Vypocet %%%%%

  % Sestaveni matice
  
  % vypocet skalarnich soucinu bazovych funkci
  % (pouzivame po castech linearni funkce, tzv. Courantovu bazi)
  for i=1:n-1
    A(i,i)   = 2;
    if(i>1)
      A(i,i-1) = -1;
    end
    if(i<n-1)
      A(i,i+1) = -1;
    end  
  end
  A = 1/h * A;
  
  % Sestaveni prave strany
  
  % vypocet skalarnich soucinu s bazovymi funkcemi
  for i=1:n-1
    b(i) = -h;
  end
  
  % okrajove podminky
  b(1) = b(1) - u0 * (-1/h);
  b(n-1) = b(n-1) - u1 * (-1/h);
  
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
  plot(x,[u0; u; u1],'-b', 'DisplayName', key);
  
  xt = linspace(0,1,1e4); % dolni mez, horni mez intervalu, pocet bodu
  plot(xt,exact(xt),'-r','DisplayName', 'exact');
  legend('Location','southeast');
  hold off;

