% Vypocet pruhybu nosniku (1D interval (0,1))
% metodou konecnych prvku
%
% -u'' = -1 v (0,1)
%  u(0) = u0, u(1) = u1
%
% u0 = poloha leveho konce nosniku
% u1 = poloha praveho konce nosniku
% n  = pocet delicich uzlu intervalu

u0 = 0;
u1 = 0;
n = 10;

% Vypoctene objekty:
% x = delici body intervalu
% u = vypocteny pruhyb
% A = matice diskretni ulohy
% b = vektor diskretni ulohy



%%%%% Vypocet %%%%%

  % Sestaveni matice
  A = zeros(n+1,n+1);
  L = zeros(2,2);
  for elem=1:n
    % lokalni matice
    L(1,1) = 1*n;
    L(1,2) = -1*n;
    L(2,1) = -1*n;
    L(2,2) = 1*n;
    % pridani do globalni matice
    index = [elem, elem+1];
    for i=1:2
      for j=1:2
        A(index(i),index(j)) = A(index(i),index(j)) + L(i,j);
      end
    end
  end

  
  % Sestaveni prave strany
  b = zeros(n+1,1);
  l = zeros(2,1);
  for elem=1:n
    % lokalni vektor
    l(1,1) = -0.5/n;
    l(2,1) = -0.5/n;
    % pridani do globalniho vektoru
    index = [elem, elem+1];
    for i=1:2
      b(index(i),1) = b(index(i),1) + l(i,1);
    end
  end
  
  % dirichletovy okrajove podminky
  A(1,:) = zeros(1,n+1);
  A(1,1) = 1;
  A(n+1,:) = zeros(1,n+1);
  A(n+1,n+1) = 1;
  b(1,1) = u0;
  b(n+1,1) = u1;
  
  % Vyreseni algebraicke soustavy
  u = A\b;
  
  % Vykresleni reseni
  x = linspace(0,1,n+1);
  key=sprintf('n=%d',n);
  hold all;
  plot(x,u,'-', 'DisplayName', key);
  legend('Location','southeast');

