% Vypocet pruhybu nosniku (1D interval (0,1))
% metodou konecnych prvku
%
% -u'' = -1 v (0,1)
%  u(0) = u0, u(1) = u1
%
% u0 = poloha leveho konce nosniku
% u1 = poloha praveho konce nosniku
% n  = pocet delicich uzlu intervalu
%
% Navratove hodnoty:
% x = delici body intervalu
% u = vypocteny pruhyb
% A = matice diskretni ulohy
% b = vektor diskretni ulohy

%%%%% Vypocet %%%%%

x0 = 0;  % leva hranice
x1 = 1;  % prava hranice
n  = 10; % pocet elementu
jac=(x1-x0)/10;

% q
bc_flux=10; 
% K
cond = 2;
% sigma
sigma = 0;
% u_f
uf = 1;
% u_0 (Dirichlet)
u0 = 10;
% epsilon - aproximace Dirichletovy podminky
epsilon=1e-3;


  % Sestaveni matice
  A = zeros(n+1,n+1);
 
  % vypocet skalarnich soucinu bazovych funkci
  % (pouzivame po castech linearni funkce, tzv. Courantovu bazi)
  i=1
  A(i,i)   = cond*jac + sigma*jac*(1.0/3.0);    
  A(i,i+1) = -cond*jac + sigma*jac*(1.0/6.0);
  
  for i=2:n    
    A(i,i)   = 2*cond*jac + 2*sigma*jac*(1.0/3.0);    
    A(i,i-1) = -cond*jac + sigma*jac*(1.0/6.0);
    A(i,i+1) = -cond*jac + sigma*jac*(1.0/6.0);
  end
  i=n+1
  A(i,i)   = cond*jac + sigma*jac*(1.0/3.0);    
  A(i,i-1) = -cond*jac + sigma*jac*(1.0/6.0);

  % dirichletovy okrajove podminky
  A(1,1) = A(1,1)+1/epsilon;  
  %A(1,1) = 10; 
  
    
  % Sestaveni prave strany
  b = zeros(n+1,1);
  % vypocet skalarnich soucinu s bazovymi funkcemi  
  i=1
  b(i,1) = sigma*uf*jac/2;
  for i=2:n
    b(i,1) = sigma*uf*jac;
  end
  i=n+1
  b(i,1) = sigma*uf*jac/2;
  % okrajove podminky  
  b(n+1,1) = b(n+1,1) + bc_flux;
  
  % Vyreseni algebraicke soustavy
  u = A\b;
  
  % Vykresleni reseni
  x = linspace(0,1,n+1);
  key=sprintf('n=%d',n);
  %hold all;
  plot(x,u,'-', 'DisplayName', key);
  legend('Location','southeast');

