% Vypocet rovnice vedení tepla (1D interval (x0,x1))
% metodou konecnych prvku
%
% -Ku'' = sigma*(u_f - u) v (x0,x1)
% 
% okrajove podminky:
%   
%   -Ku'(x1) = -bc_flux    (tok dovnitř je pro bc_flux>0)
%   u(x0) = u0, Dirichlet
%  
%   aproximováno pomocí:
%   Ku'(x0) = (u0 - u)/epsilon
%
%  Úkoly:
%   - zkusit výpočet pro:
%       - bc_flux=0
%       - změnit uf
%       - měnit sigma (sigma = 0.01)
%       - nastavit výrazně větší epsilon (epsilon=100)
%       - měnit počet elementů (opravte výpočet jakobianu)


%%%%% Vypocet %%%%%

x0 = 0;  % leva hranice
x1 = 1;  % prava hranice
n  = 10; % pocet elementu
jac=(x1-x0)/10;

% q
bc_flux=1; 
% K
cond = 2;
% sigma
sigma = 1;
% u_f
uf = 1;
% u_0 (Dirichlet)
u0 = 10;
% epsilon - aproximace Dirichletovy podminky
epsilon=1e-3;


  % Sestaveni matice
  A = zeros(n+1,n+1);  
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

  % Dirichletova okrajová podminka pro u0=0, 
  % aproximovaná pomocí Robinovy podmínky
  A(1,1) = A(1,1)+1/epsilon;    
      
  
  
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
  
  % Neumannova okrajová podmínka
  b(n+1,1) = b(n+1,1) + bc_flux;
  
  
  
  % Vyreseni algebraicke soustavy
  u = A\b;
  
  % Vykresleni reseni
  x = linspace(0,1,n+1);
  key=sprintf('n=%d',n);
  %hold all;
  plot(x,u,'-', 'DisplayName', key);
  legend('Location','southeast');

