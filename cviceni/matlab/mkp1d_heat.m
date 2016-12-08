% Vypocet rovnice vedeni tepla (1D interval (x0,x1))
% metodou konecnych prvku
%
% -K*u'' = sigma*(u-u_f) v (x0,x1)
%
% okrajove podminky:
%
% -K*u'(x1) = q   (Robin)
%  u(x0)    = u_d (Dirichlet)
%
% Dirichletova podminka aproximovana pomoci:
%  K*u'(x0) = (u-u_d)/epsilon
%
% u_d = teplota v bode x0
% q   = tepelny tok v bode x1
% n   = pocet delicich uzlu intervalu
%
% x = delici body intervalu
% u = vypoctena teplota
% A = matice diskretni ulohy
% b = vektor diskretni ulohy
%
% Ukoly:
%   - zkusit vypocet pro:
%       - bc_flux=0
%       - zmenit uf
%       - menit sigma (sigma = 1000)
%       - nastavit vyrazne vetsi epsilon (epsilon=100)
%       - menit pocet elementu (opravte vypocet jakobianu)
%
%   - aplikace nenulove Dirichletovy okrajove podminky
%   - prepis na sestaveni A a b cyklem přes elementy

% Vstupni parametry

x0 = 0;  % leva hranice
x1 = 1;  % prava hranice
% q = -bc_flux (tok dovnitr je pro bc_flux>0)
bc_flux=1; 
% K
cond = 2;
% sigma
sigma = 1;
% u_f
uf = 1;
% u_d (Dirichlet)
u0 = 10;
% epsilon - aproximace Dirichletovy podminky
epsilon=1e-3;


%%%%% Vypocet %%%%%
n  = 10; % pocet elementu

jac=(x1-x0)/10; % jakobian (velikost elementu)

  % Sestaveni matice
  A = zeros(n+1,n+1);  
  i=1
  A(i,i)   = (cond/jac^2 + sigma*(1.0/3.0))*jac;
  A(i,i+1) = (-cond/jac^2 + sigma*(1.0/6.0))*jac;
  
  for i=2:n
    A(i,i)   = (2*cond/jac^2 + 2*sigma*(1.0/3.0))*jac;
    A(i,i-1) = (-cond/jac^2 + sigma*(1.0/6.0))*jac;
    A(i,i+1) = (-cond/jac^2 + sigma*(1.0/6.0))*jac;
  end
  i=n+1
  A(i,i)   = (cond/jac^2 + sigma*(1.0/3.0))*jac;
  A(i,i-1) = (-cond/jac^2 + sigma*(1.0/6.0))*jac;

  % Dirichletova okrajova podminka pro u0=0, 
  % aproximovana pomoci Robinovy podminky
  A(1,1) = A(1,1) + 1/epsilon;    
      
  
  
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

