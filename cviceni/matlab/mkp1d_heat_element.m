% Vypocet rovnice vedeni tepla (1D interval (x0,x1))
% metodou konecnych prvku
%
% -K*u'' = sigma*(u-u_f) v (x0,x1)
%
% okrajove podminky:
%
% -K*u'(x1) = -q  (Robin)
%  u(x0)    = u_d (Dirichlet)
%
% Dirichletova podminka aproximovana pomoci:
%  K*u'(x0) = (u-u_d)/epsilon
%
% u_d = teplota v bode x0
% q   = tepelny tok do oblasti v bode x1
% n   = pocet delicich uzlu intervalu
%
% x = delici body intervalu
% u = vypoctena teplota
% A = matice diskretni ulohy
% b = vektor diskretni ulohy
%
% Ukoly:
%   - vypocet lokalni matice numerickou integraci

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

jac=(x1-x0)/n; % jakobian (velikost elementu)

  % Sestaveni matice
  A = zeros(n+1,n+1);  
  L = zeros(2,2);
  for elem=1:n
    % lokalni matice
    L(1,1) = (cond/jac^2 + sigma*(1/3))*jac;
    L(1,2) = (-cond/jac^2 + sigma*(1/6))*jac;
    L(2,1) = (-cond/jac^2 + sigma*(1/6))*jac;
    L(2,2) = (cond/jac^2 + sigma*(1/3))*jac;
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
    l(1,1) = sigma*uf*0.5*jac;
    l(2,1) = sigma*uf*0.5*jac;
    % pridani do globalniho vektoru
    index = [elem, elem+1];
    for i=1:2
      b(index(i),1) = b(index(i),1) + l(i,1);
    end
  end
  
  % Dirichletova okrajova podminka pro u0=0, 
  % aproximovana pomoci Robinovy podminky
  A(1,1) = A(1,1) + 1/epsilon;
  b(1,1) = b(1,1) + u0/epsilon;
  
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

