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

  % jakobian (velikost elementu)
  jac=(x1-x0)/n;
  % bazove funkce na ref. elementu
  nphi = 2;                  % pocet bazovych funkci
  phi_ref    = @(x)[x, 1-x]; % funkce
  phidot_ref = @(x)[1, -1];  % derivace
  ndof = n+1; % pocet stupnu volnosti na cele oblasti
  index_elem = @(e)[e, e+1]; % globalni cisla bazovych funkci na elementu e
  [qx,qw] = gauss_1d(2,0,1); % body a vahy Gaussovy kvadratury

% Sestaveni matice
  A = zeros(ndof,ndof);
  for elem=1:n
    % lokalni matice
    L = zeros(nphi,nphi);
    for k=1:length(qx)
      phi = phi_ref(qx(k));
      phidot = phidot_ref(qx(k))/jac;
      for i=1:nphi
        for j=1:nphi
          L(i,j) = L(i,j) + (cond*phidot(j)*phidot(i) + sigma*phi(j)*phi(i))*jac*qw(k);
        end
      end
    end
    % pridani do globalni matice
    index = index_elem(elem);
    for i=1:nphi
      for j=1:nphi
        A(index(i),index(j)) = A(index(i),index(j)) + L(i,j);
      end
    end
  end
      
  % Sestaveni prave strany
  b = zeros(ndof,1);
  for elem=1:n
    % lokalni vektor
    l = zeros(nphi,1);
    for k=1:length(qx)
      phi = phi_ref(qx(k));
      for i=1:nphi
        l(i,1) = l(i,1) + sigma*uf*phi(i)*jac*qw(k);
      end
    end
    % pridani do globalniho vektoru
    index = index_elem(elem);
    for i=1:nphi
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

