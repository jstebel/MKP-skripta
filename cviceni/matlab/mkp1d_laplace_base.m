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

  % jakobian (velikost elementu)
  jac=1/n;
  % bazove funkce na ref. elementu
  nphi = 2;                  % pocet bazovych funkci
  phi_ref    = @(x)[x, 1-x]; % funkce
  phidot_ref = @(x)[1, -1];  % derivace
  ndof = n+1; % pocet stupnu volnosti na cele oblasti
  index_elem = @(e)[e, e+1]; % globalni cisla bazovych funkci na elementu e
  %nphi = 3;
  %phi_ref = @(x)[2*x*x-3*x+1, -4*x*x+4*x, 2*x*x-x];
  %phidot_ref = @(x)[4*x-3, -8*x+4, 4*x-1];
  %ndof = 2*n+1;
  %index_elem = @(e)[2*e-1, 2*e, 2*e+1];
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
          L(i,j) = L(i,j) + jac*qw(k)*(phidot(j)*phidot(i));
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
        l(i,1) = l(i,1) + jac*qw(k)*(-phi(i));
      end
    end
    % pridani do globalniho vektoru
    index = index_elem(elem);
    for i=1:nphi
      b(index(i),1) = b(index(i),1) + l(i,1);
    end
  end
  
  % Dirichletovy okrajove podminky
  A(1,:) = zeros(1,ndof);
  A(1,1) = 1;
  A(ndof,:) = zeros(1,ndof);
  A(ndof,ndof) = 1;
  b(1,1) = u0;
  b(ndof,1) = u1;
  
  % Vyreseni algebraicke soustavy
  u = A\b;
  
  % Vykresleni reseni
  x = linspace(0,1,ndof);
  key=sprintf('n=%d',n);
  hold all;
  plot(x,u,'-', 'DisplayName', key);
  legend('Location','southeast');

