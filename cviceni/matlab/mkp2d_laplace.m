% Vypocet pruhybu membrany (ctverec (0,1)x(0,1))
% metodou konecnych prvku
%
% - laplace u = -1 v Omega = (0,1)x(0,1)
%  u = u0 na hranici Omega
%
% u0 = pruhyb membrany na hranici
% n  = pocet elementu v 1 smeru

u0 = 0;
n = 10;

% Vypoctene objekty:
% x,y = delici body intervalu
% u = vypocteny pruhyb
% A = matice diskretni ulohy
% b = vektor diskretni ulohy



%%%%% Vypocet %%%%%
  % souradnice uzlu site
  x = linspace(0,1,n+1);
  y = linspace(0,1,n+1);;
  [xx,yy] = meshgrid(x,y);
  % pocet elementu
  nelem = n*n;
  % indexy hranicnich uzlu (0=vnitrni, 1=vlevo, 2=vpravo, 3=dole, 4=nahore)
  boundary_index = (xx==0) + 2*(xx==1) + 3*(yy==0 & xx>0 & xx<1) + 4*(yy==1 & xx>0 & xx<1);

  % matice transformace souradnic
  J = [1/n 0; 0 1/n];
  Jinv = inv(J);
  % jakobian (velikost elementu)
  jac=det(J);
  
  % bazove funkce na ref. elementu (bilinearni prvky)
  % pocet bazovych funkci
  nphi = 4;
  % funkce
  phi_ref    = @(x,y)[1-x-y+x*y, x-x*y, y-x*y, x*y];
  % derivace
  phidot_ref = @(x,y)[[-1+y,-1+x]; [1-y,-x]; [-y,1-x]; [y,x]];
  % pocet stupnu volnosti na cele oblasti
  ndof = (n+1)^2;
  % globalni cisla bazovych funkci na elementu e
  index_elem = @(e)[(n+1)*floor((e-1)/n)+mod((e-1),n)+1,...
                    (n+1)*floor((e-1)/n)+mod((e-1),n)+2,...
                    (n+1)*(1+floor((e-1)/n))+mod((e-1),n)+1,...
                    (n+1)*(1+floor((e-1)/n))+mod((e-1),n)+2];
  
  % body a vahy Gaussovy kvadratury
  [qx,qw] = gauss_1d(2,0,1);

  % Sestaveni matice
  A = zeros(ndof,ndof);
  for elem=1:nelem
    % lokalni matice
    L = zeros(nphi,nphi);
    for kx=1:length(qx)
      for ky=1:length(qx)
        phi = phi_ref(qx(kx), qx(ky));
        phidot = phidot_ref(qx(kx), qx(ky))*Jinv';
        for i=1:nphi
          for j=1:nphi
            L(i,j) = L(i,j) + dot(phidot(j,:),phidot(i,:))*jac*qw(kx)*qw(ky);
          end
        end
      end
    end
    % pridani do globalni matice
    index = index_elem(elem);
    for i=1:nphi
      if (boundary_index(index(i)) > 0)
        % stupen volnosti lezi na hranici - Dirichletova podminka
        A(index(i), index(i)) = 1;
      else
        for j=1:nphi
          A(index(i),index(j)) = A(index(i),index(j)) + L(i,j);
        end
      end
    end
  end

  
  % Sestaveni prave strany
  b = zeros(ndof,1);
  for elem=1:nelem
    % lokalni vektor
    l = zeros(nphi,1);
    for kx=1:length(qx)
      for ky=1:length(qx)
        phi = phi_ref(qx(kx), qx(ky));
        for i=1:nphi
          l(i,1) = l(i,1) + (-phi(i))*jac*qw(kx)*qw(ky);
        end
      end
    end
    % pridani do globalniho vektoru
    index = index_elem(elem);
    for i=1:nphi
      if (boundary_index(index(i)) > 0)
        % stupen volnosti lezi na hranici - Dirichletova podminka
        b(index(i),1) = u0;
      else
        b(index(i),1) = b(index(i),1) + l(i,1);
      end
    end
  end
  
  % Vyreseni algebraicke soustavy
  u = sparse(A)\b;
  
  % Vykresleni reseni
  surf(xx,yy,reshape(u,size(xx)));

