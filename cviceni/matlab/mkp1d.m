function [x,u,A,b] = mkp1d(u0,u1,n)
% Vypocet pruhybu nosniku (1D interval (0,1))
% metodou konecnych prvku
%
% k  = funkce definujici tuhost nosniku
% f  = pusobici sila (zatizeni)
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

  % Sestaveni matice
  A = zeros(n+1,n+1);
  % dirichletovy okrajove podminky
  A(1,1) = 1;
  A(n+1,n+1) = 1;
  % vypocet skalarnich soucinu bazovych funkci
  % (pouzivame po castech linearni funkce, tzv. Courantovu bazi)
  for i=2:n
    A(i,i)   = integral((i-2)/n,(i-1)/n,@(x)(k(x)*n^2)) + integral((i-1)/n,i/n,@(x)(k(x)*n^2));
    A(i,i-1) = integral((i-2)/n,(i-1)/n,@(x)(-k(x)*n^2));
    A(i,i+1) = integral((i-1)/n,i/n,@(x)(-k(x)*n^2));
  end
  
  % Sestaveni prave strany
  b = zeros(n+1,1);
  % okrajove podminky
  b(1,1) = u0;
  b(n+1,1) = u1;
  % vypocet skalarnich soucinu s bazovymi funkcemi
  for i=2:n
    b(i,1) = integral((i-2)/n,(i-1)/n,@(x)(f(x,n)*(x-(i-2)/n)*n)) + integral((i-1)/n,i/n,@(x)(f(x,n)*(i/n-x)*n));
  end
  
  % Vyreseni algebraicke soustavy
  u = A\b;
  
  % Vykresleni reseni
  x = linspace(0,1,n+1);
  key=sprintf('n=%d',n);
  hold all;
  plot(x,u,'-', 'DisplayName', key);
  legend('Location','southeast');
end


% Pusobici sila (zatizeni)
function r=f(x,n)
  r = -1;
end


% Tuhost nosniku
function r=k(x)
  r=1;
end

% Numericky vypocet integralu (lichobeznikove pravidlo)
function r = integral(a,b,func)
  r = (b-a)*0.5*(func(a)+func(b));
end