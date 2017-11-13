function [x,u,e] = mkp1d_2_func(startPoint,endPoint,n)

% MKP na intervalu [startPoint, endPoint]
%
% Resime obecne rovnice:
%                     -p1(x) u'(x) = q(x)
% q'(x) + p2(x) u'(x) + p3(x) u(x) = f(x)
% 
% tedy
%
% -(p1(x)u'(x))'(x) + p2(x) u'(x) + p3(x) u(x) = f(x)
%
% with weak form:
% 
% [-p1 u']_a^b + \int{p1 u'v' + p2 u'v + p3 u v}_\Omega  = \int{f v}_\Omega
%
% p1 ... vodivost
% p2 ... rychlost advekce
% f ... zdrojovy clen
%
% n  = pocet intervalu
%
% Vypoctene objekty:
% x = delici body intervalu
% u = vypoctena teplota
% A = matice diskretni ulohy
% b = vektor diskretni ulohy

if(startPoint >= endPoint) 
  error('Invalid interval!');
end

if(n <= 2)
  error('Invalid number of elements!');
end

if(startPoint != 0)
  error('Need to change analytic solution.')
end

% diskretizace intervalu
L = endPoint-startPoint;
h = L/n;
x = linspace(startPoint, endPoint, n+1);

% definice funkci v rovnici
K = 3;
d = 2.5;
p1 = @(x) 1.0;
p2 = @(x) 0.0;
p3 = @(x) K;
f = @(x) d*x;

% definice OKP
uB = -2;
qN = 0;


% definice presneho reseni pro porovnani
sqrtK = sqrt(K);
expK = @(x) exp(sqrtK*x);
C2 = (d*sqrtK^(-3) + (uB-d*endPoint/K)*expK(-endPoint)) / (1+expK(-2*endPoint));
C1 = C2 - d*sqrtK^(-3);
exact = @(x) [C1*expK(x) + C2*expK(-x) + d/K*x;      sqrtK*(C1*expK(x) - C2*expK(-x)) + d/K];



% nastaveni konecnych prvku, bazovych fci, num. integrace

%% jakobian zobrazeni na ref. element (velikost elementu)
jac=h;

%% LINEARNI KP
%% bazove funkce na ref. elementu
nphi = 2;                  % pocet bazovych funkci
phi_ref    = @(x)[1-x; x]; % bazove funkce
phidot_ref = @(x)[-1; 1];  % derivace bazove fce

ndof = n+1; % pocet stupnu volnosti na cele oblasti
index_elem = @(e)[e, e+1]; % globalni cisla bazovych funkci na elementu e

% KVADRATICKE KP
%nphi = 3;
%phi_ref = @(x)[2*x.*x-3*x+1; -4*x.*x+4*x; 2*x.*x-x];
%phidot_ref = @(x)[4*x-3; -8*x+4; 4*x-1];
%ndof = 2*n+1;
%index_elem = @(e)[2*e-1, 2*e, 2*e+1];

% gauss_1d(nq, a, b) na intervalu [a,b], nq - pocet kv. bodu
nq = 3;
[qx,qw] = gauss_1d(nq,0,1); % body a vahy Gaussovy kvadratury na ref. elementu
% funkce, ktera mapuje bod x z intervalu [0,1] na interval [a,b]
real_point = @(a,b,x) a + (b-a)*x;


% linearni system
A = zeros(ndof);
b = zeros(ndof,1);


%%%%% Vypocet %%%%%

  % Sestaveni matice
  
  % Cyklus pres elementy
  for elem = 1:n
    % sestaveni lokalni matice, rhs
    M = zeros(nphi,nphi);
    m = zeros(nphi,1);
    for q=1:nq  % pres kvadraturni body
      % qx(q) je bod na ref. elementu, xx je bod na realnem elementu
      xx = real_point(x(elem), x(elem+1), qx(q));
      phi = phi_ref(qx(q));
      phidot = phidot_ref(qx(q))/jac;
      % naplneni lok. matice a vektoru
      for i=1:nphi
        for j=1:nphi
          integrand = p1(xx)*phidot(i)*phidot(j) + p2(xx)*phi(i)*phidot(j) + p3(xx)*phi(i)*phi(j); 
          M(i,j) = M(i,j) + integrand * qw(q)*jac;
        end
        m(i) = m(i) + f(xx)*phi(i) * qw(q)*jac;
      end
    end
 
    % odpovidajici stupne volnosti
    idx = index_elem(elem);
    % pridani do globalni matice
    A(idx,idx) = A(idx,idx) + M;
    % pridani do globalni prave strany
    b(idx) = b(idx) + m;
  end
  
  % Neumann BC
  b(1) = b(1) - qN;
  
  % Dirichlet BC
  b(ndof) = uB;
  b(ndof-1) = b(ndof-1) - A(ndof-1,ndof)*uB; 
  A(ndof-1:ndof,ndof-1:ndof) = [A(ndof-1,ndof-1) 0; 0 1];
  
  
  
  % Zobraz soustavu
  %output_precision(1) %Octave
  %[A b]
  
  % Vyreseni algebraicke soustavy
  u = A\b;
  
%  % compute error
%%  e = compute_error_1d(x,u,exact);
%  e = compute_error_der_1d(x,u,exact,exact_der);
  e = compute_error(x,phi_ref, phidot_ref, index_elem, u, exact);
  
%  % Vykresleni reseni
%  key=sprintf('n=%d',n);
%  figure;
%  hold on;
%  % Vykresleni reseni
%  plot(x,u,'-b', 'DisplayName', key);
%  
%  xt = linspace(startPoint,endPoint,1e4);
%  plot(xt,exact(xt)(1,1:length(xt)),'-r','DisplayName', 'exact');
%  plot(xt,exact(xt)(2,1:length(xt)),'-g','DisplayName', 'exact_{der}');
%  legend('Location','southeast');
%  hold off;
end

