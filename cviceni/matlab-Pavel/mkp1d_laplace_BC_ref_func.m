function [x,u,e] = mkp1d_laplace_BC_ref_func(L,n)

% Poissonova rovnice (1D interval (0,L))
% metodou konecnych prvku
%
% Robinova (Newtonova) okr. podminka v x=0
% a nenulova (nehomogenni) Neumannova okr. podminka v x=L
% K ... vodivost
% f ... zdrojovy clen
%
% -Ku'' = f v (0,L)
%  -(-Ku'(0)) = sigma(u(0) - uD)               % ma vyznam toku ven ("doleva")
%  -Ku'(L) = qN                                % ma vyznam toku ven ("doprava")
%
% fyzikalni motivace - uloha vedeni tepla:
% u ... teplota
% uD ... dana prilehla teplota na levem konci
% sigma ... prechodovy koeficient (velky -> teploty se blizi, maly ->teploty vice rozdilne)
% qN ... dany tepelny tok na prave hranici
%
% n  = pocet intervalu

K = 0.3;
uD = 0.1;
qN = -1;
sigma = 10;

h = L/n;

% definice zdrojoveho clenu f
% f = -1
%f = - ones(n+1,1);
f = @(x) sin(2*x);
% f = 1 na [0,L/2], -1 na [L/2,L]
%f(1:n/2) = -f(1:n/2);

% definice presneho reseni pro porovnani
C1 = qN + 0.5*cos(2*L);
C2 = K * (C1 - 0.5 - sigma*uD)/sigma;
exact = @(x) [1/K * (0.25*sin(2*x)- C1*x - C2);   1/K * (0.5*cos(2*x) - C1)];
%exact_der = @(x) 1/K * (0.5*cos(2*x) - C1);

% nastaveni konecnych prvku, bazovych fci, num. integrace

%% jakobian zobrazeni na ref. element (velikost elementu)
jac=h;

%% LINEARNI KP
%% bazove funkce na ref. elementu
%nphi = 2;                  % pocet bazovych funkci
%phi_ref    = @(x)[1-x; x]; % bazove funkce
%phidot_ref = @(x)[-1; 1];  % derivace bazove fce
%
%ndof = n+1; % pocet stupnu volnosti na cele oblasti
%index_elem = @(e)[e, e+1]; % globalni cisla bazovych funkci na elementu e

% KVADRATICKE KP
nphi = 3;
phi_ref = @(x)[2*x.*x-3*x+1; -4*x.*x+4*x; 2*x.*x-x];
phidot_ref = @(x)[4*x-3; -8*x+4; 4*x-1];
ndof = 2*n+1;
index_elem = @(e)[2*e-1, 2*e, 2*e+1];

% gauss_1d(nq, a, b) na intervalu [a,b], nq - pocet kv. bodu
nq = 3;
[qx,qw] = gauss_1d(nq,0,1); % body a vahy Gaussovy kvadratury na ref. elementu
% funkce, ktera mapuje bod x z intervalu [0,1] na interval [a,b]
real_point = @(a,b,x) a + (b-a)*x;

% Vypoctene objekty:
% x = delici body intervalu
% u = vypoctena teplota
% A = matice diskretni ulohy
% b = vektor diskretni ulohy

x = linspace(0,L,n+1);
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
          M(i,j) = M(i,j) + K*phidot(i)*phidot(j) *qw(q)*jac;
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
  
  % Robin (Newton) BC
  A(1,1) = A(1,1) + sigma;
  b(1) = b(1) + sigma * uD;
  
  % Neumann BC
  b(ndof) = b(ndof) - qN;
  
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
%  xt = linspace(0,L,1e4);
%  plot(xt,exact(xt)(1,1:length(xt)),'-r','DisplayName', 'exact');
%  plot(xt,exact(xt)(2,1:length(xt)),'-g','DisplayName', 'exact_{der}');
%  legend('Location','southeast');
%  hold off;
end

