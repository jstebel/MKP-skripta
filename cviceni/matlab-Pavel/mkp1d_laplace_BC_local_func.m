function [x,u,e] = mkp1d_laplace_BC_local_func(L,n)

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
% n  = pocet delicich uzlu intervalu

K = 0.3;
uD = 0.1;
qN = -1;
sigma = 10;

h = L/n;

% definice zdrojoveho clenu f
% f = -1
f = - ones(n+1,1);
% f = 1 na [0,L/2], -1 na [L/2,L]
%f(1:n/2) = -f(1:n/2);

% definice presneho reseni pro porovnani
C1 = -qN-L;
C2 = K/sigma*(sigma*uD-qN-L);
exact = @(x) 1/K * (0.5*x.*x + C1*x + C2);
exact_der = @(x) 1/K * (x + C1);

% Vypoctene objekty:
% x = delici body intervalu
% u = vypoctena teplota
% A = matice diskretni ulohy
% b = vektor diskretni ulohy

A = zeros(n+1);
b = zeros(n+1,1);

%%%%% Vypocet %%%%%

  % Sestaveni matice
  
  % Lokalni matice - na jednom intervalu 2 stupne volnosti
  nM = 2;
  M = zeros(nM,nM);
  m = zeros(nM,1);
  
  % Cyklus pres elementy
  for elem = 1:n
    % sestaveni lokalni matice
    M = K/h * [1 -1; -1 1];
    % sestaveni lokalni prave strany
    m = f(elem) * h * 0.5 * [1; 1];
    
    % odpovidajici stupne volnosti
    idx = [elem, elem+1];
    % pridani do globalni matice
    A(idx,idx) = A(idx,idx) + M;
    % pridani do globalni prave strany
    b(idx) = b(idx) + m;
  end
  
  % Robin (Newton) BC
  A(1,1) = A(1,1) + sigma;
  b(1) = b(1) + sigma * uD;
  
  % Neumann BC
  b(n+1) = b(n+1) - qN;
  
  % Zobraz soustavu
  %output_precision(1) %Octave
  %[A b]
  
  % Vyreseni algebraicke soustavy
  u = A\b;
  
%  % compute error
  x = linspace(0,L,n+1);
%%  e = compute_error_1d(x,u,exact);
  e = compute_error_der_1d(x,u,exact,exact_der);
  
  % Vykresleni reseni
%  x = linspace(0,L,n+1);
%  key=sprintf('n=%d',n);
%  figure;
%  hold on;
%  % Vykresleni reseni
%  plot(x,u,'-b', 'DisplayName', key);
%  
%  xt = linspace(0,L,1e4);
%  plot(xt,exact(xt),'-r','DisplayName', 'exact');
%  legend('Location','southeast');
%  hold off;
end

