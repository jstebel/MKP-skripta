function [x,u,e] = mkp1d_laplace_BC_func(L,n)

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
  
  % vypocet skalarnich soucinu bazovych funkci
  % (pouzivame po castech linearni funkce, tzv. Courantovu bazi)
 for i=1:n+1
    A(i,i)   = 2;
    if(i>1)
      A(i,i-1) = -1;
    end
    if(i<n+1)
      A(i,i+1) = -1;
    end  
  end
  A = K/h * A;
  
  % Sestaveni prave strany
  
  % vypocet skalarnich soucinu s bazovymi funkcemi, krome okrajovych
  for i=2:n
    b(i) = f(i)*h;
  end
  
  % Robin (Newton) BC
  A(1,1) = K/h + sigma;
  b(1) = f(1)*h/2 + sigma * uD;
  
  % Neumann BC
  A(n+1,n+1) = K/h;
  b(n+1) = f(n+1)*h/2 - qN;
  
  % Zobraz soustavu
  %output_precision(1) %Octave
  %[A b]
  
  % Vyreseni algebraicke soustavy
  u = A\b;
  
  % compute error
  x = linspace(0,L,n+1);
%  e = compute_error_1d(x,u,exact);
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
endfunction

