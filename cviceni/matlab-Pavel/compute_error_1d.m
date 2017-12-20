function err = compute_error_1d(x,u,p)
  % pocet bodu
  n = length(x);
  % vysledna chyba
  err = 0;
  
  for j = 2:n
    % vytvor kvadraturu 3.radu na intervalu Ej = [x(j-1), x(j)]
    % xx - kvadraturni body, ww - kvadraturni vahy
    [xx,ww] = gauss_1d(3,x(j-1),x(j));
    % pocet bodu kvadratury
    nq = length(xx);
    
    % chyba na aktualnim elementu
    et = 0;
    for q = 1:nq
      % hodnota linearni aproximace v kvadraturnim bode xx(q)
      uxx = u(j-1) + (u(j)-u(j-1)) * (xx(q) - x(j-1)) / (x(j)-x(j-1));
      % hodnota chyby
      exx = (uxx - p(xx(q)));
      exx = exx*exx; % druha mocnina - L2 norma
      % prispevek do integralu
      et = et + exx*ww(q);
    end
    
    err = err + et;
  end
  err = sqrt(err);  % L2 norma
end