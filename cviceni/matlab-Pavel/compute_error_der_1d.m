function errors = compute_error_der_1d(x,u,p, p_der)
  % pocet bodu
  n = length(x);
  % vysledna chyba
  err = 0;
  err_der = 0;
  
  for j = 2:n
    % vytvor kvadraturu 3.radu na intervalu Ej = [x(j-1), x(j)]
    % xx - kvadraturni body, ww - kvadraturni vahy
    [xx,ww] = gauss_1d(3,x(j-1),x(j));
    % pocet bodu kvadratury
    nq = length(xx);
    
    % chyba na aktualnim elementu
    et = 0;
    et_der = 0;
    for q = 1:nq
      % hodnota linearni aproximace v kvadraturnim bode xx(q)
      uxx = u(j-1) + (u(j)-u(j-1)) * (xx(q) - x(j-1)) / (x(j)-x(j-1));
      % hodnota chyby
      exx = (uxx - p(xx(q)));
      exx = exx*exx; % druha mocnina - L2 norma
      % prispevek do integralu
      et = et + exx*ww(q);
      
      % hodnota derivace linearni aproximace v kvadraturnim bode xx(q)
      uxx_der = (u(j)-u(j-1)) / (x(j)-x(j-1));
      % hodnota chyby
      exx_der = (uxx_der - p_der(xx(q)));
      exx_der = exx_der*exx_der; % druha mocnina - L2 norma
      % prispevek do integralu
      et_der = et_der + exx_der*ww(q);
      
    end
    
    err = err + et;
    err_der = err_der + et_der;
  end
  err = sqrt(err);  % L2 norma
  err_der = sqrt(err_der);  % L2 norma
  
  errors = [err; err_der];
end