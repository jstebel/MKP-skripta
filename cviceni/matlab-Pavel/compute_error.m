function err = compute_error(x, basis, basis_der, index_elem, u, p)
  % pocet elementu
  n = length(x)-1;
  % pocet bazovych funkci
  nphi = length(basis(0));
  % vysledna chyba
  err = [0;0];
  
  % body a vahy Gaussovy kvadratury na ref. elementu
  nq = 3;
  [qx,qw] = gauss_1d(nq,0,1);
  % funkce, ktera mapuje bod x z intervalu [0,1] na interval [a,b]
  real_point = @(a,b,x) a + (b-a)*x;

  for elem = 1:n
    % chyba na aktualnim elementu
    et = 0;
    % odpovidajici stupne volnosti
    idx = index_elem(elem);
    
    jac = x(elem+1) - x(elem);
    
    for q = 1:nq
      xx = real_point(x(elem), x(elem+1), qx(q));
      phi = basis(qx(q));
      phidot = basis_der(qx(q))/jac;
      
      uxx = [0;0];
      for j=1:nphi
        uxx(1) = uxx(1) + u(idx(j)) * phi(j);
        uxx(2) = uxx(2) + u(idx(j)) * phidot(j);
      end
      
      % hodnota chyby
      exx = (uxx - p(xx));
      exx = exx.*exx; % druha mocnina - L2 norma
      % prispevek do integralu
      et = et + exx*qw(q)*jac;
    end
    
    err = err + et;
  end
  err = sqrt(err);  % L2 norma
end