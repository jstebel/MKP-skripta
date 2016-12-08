function [xx, ww] = GaussQuad(n, a, b)

   u = 1 : n-1;
   u = u ./ sqrt(4*u.^2 - 1);

   % Same as A = diag(u, -1) + diag(u, 1), but faster (no addition).
   A = zeros(n, n);
   A( 2 : n+1 : n*(n-1) ) = u;
   A( n+1 : n+1 : n^2-1 ) = u;

   % Find the base points X and weight factors W for the interval [-1,1].
   [v, x] = eig(A);
   [x, k] = sort(diag(x));
   w = 2 * v(1,k)'.^2;

   % Linearly transform from [-1,1] to [a,b].
   xx = (b - a) / 2 * x + (a + b) / 2;   % transform base points X
   ww = (b - a) / 2 * w;                 % adjust weigths