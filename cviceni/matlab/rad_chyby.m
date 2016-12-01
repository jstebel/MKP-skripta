% Vypocet radu chyby porovnanim s referencnim resenim

% vypocet referencniho reseni
% a) numericke reseni na velmi jemne siti
%nref = 10240;
%[xref,uref] = mkp1d(0,0,nref);
%
% b) analyticke reseni
uref = @(x)(-0.5*x*(1-x));

fprintf('n    chyba       rad chyby\n');
difference = 0;
for n = [10 20 40 80 160 320 640 1280]
  % postupne vypocet pro ruzne jemna deleni intervalu
  [x,u] = mkp1d(0,0,n);
  diff_old = difference;
  difference = 0;
  for j = 2:n+1
% varianta a)
%    difference = difference + (u(j)-uref(j*(nref/n)))^2/n;
% varianta b)
    x = (j-1.5)/n;
    difference = difference + ((u(j)+u(j-1))/2-uref(x))^2/n;
  end
  if (diff_old > 0)
    fprintf('%4d %8g %8g\n', n, difference, log(diff_old/difference));
  else
    fprintf('%4d %8g\n', n, difference);
  end
end