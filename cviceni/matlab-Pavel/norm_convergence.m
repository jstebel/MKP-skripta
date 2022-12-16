clear all;
close all;
clc;

func_1 = @(x,n) x.^n;
func_2 = @(x,n) x.^2 + x/n;

% func_3 = @(x,n) ifelse( x<1.0/(n^3) , n, x.^(-1.0/3.0));


func_4 = @(x,n) exp(-x).*cos(x/n);


% plot_convergence(func_1,0,1,20);
% plot_convergence(func_2,-1,1,20);
plot_convergence(@func_3,0,1,5);
% plot_convergence(func_4,-1,1,20);

function y = plot_convergence(func,a,b,N)
  x = linspace(a,b,1000);
  figure;
  hold on;
  for n = 1:N
    plot(x,func(x,n),'DisplayName', sprintf('func_{%d}',n));
  end
  hold off;  
end

function y = func_3(x,n)
 if(length(x) == 1)
   if(x <= 1/n)
     y = -2*n;
   else
     y = -2*x^(-1.0);
   end
 else
 y=1
   y = zeros(length(x));
   for i=1:length(x)
     if(x(i) <= 1/n)
       y(i) = -2*n;
     else
       y(i) = -2*x(i)^(-1.0);
     end
   end
 end
end
