function [ d ] = splineDerivative( x, y )
%SPLINEDERIVATIVE Summary of this function goes here
%   Detailed explanation goes here
    
    
    h = x(2:end) - x(1:end-1) ;
    N = length(h)+1;
    del = ( y(2:end) - y(1:end-1) )./ (h);
    
    for k = 1:length(del)
        if isnan(del(k)) == 1
            del(k) = 0;
        end
    end
    e1 = [h(2:end); h(end) + h(end-1)];
    e2 = [h(2); 2*(h(1:end-1) + h(2:end)); h(end-1)];
    e3 = [h(2) + h(1); h(1:end-1)];

    r = [((h(1)+2*e3(1))*h(2)*del(1)+ h(1)^2*del(2))/e3(1) ; ...
      3*(h(2:N-1).*del(1:N-2)+ h(1:N-2).*del(2:N-1)) ; ...
      (h(N-1)^2*del(N-2)+(2*e1(N-1)+h(N-1))*h(N-2)*del(N-1))/e1(N-1) ];
  
   d = tridisolve(e1,e2,e3,r);
end

