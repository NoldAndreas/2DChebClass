function rho = NewtonMethod(rho,f,TolFun,MaxIter)     
     if(nargin == 2)
         TolFun  = 10^(-10);
         MaxIter = 2000;
     end     

     [v,J] = f(rho);
     err = max(abs(v));
     i = 0;
     while(err > TolFun && i < MaxIter)            
         rho    = rho -  J\v;
         [v,J]  = f(rho);
         err    = norm(v);
         %display(num2str(err));             
         i  = i+1;
     end
     display(['Error: ',num2str(err)]);
     display(['No of iterations: ' , num2str(i)]);
end