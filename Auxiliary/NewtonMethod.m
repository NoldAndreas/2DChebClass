function rho = NewtonMethod(rho,f,TolFun,MaxIter,lambda)     
     if(nargin == 2)
         TolFun  = 10^(-10);
         MaxIter = 2000;
     end     
     if(nargin < 5)
         lambda = 1;
     end

     [v,J] = f(rho);
     err = max(abs(v));
     i = 0;
      
	while(err > TolFun && i < MaxIter)            
         no = fprintf(['Newton iteration no ',num2str(i),' , error: ',num2str(err)]);        
         
         rho    = rho - lambda*(J\v);
         [v,J]  = f(rho);
         err    = norm(v);
         %display(num2str(err));             
         i  = i+1;
         
         for ih = 1:no
            fprintf('\b');
        end 
     end
     display(['Error: ',num2str(err)]);
     display(['No of iterations: ' , num2str(i)]);
end