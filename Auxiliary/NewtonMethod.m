function [rho,errorHistory] = NewtonMethod(rho,f,TolFun,MaxIter,lambda,opts)

    if(nargin < 6)
        opts = {};
    end

     if(nargin == 2)
         TolFun  = 10^(-10);
         MaxIter = 2000;
     end     
     if(nargin < 5)
         lambda = 1;
     end

     [v,J] = f(rho);
     err   = max(abs(v));
     i     = 1;
     errorHistory = [];
          
      
	while(err > TolFun && i < MaxIter && ...
            ((i < 5) || (errorHistory(i-3) < errorHistory(i-4)) || ...
                        (errorHistory(i-2) < errorHistory(i-3)) || ...
                        (errorHistory(i-1)< errorHistory(i-2))))
         no = PrintErrorPos(err,['Newton iteration no ',num2str(i)]);         
         
         rho    = rho - lambda*(J\v);
         [v,J]  = f(rho);
         err    = norm(v);
         errorHistory(i) = err;
         %display(num2str(err));             
         i  = i+1;
         
        % for ih = 1:no
         %   fprintf('\b');
        %end                 
    end
    
    if((err > TolFun) && ~IsOption(opts,'returnLastIteration'))
        rho = [];
    end
    
    display(['Error: ',num2str(err)]);
    display(['No of iterations: ' , num2str(i)]);
end