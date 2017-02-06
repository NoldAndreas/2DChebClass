function [x] = BasicContinutation(func,parameterRange,xInitial)

     h = waitbar(0,'Computing Continuation...');

    opts   = optimset('Display','off');
    param  = parameterRange(1);
    x1     = fsolve(@f,[param;xInitial],opts);
    x      = zeros(length(x1),length(parameterRange));
    x(:,1) = x1;       
    
    for i=2:length(parameterRange)        
        waitbar(i/length(parameterRange),h);
        param  = parameterRange(i);
        x(:,i) = fsolve(@f,[param;x(2:end,i-1)],opts);
    end
    
    function y = f(x)        
        y = [x(1) - param;func(x(2:end),x(1))];
    end
    
end