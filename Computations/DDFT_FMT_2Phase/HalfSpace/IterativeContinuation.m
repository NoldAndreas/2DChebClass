function x = IterativeContinuation(func,Nosteps,theta,xInitial,params,user_scalarProduct,userXTransform)    

    thetaMin    = theta/32;    
    thetaMax    = theta*32;
    maxDistance = theta*5;
    MaxLoopsWithoutProblems = 10;
    
    LoopsWithoutProblems = 0;
    
    if((nargin < 6) || isempty(user_scalarProduct))
        user_scalarProduct = @times;
    end

    
    if(isfield(params,'Algorithm'))
        opts     = optimset('Display','off','Algorithm',params.Algorithm);   
    else 
        opts     = optimset('Display','off');   
    end
    
    xTTangent = 0*xInitial';    
    xTTangent(1) = 1;
    xTLast    = userXTransform(xInitial);    
    
    dirname = strrep(['IterativeContinuation',func2str(func)], '/','_');    
    x       = DataStorageLoop(dirname,@SolveIterativeStep,params,xInitial,[],Nosteps);            

    function [x,stop] = SolveIterativeStep(params,x_ic)
        stop = false;
        if(params.i == 1)
            x  = fsolve(@fPoint,x_ic,opts);    
            xT = userXTransform(x);
            xTLast = xT;
        else
            repeat = true;
            repeatNo = 0;
            while(repeat)
                x = fsolve(@fStep,x_ic,opts);
                xT       = userXTransform(x);
                distance = usernorm(xT-xTLast);
                if((distance > maxDistance) && (theta/2 >= thetaMin))
                    theta   = theta/2;
                    repeat  = true;
                    repeatNo = repeatNo + 1;
                    LoopsWithoutProblems = 0;
                    disp(['theta decreased to ',num2str(theta)]);
                elseif((distance > maxDistance) && (theta/2 < thetaMin))
                    stop = true;                    
                    LoopsWithoutProblems = 0;
                    break;
                else
                    %either the next step x is not too far from x_ic
                    % or theta is already as small as possible and a jump
                    % is accepted
                    repeat = false;
                end
            end
            if((repeatNo == 0) && (theta*2 <= thetaMax))
                LoopsWithoutProblems = LoopsWithoutProblems +1;                
            end            
            if(LoopsWithoutProblems == MaxLoopsWithoutProblems)
                theta = theta*2;
                disp(['theta increased to ',num2str(theta)]);
                LoopsWithoutProblems = 0;
            end
                    
            xTTangent = (xT-xTLast);
            xTTangent = xTTangent'/usernorm(xTTangent);    
            xTLast    = xT;
        end              
    end
    
    function y = fStep(x)    
        xT = userXTransform(x);
        hs = user_scalarProduct(xTTangent,xT-xTLast) - theta;
        y  = [hs;func(x(2:end),x(1))];
    end

    function y = fPoint(x)
        y = [x(1) - xInitial(1);func(x(2:end),x(1))];
    end

    function n = usernorm(xT)
        n = sqrt(user_scalarProduct(xT,xT));
    end
    
end