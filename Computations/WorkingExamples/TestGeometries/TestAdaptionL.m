function TestAdaptionL

    N = 20;
    L = 1;
    
    [x,w]        = ClenCurtFlip(N-1);
    [L,errorL]   = FindOptimalL(@SqrtMap,x,@f,10,true);
        
    [y,dy]       = SqrtMap(x,L,inf);
    Int          = w.*dy';
    Int(dy==inf) = 0;
    
    yy = (-5:0.1:5)';
    xx = InvSqrtMap(yy,L,inf);    
    InterPol = barychebevalMatrix(x,xx);
    
    xxInf       = (-1:0.01:1)';
    yyInf       = SqrtMap(xxInf,L,inf);
    InterPolInf = barychebevalMatrix(x,xxInf);
    
    
    [val,checkInt] = f(y);
    
    figure;
    subplot(1,2,1);
    plot(y,val,'o'); hold on;
    plot(yy,InterPol*val,'b','LineWidth',1.5);
    plot(yy,f(yy),'g--');
    xlim([min(yy) max(yy)]);
    
    subplot(1,2,2);
    plot(x,val,'o'); hold on;
    plot(xxInf,InterPolInf*val,'b','LineWidth',1.5);
    plot(xxInf,f(yyInf),'g--');
    
    disp(['Error of Integration: ',num2str(Int*val-checkInt)]);
        
    function [error] = GetErrorOfIntegration(Ls)
        [xs,ws]  = ClenCurtFlip(N-1);
        [ys,dys] = SqrtMap(xs,Ls,inf);
        Ints    = ws.*dys';
        Ints(dys==inf) = 0;
        
        [vals,checkInts] = f(ys);
        error = abs(Ints*vals-checkInts);
    
    end

    function [y,integral] = f(x)
        y        = exp(-x.^2);        
        integral = sqrt(pi);
    end

end