function TestPadeApproximation

    N      = 20;
    npoles = 2;
    
    
    xC  = ClenCurtFlip(N);
    xx  = (-1:0.01:1)';
     
    Interp = barychebevalMatrix(xC,xx);
        
    fx = f(xC);
    
    
	plot(xC,fx,'o','MarkerFace','g','MarkerSize',10); hold on;
    plot(xx,Interp*fx);    hold on;
    
    [d,ep]       = GetPolesPadeApproximation(fx,npoles);
    
    %Plot Improved Interpolation using pole-information    
    
    [xTee,dxTee] = Tref(xC,d,ep);
    xxTee        = InvTref(xx,d,ep);
    InterpTee    = barychebevalMatrix(xC,xxTee);
    
    
    plot(xTee,f(xTee),'o','MarkerFace','r','MarkerSize',5);    
    plot(xx,InterpTee*f(xTee),'r');
    
    plot(xx,f(xx),'--g');
    
    plot([d,d],[min(f(xTee)),max(f(xTee))],'-.k');
    
    function y = f(x)
        %y  = exp(-(xC*10).^2);        
        %y  = exp(-(x*10).^2);    
        %y  = 1./(0.001+x.^2);
        
        sP  = 0.5 + 0.01*1i;
        y  = 1./((x-sP).*(x-conj(sP)));
    end
    
end