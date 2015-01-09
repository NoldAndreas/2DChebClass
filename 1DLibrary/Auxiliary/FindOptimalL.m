function [L,errorL] = FindOptimalL(Map,x,fAna,Lmax,boolPlot)
    
    xx       = (-1:0.01:1)';
    InterPol = barychebevalMatrix(x,xx);
    
    [L,errorL] = fminbnd(@GetErrorOfInterpolation,0,Lmax);
    disp(['Optimal Parameter: ',num2str(L),' , Error:',num2str(errorL)]);
    
    if((nargin == 5) && boolPlot)
        Lplot = 0:0.01:Lmax;
        errorL = zeros(size(Lplot));
        for i = 1:length(Lplot)
           errorL(i) =  GetErrorOfInterpolation(Lplot(i));
        end
        
        figure;
        subplot(2,1,1);
        semilogy(Lplot,errorL);
        xlabel('L');
        ylabel('Error');    
        
        subplot(2,1,2);
        ys          = Map(x,L);
        yy          = Map(xx,L);                               
        valsYY      = fAna(yy);
        plot(ys,fAna(ys),'o'); hold on;
        plot(yy,fAna(yy));
        plot(yy,InterPol*fAna(ys));
        xlim([0 (2*L)]);
    end
    
    function [error] = GetErrorOfInterpolation(Ls)        
        ys          = Map(x,Ls);
        yy          = Map(xx,Ls);               
        
        valsY       = fAna(ys);
        valsYY      = fAna(yy);
        
        error       = max(abs(valsYY-InterPol*valsY));
        %error       = norm(valsYY-InterPol*valsY);
        %[ys,dys]    = Map(x,Ls);
        %Ints             = intX.*dys';
        %Ints(dys==inf)   = 0;
        
        %[vals,checkInts] = fAna(ys);
        %error            = abs(Ints*vals - checkInts);    
    end
end