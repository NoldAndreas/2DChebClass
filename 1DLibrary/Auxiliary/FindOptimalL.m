function [L,errorL] = FindOptimalL(Map,x,fAna,Lmax,boolPlot)
    
    xx = (-1:0.01:1)';
    
    
    [L,errorL] = fminbnd(@GetErrorOfIntegration,0,Lmax);
    disp(['Optimal Parameter: ',num2str(L),' , Error:',num2str(errorL)]);
    
    if((nargin == 5) && boolPlot)
        Lplot = 0:0.01:Lmax;
        errorL = zeros(size(Lplot));
        for i = 1:length(Lplot)
           errorL(i) =  GetErrorOfIntegration(Lplot(i));
        end
        
        semilogy(Lplot,errorL);
        xlabel('L');
        ylabel('Error');    
    end
    
    function [error] = GetErrorOfIntegration(Ls)        
        ys          = Map(x,Ls);
        yy          = Map(xx,Ls);
                
        InterPol    = barychebevalMatrix(x,xx);
        
        valsY       = fAna(ys);
        valsYY      = fAna(yy);
        
        error       = max(abs(valsYY-InterPol*valsY));
        %[ys,dys]    = Map(x,Ls);
        %Ints             = intX.*dys';
        %Ints(dys==inf)   = 0;
        
        %[vals,checkInts] = fAna(ys);
        %error            = abs(Ints*vals - checkInts);    
    end
end