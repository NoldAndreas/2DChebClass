function [d,e] = GetPolesPadeApproximation(fx,npoles)
        
    N = length(fx-1);
    a = chebpoly(fx);
        
    m = find(abs(a)/max(abs(a))>1e-6,1,'last');
    %if(m< length(a))
    %    disp(['m = ',num2str(m)]);
    %end
    m = min(round(m/2),N-2*npoles);
    
    
    poles = chebpade(m,npoles,a);    
    d     = -real(poles(1:2:npoles)); % REAL PARTS OF POLES
    
    %The minus accounts for the fact that fx is in ordered w.r.t chebychev
    %points x = -1...1, 
    %while chebpade assumes that x = 1...-1
    
    e     = abs(imag(poles(1:2:npoles)));% IMAGINARY PARTS OF POLES
    
    if(e == 0)
        disp('GetPolesPadeApproximation: Imaginary part of pole = 0.');
    end
    
end