function HI = RotnePrager12(r,optsPhys)

    alpha   = optsPhys.alpha;
    sigmaH  = optsPhys.sigmaHS;
    
    r = abs(r);
    
    rInv    = r.^(-1);
    
    HI      = 3/4*sigmaH.*rInv - 1/4*alpha*sigmaH^3.*rInv.^3;

    HI(r==inf) = 0;
    
end

