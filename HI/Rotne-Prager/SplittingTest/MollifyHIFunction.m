function HI = MollifyHIFunction(x,y,optsPhys)

    HIfn = optsPhys.HIfn;
    HI   = HIfn(x,y,optsPhys);
    
    r = sqrt(x.^2+y.^2);
    
    R = optsPhys.sigmaHS;
    
    weight       = 1 - exp(-(r./(R-r)).^2);
    weight(r>R)  = 1;
    
%     weight = 1;
    
    n1 = size(HI,2);
    n2 = size(HI,3);
    for i1 = 1:n1
        for i2 = 1:n2
            HI(:,i1,i2) = HI(:,i1,i2).*weight;
        end
    end
    
    HI(r == inf)  = 0;
    HI(r == 0)    = 0;

end