function w12 = RP12(optsPhys,domain)

    sigmaHij = optsPhys.sigmaHij;
    
    w12Temp = domain.ComputeConvolutionMatrix(@RP,optsPhys.HIShapeParams,false);
    
     w12 = [w12Temp(:,:,1,1), w12Temp(:,:,1,2) ; ...
            w12Temp(:,:,2,1), w12Temp(:,:,2,2) ];


    function rr = roxr(x,y)
        % x and y are the kron products
        rr = zeros([size(x,1),2,2]);
        rm2 = (x.^2 + y.^2).^(-1);
        rr(:,1,1) = x.*x.*rm2;
        rr(:,1,2) = x.*y.*rm2;
        rr(:,2,1) = y.*x.*rm2;
        rr(:,2,2) = y.*y.*rm2;        
    end

    function id = IoxI(N)
        id = zeros(N,2,2);
        id(:,1,1) = 1; id(:,2,2) = 1;
    end

    function rInv = rInverse(x,y)
        rInv = (x.^2 + y.^2).^(-1/2);
        rInv = rInv(:,ones(4,1));
        rInv = reshape(rInv,[],2,2);
    end

    function w12 = RP(x,y)
        N = length(x);
        id = IoxI(N);
        rr = roxr(x,y);
        rInv = rInverse(x,y);
        
        w12 = 3/8*sigmaHij*rInv.*(id + rr) + 1/16*sigmaHij^3*rInv.^3.*(id-3*rr);
        
        w12(isnan(rInv) | rInv == 0)=0;
    end
        

end

