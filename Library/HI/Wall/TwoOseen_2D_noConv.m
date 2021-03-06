function HI = TwoOseen_2D_noConv(x,y,optsPhys)

    sigmaH = optsPhys.sigmaHS;

    x0 = optsPhys.pts.x0;
    y0 = optsPhys.pts.y0;
    
    x   = x0 - x;
    y1  = y0 - y;
    y2  = y0 + y; 
    
    N = length(x);
    id = IoxI(N);

    Oseen1 = oseen(x,y1,sigmaH);
    Oseen2 = oseen(x,y2,sigmaH);

    HI = Oseen1 - Oseen2;

%     HI(~isfinite(HI)) = 0;

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

   function Oseen = oseen(x,y,sigmaH)
        rr    = roxr(x,y);
        rInv  = rInverse(x,y);
        Oseen = 3/8*sigmaH*rInv.*(id + rr);
        Oseen(~isfinite(Oseen)) = 0;
    end

end