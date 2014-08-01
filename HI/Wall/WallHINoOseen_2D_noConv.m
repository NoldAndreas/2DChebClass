function HI = WallHINoOseen_2D_noConv(x,y,optsPhys)

    N = length(x);
    id = IoxI(N);

    sigmaH = optsPhys.sigmaHS;

    x0 = optsPhys.pts.x0;
    y0 = optsPhys.pts.y0;
    
    xShift = x0 - x;
    
    % second coordinate is the one over which we integrate, i.e. x and y
    
    y1 = y0 - y;
    y2 = y0 + y; % mirrored y
    
    %Oseen1 = oseen(xShift,y1,sigmaH);
    Oseen2 = oseen(xShift,y2,sigmaH);
    DeltaT = deltaT(x0,y0,x,y,sigmaH);
    
    %HI = Oseen1 - Oseen2 + DeltaT;
   
    HI = - Oseen2 + DeltaT;
%    HI(~isfinite(HI)) = 0;
    
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
        rr = roxr(x,y);
        rInv = rInverse(x,y);
        Oseen = 3/8*sigmaH*rInv.*(id + rr);
        Oseen(~isfinite(Oseen)) = 0;
    end

    function dT = deltaT(x1,y1,x2,y2,sigmaH)
        dT = zeros([N,2,2]);
        sInv  = ((x1-x2).^2  + (y1+y2).^2).^(-1/2); % distance with reflected point
        sInv3 = sInv.^3;
        sInv5 = sInv.^5;
        dT(:,1,1) = -2*y1.*y2.*(sInv3 - 3*(x1-x2).^2.*sInv5);
        dT(:,2,2) =  2*y1.*y2.*(sInv3 - 3*(y1+y2).^2.*sInv5);
        dT(:,1,2) =  2*(x1-x2).*(y2.*sInv3 - 3*y1.*y2.*(y1+y2).*sInv5);
        dT(:,2,1) =  2*(x1-x2).*(y2.*sInv3 + 3*y1.*y2.*(y1+y2).*sInv5);
        dT =  3/8*sigmaH*dT;
        dT(~isfinite(dT)) = 0;
    end


end 