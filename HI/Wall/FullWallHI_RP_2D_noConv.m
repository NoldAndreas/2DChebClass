function HI = FullWallHI_RP_2D_noConv(x,y,optsPhys)

% See von Hansen, Hinczewski, Netz, JP. 134. 235102 (2011)
% y takes the role of z

    N = length(x);
    id = IoxI(N);

    sigmaH = optsPhys.sigmaHS;

    x0 = optsPhys.pts.x0;
    y0 = optsPhys.pts.y0;
    
    xShift = x0 - x;
    
    % second coordinate is the one over which we integrate, i.e. x and y
    
    y1 = y0 - y;
    y2 = y0 + y; % mirrored y
    
    RP1 = RotnePrager(xShift,y1,sigmaH);
    RP2 = RotnePrager(xShift,y2,sigmaH);
    DeltaD = deltaD(x0,y0,x,y,sigmaH);
    
    HI = RP1 - RP2 + DeltaD;
    
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

    function RP = RotnePrager(x,y,sigmaH)
        rr = roxr(x,y);
        rInv = rInverse(x,y);
        RP = 3/8*sigmaH*rInv.*(id + rr) + 1/16*sigmaH^3*rInv.^3.*(id-3*rr);
        RP(~isfinite(RP)) = 0;
    end

    function dD = deltaD(x1,y1,x2,y2,sigmaH)
        deltaD_O = zeros([N,2,2]);
        deltaD_RP = zeros([N,2,2]);
        
        Rx = x1-x2;
        Ry = y1+y2;
        RInv  = (Rx.^2  + Ry.^2).^(-1/2); % distance with reflected point
        RInv3 = RInv.^3;
        RInv5 = RInv.^5;
        RInv7 = RInv.^7;
        
        deltaD_O(:,1,1) = -y1.*y2.*(RInv3 - 3*Rx.^2.*RInv5);
        deltaD_O(:,2,2) =  y1.*y2.*(RInv3 - 3*Ry.^2.*RInv5);
        deltaD_O(:,1,2) =  Rx.*y2.*(RInv3 - 3*y1.*Ry.*RInv5);
        deltaD_O(:,2,1) =  Rx.*y2.*(RInv3 + 3*y1.*Ry.*RInv5);
        
        deltaD_RP(:,1,1) =  Ry.^2.*(RInv5 - 5*Rx.^2.*RInv7);
        deltaD_RP(:,2,2) = -Ry.^2.*(3*RInv5 - 5*Ry.^2.*RInv7);
        deltaD_RP(:,1,2) = -Rx.*Ry.*(2*RInv5 - 5*Ry.^2.*RInv7);
        deltaD_RP(:,2,1) = -5*Rx.*Ry.^3.*RInv7;
        
        dD =  3/4*sigmaH*deltaD_O + 3/16*sigmaH^3*deltaD_RP;
        dD(~isfinite(dD)) = 0;
    end


end 