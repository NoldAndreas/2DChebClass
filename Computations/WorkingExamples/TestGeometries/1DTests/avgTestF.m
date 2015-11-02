function avgTestF()

    mainGeom.N = 30; mainGeom.yMin = -5; mainGeom.yMax = 5;

    mLine = SpectralLine(mainGeom);

    y = mLine.Pts.y;

    sigma = 1;
    N = 10;

    F = zeros(mainGeom.N,mainGeom.N);

%    fFull = shiftQuad(mLine.Pts.y,0);
    
    for iY = 1:length(y)
    
        % forwards integration
        geom.yMin = y(iY); geom.yMax = min(y(iY)+sigma,mainGeom.yMax); geom.N = N;

        aLine = SpectralLine(geom);

        Int  = aLine.ComputeIntegrationVector;
        
        x = mLine.CompSpace(aLine.Pts.y);
        
        Interp = mLine.ComputeInterpolationMatrix(x,false);
        
        InterPol = Interp.InterPol;
        
%         fFin = InterPol*fFull;
%         
%         fFinExact = shiftQuad(aLine.Pts.y,0);
%        
%         hold off
%         plot(aLine.Pts.y,fFin,'-or');
%         hold on
%         plot(aLine.Pts.y,fFinExact,'ob');
% 
%         pause

        F(iY,:) = Int*InterPol;
        
    end
    
    Eta.F = F;
    
    shift = 2;
    
    f = shiftQuad(y,shift);
    
    fInt = shiftQuadInt(y,shift);
    
    test = Eta.F*f;
    
    
    [fInt,test]
    
    function f = shiftQuad(y,a)
        f = (y-a).^2;
    end
    
    function f = shiftQuadInt(y,a)
        % f(x) = (x-a)^2; integral between y and y+b:        
        b = min(sigma, mainGeom.yMax - y);
        
        f = a^2*b - a*b.^2 + b.^3/3 - 2*a*b.*y + b.^2.*y + b.*y.^2;
    end
    
end