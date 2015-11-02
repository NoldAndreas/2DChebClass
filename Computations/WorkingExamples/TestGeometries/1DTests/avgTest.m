function avgTest()

    infGeom.N = 100; infGeom.L = 5;

    aInfLine = InfSpectralLine(infGeom);

    y = aInfLine.Pts.y;

    sigma = 1;
    N = 100;

    F = zeros(infGeom.N,infGeom.N);

%    fFull = Gaussian(aInfLine.Pts.y,0);
    
    for iY = 1:length(y)
        
        % forwards integration
        geom.yMin = y(iY); geom.yMax = y(iY)+sigma; geom.N = N;

        aLine = SpectralLine(geom);
        
        Int  = aLine.ComputeIntegrationVector;
        
        x = aInfLine.CompSpace(aLine.Pts.y);
        
        Interp = aInfLine.ComputeInterpolationMatrix(x,false);
        
        InterPol = Interp.InterPol;
        
        InterPol(InterPol==inf)=0;
        
%         fFin = InterPol*fFull;
%         
%         fFinExact = Gaussian(aLine.Pts.y,0);
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
    
    f = Gaussian(y,shift);
    
    fInt = GaussianInt(y,shift);
    
    test = Eta.F*f;

    %[fInt, test]
    
    plot(y,fInt,'r');
    hold on
    plot(y,test,'ob');
    
    function f = Gaussian(y,a)
        f = exp(-(y-a).^2);
    end


    function f = GaussianInt(y,a)
        b = sigma;
        f = sqrt(pi)/2 * ( erf(a-y) - erf(a-b-y) );
    end
    
end