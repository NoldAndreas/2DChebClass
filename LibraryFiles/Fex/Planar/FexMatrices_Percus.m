function Eta = FexMatrices_Percus(optsPhys,IDC)
    
    y = IDC.Pts.y;

    sigma = optsPhys.V2.sigmaS;
    sigma = sigma(1); % FIX THIS FOR MULTIPLE SPECIES

    if(isfield(optsPhys,'N'))
        N = opts.N;
    else
        N = 100;
    end

    F = zeros(IDC.N,IDC.N);
    B = zeros(IDC.N,IDC.N);

    if(isprop(IDC,'yMin'))
        yMin = IDC.yMin;
        yMax = IDC.yMax;
    else
        yMin = -inf;
        yMax = inf;
    end
    
    
    for iY = 1:length(y)

        % forwards integration
        %geom.yMin = y(iY); geom.yMax = y(iY)+sigma; geom.N = N;
        geom.yMin = y(iY); geom.yMax = min(y(iY)+sigma,yMax); geom.N = N;

        aLine = SpectralLine(geom);
        Int  = aLine.ComputeIntegrationVector;

        x = IDC.CompSpace(aLine.Pts.y);

        Interp   = IDC.ComputeInterpolationMatrix(x,false);
        InterPol = Interp.InterPol;
        InterPol(InterPol==inf) = 0;

        F(iY,:) = Int*InterPol;

        % backwards integration
        %geom.yMin = y(iY)-sigma; geom.yMax = y(iY); geom.N = N;
        geom.yMin = max(y(iY)-sigma, yMin); geom.yMax = y(iY); geom.N = N;

        aLine = SpectralLine(geom);
        Int  = aLine.ComputeIntegrationVector;

        x = IDC.CompSpace(aLine.Pts.y);

        Interp = IDC.ComputeInterpolationMatrix(x,false);
        InterPol = Interp.InterPol;
        InterPol(InterPol==inf) = 0;

        B(iY,:) = Int*InterPol;

    end

    F(isnan(F)) = 0;
    B(isnan(B)) = 0;

    Eta.F = F;
    Eta.B = B;

end