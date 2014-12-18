function ConvolutionTest

    geom.N = 50;
    geom.L = 4;  % optimize by looking at a derivative, or do parameter testing
    
    infLine = InfSpectralLine(geom);
    
    PlotArea = struct('N',200,'yMin',-10,'yMax',10);
    
    [Pts,Diff,Int,Ind,Interp] = infLine.ComputeAll(PlotArea);
    
    mu1 = 0;
    mu2 = 1;
    
    sigma1 = 0.25;
    sigma2 = 6;
    
    muConv = mu1+mu2;
    sigmaConv = sqrt(sigma1^2 + sigma2^2);
    
    y = Pts.y;
    
    yTest = y(geom.N-20);
    shiftedKernel = Gaussian1(y-yTest);
    figure
    plot(Interp.pts,Interp.InterPol*shiftedKernel,'r');
    hold on
    plot(y,shiftedKernel,'or');

    yTest = y(geom.N-15);
    shiftedKernel = Gaussian1(y-yTest);
    plot(Interp.pts,Interp.InterPol*shiftedKernel,'b');
    hold on
    plot(y,shiftedKernel,'ob');

    
    return
    
    g1 = Gaussian1(y);
    g2 = Gaussian2(y);
    
    gConv = GaussianConv(y);
    
    shapeParams.N  = 200;
    shapeParams.L  = 2;
    
    % false gives pointwise convolution
    convMatrixPtwise = infLine.ComputeConvolutionMatrix(@Gaussian1,shapeParams,false);
    
    % true uses naive convolution
    convMatrixBasic = infLine.ComputeConvolutionMatrix(@Gaussian1,shapeParams,true);

    gPtwise = convMatrixPtwise*g2;
    gBasic = convMatrixBasic*g2;
    
    figure
    plot(y,g1,'r');
    hold on
    plot(y,g2,'b');
    
    figure
    plot(y,gPtwise,'r');
    hold on
    plot(y,gBasic,'--b');
    plot(y,gConv,'og');
    
    ptwiseErr = L2norm(gPtwise-gConv)/L2norm(gConv)
    basicErr  = L2norm(gBasic-gConv)/L2norm(gConv)
    
    
    
    
    
    
    function g = Gaussian1(x)
        g = 1/sqrt(2*pi)/sigma1 * exp( -(x-mu1).^2/(2*sigma1^2) );
    end

    function g = Gaussian2(x)
        g = 1/sqrt(2*pi)/sigma2 * exp( -(x-mu2).^2/(2*sigma2^2) );
    end

    function g = GaussianConv(x)
        g = 1/sqrt(2*pi)/sigmaConv * exp( -(x-muConv).^2/(2*sigmaConv^2) );
    end

    function Norm = L2norm(f)
        Norm = Int*(f.^2);
    end


end