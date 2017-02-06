function TestBarkerHenderson1DConvolution()    

    shapeParamsHalfLim.N = 80;
    shapeParamsHalfLim.L = -10;
    shapeParamsHalfLim.yMin = -1;
    
    PlotArea = struct('N',100,'yMin',1,'yMax',10);
    
    HIS = HalfInfSpectralLine(shapeParamsHalfLim);     
    %Ind  = HIS.ComputeIndices();            
	Int  = HIS.ComputeIntegrationVector();
    HIS.InterpolationPlot(PlotArea,true);  
                        
   % HIS.ComputeAll(PlotArea);
    y = HIS.Pts.y;
    x = HIS.Pts.x;
    
    [V,VInt] = BH1(y);
    %V        = Gaussian1(y);
    
    [th,dth] = HIS.PhysSpace(HIS.Pts.x);
    
    %dth = y./(1-x.^2);
    plot(HIS.Pts.x,V.*dth);
    
    
    
    HIS.plot(V.*dth);    
    
    PrintErrorPos(HIS.Int*V-VInt,'Integration of 1/x^4');
    
    function [g,VInt] = BH1(x)
        %x = abs(x);
        %g = x.*BarkerHenderson_2D(x); g(x==inf) = 0; g(x==-inf) = 0;
        %g(x <= 1) = BarkerHenderson_2D(1);
        g    = 1./x.^4;
        VInt = 1/3;
    end

    function g = Gaussian1(x)
        g = exp( -x.^2);
    end
end