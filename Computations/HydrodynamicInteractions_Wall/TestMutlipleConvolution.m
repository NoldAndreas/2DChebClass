function TestMutlipleConvolution

    PhysArea = struct('N',[20,20],...
                      'L1_Skewed',4,...
                      'L2_Skewed',2,...
                      'L2_AD_Skewed',2.,...
                      'y2wall',0.,...
                      'N2bound',14,'h',1,...
                      'alpha_deg',90);
                  
	PlotArea = struct('y1Min',-5,'y1Max',15,'y2Min',1,'y2Max',15,'N1',80,'N2',80);

    HS                 = HalfSpaceSkewed_FMT(PhysArea,1);
	[Pts,Diff,Int,Ind] = HS.ComputeAll();
	PtsCart            = HS.GetCartPts();
    HS.InterpolationPlotCart(PlotArea,true);
    
    %% Initialize parameters of annulus
	shapeAnn.RMin = 1;
    shapeAnn.L    = 1;
    shapeAnn.N    = [20,20];
    area          = InfAnnulus(shapeAnn);
    
    Conv          = HS.ComputeConvolutionFiniteSupport(area,{'RotnePragner_f1','RotnePragner_f2','RotnePragner_f0'},HS.Pts);
        
    %% Test
    ft   = f(PtsCart);
    
    subplot(2,2,1);    
    HS.doPlots(ft);
    subplot(2,2,2);    
    HS.doPlots(Conv(:,:,2)*ft);
    subplot(2,2,3);    
    HS.doPlots(Conv(:,:,3)*ft);
    subplot(2,2,4);    
    HS.doPlots(Conv(:,:,4)*ft);
    
    
    
    %% Auxiliary function 
    
    function y = f(pts)
        y = exp(- ((pts.y1_kv).^2 + (pts.y2_kv).^2)/20);
    end
    
end