function Test_MollifyAndCompare_InfSpace_Planar()


    
    %************************************************
    %***************  Initialization ****************
    %************************************************
    
	PhysArea    = struct('N',[20,20],'L1',2,'L2',2);
    PlotArea    = struct('y1Min',-4,'y1Max',4,'N1',100,...
                       'y2Min',-4,'y2Max',4,'N2',100);
    shapeParams = struct('R',1,'N',[20,10]);
    
    %************************************************
    %****************  Preprocess  ****************
    %************************************************       
    
    tic                  
    IS                        = InfSpace(PhysArea);
    IS.ComputeAll(PlotArea);         
    
    
    AD = DataStorage('MollifierTest',@ComputeConvMatrix,shapeParams,IS);
    
    z = TestFunction(IS.Pts);
    subplot(1,2,1);
    IS.plot(z,'SC');
    subplot(1,2,2);
    IS.plot(AD*z,'SC');
    
    
    
    function AD = ComputeConvMatrix(shapeParams,is)
        AD = is.ComputeConvolutionMatrix(@MollifierXY,shapeParams);
    end
    
    function z = TestFunction(pts)
        z = exp(-(pts.y1_kv.^2 + pts.y2_kv.^2));
    end

    function z = MollifierXY(x,y)                
        z       = Mollifier(sqrt(x.^2 + y.^2));
        %z(r<1)  = exp(-1./(1-(r(r<1)).^2));
        %z(r>=1) = 0;
    end    	
    
end