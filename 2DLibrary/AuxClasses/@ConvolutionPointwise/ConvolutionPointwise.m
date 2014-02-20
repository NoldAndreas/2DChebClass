classdef (Abstract) ConvolutionPointwise < handle
    
    methods (Access = public) 
        M_conv  = ComputeConvolutionMatrix_Pointwise(this,f);  
        M_conv  = ComputeConvolutionMatrix_PointwiseX(this,f)
        TestConvolutionMatrix(this,ptsCheck,f,VisualOutput);        
    end
    
    methods (Abstract = true,Access = public)         
         [y1D,dy1D] = PhysSpace1_Sub(y10,y20,x,xCheb,f);
         [x]        = CompSpace1_Sub(y10,y20,y1D);
         
         [y2D,dy2D,d_Pade,ep_Pade] = PhysSpace2_Sub(y10,y20,x,xCheb,f);
         [x]                       = CompSpace2_Sub(y10,y20,y2D,d_Tee,ep_Tee);
         
         ptsCart = GetCartPts(this,pts_y1,pts_y2);         
         d       = GetDistance(this,pts_y1,pts_y2);

    end  
end