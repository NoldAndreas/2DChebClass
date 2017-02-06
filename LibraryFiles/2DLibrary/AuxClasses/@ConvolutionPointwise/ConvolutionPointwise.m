classdef (Abstract) ConvolutionPointwise < handle
    
    methods (Access = public) 
        M_conv  = ComputeConvolutionMatrix_Pointwise(this,f);  
        M_conv  = ComputeConvolutionMatrix_PointwiseX(this,f)
        TestConvolutionMatrix(this,ptsCheck,f,VisualOutput);       
        
        function PlotSubGrid(this,y10,y20,f)
            n1  = this.ConvN(1);  
            n2  = this.ConvN(2);

            [x1,wInt1]  = ClenCurtFlip(n1-1);  
            [x2,wInt2]  = ClenCurtFlip(n2-1); 	
            
            [y1D,dy1D]     = PhysSpace1_Sub(this,y10,y20,x1,x1,f);        
            [y2D,dy2D]     = PhysSpace2_Sub(this,y10,y20,x2,x2,f);
            
            Pts.y1_kv      = kronecker(y1D,ones(size(y2D)));
            Pts.y2_kv      = kronecker(ones(size(y1D)),y2D);
            
            h = scatter(Pts.y1_kv,Pts.y2_kv,'ob');              
            set(h,'MarkerEdgeColor','k','MarkerFaceColor','g');
            hold on;
            
            %Plot grid lines
            xI = (-1:0.01:1)';
            O = ones(size(xI));

            yG1 = PhysSpace1_Sub(this,y10,y20,xI,x1,f);  
            yG2 = PhysSpace2_Sub(this,y10,y20,xI,x2,f);  
            
            %(1) Plot x1-isolines
            for i1=1:n1                
                plot(y1D(i1)*O,yG2); hold on;
            end
            
            %(2) Plot x2-isolines
            for i2=1:n2         
                plot(yG1,y2D(i2)*O); hold on;
            end
             
             plot(yG1,O*PhysSpace2_Sub(this,y10,y20,0,x2,f),'r','linewidth',2); hold on;            
             plot(yG1,O*PhysSpace2_Sub(this,y10,y20,-1/sqrt(2),x2,f),':k','linewidth',1.5); hold on;
             plot(yG1,O*PhysSpace2_Sub(this,y10,y20,1/sqrt(2),x2,f),':k','linewidth',1.5); hold on;
             
             plot(O*PhysSpace1_Sub(this,y10,y20,0,x1,f),yG2,'r','linewidth',2); hold on;
             plot(O*PhysSpace1_Sub(this,y10,y20,1/sqrt(2),x1,f),yG2,':k','linewidth',1.5); hold on;
             plot(O*PhysSpace1_Sub(this,y10,y20,-1/sqrt(2),x1,f),yG2,':k','linewidth',1.5); hold on;
             
        end
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