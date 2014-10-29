classdef HalfInfCapillary < SpectralSpectral% & ConvolutionPointwise

    properties        
        L1       
        y2Min,y2Max                
    end
    
    methods        
        function this = HalfInfCapillary(Geometry)

            this@SpectralSpectral(Geometry.N(1),Geometry.N(2));

            this.L1      = Geometry.L1; 
            this.y2Min   = Geometry.y2Min; 
            this.y2Max   = Geometry.y2Max;             

            InitializationPts(this);                                    
            this.polar = 'cart';
        end
    end
    
    methods (Access = public)
         function [y1,dy1,dx,ddx,dddx,ddddx] = PhysSpace1(this,x1)
             [y1,dy1,dx,ddx,dddx,ddddx] = QuotientMap(x1,this.L1,0);                          
         end
         function xf = CompSpace1(this,y1)
             xf  = InvQuotientMap(y1,this.L1,0);%  PosRay(z,L1);             
         end
        function [th,dth,dx,ddx,dddx,ddddx] = PhysSpace2(this,x2)                
            [th,dth,dx,ddx,dddx,ddddx] = LinearMap(x2,this.y2Min,this.y2Max);
        end
        function x2 = CompSpace2(this,y2)                        
            x2  = InvLinearMap(y2,this.y2Min,this.y2Max);
        end              
        
        function M_conv = ComputeConvolutionMatrix(this,f,saveBool,ptsCheck)
            if(nargin>=4)
                M_conv  = ComputeConvolutionMatrix_Pointwise(this,f,ptsCheck);
            else
                M_conv  = ComputeConvolutionMatrix_Pointwise(this,f);
            end
            
            if((nargin == 3) && islogical(saveBool) && saveBool)
                this.Conv = M_conv;
            end            
        end 
                     
    end
end