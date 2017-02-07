classdef TrefHalfSpace < TrefSpectralSpectral

    properties                
        y2Min = 0;
        L1,L2
        y10 = 0;
        LConv=[];
        L2Conv=[];
    end
    
    methods        
        function this = TrefHalfSpace(Geometry)
            this@TrefSpectralSpectral(Geometry.N(1),Geometry.N(2));
                        
            this.L1    = Geometry.L1;
            this.L2    = Geometry.L2;
            if(isfield(Geometry,'y10'))
                this.y10 = Geometry.y10;
            end
            if(isfield(Geometry,'y2Min'))
                this.y2Min = Geometry.y2Min;
            end
            if(isfield(Geometry,'Conv'))
                this.LConv = Geometry.Conv.L;
                this.L2Conv = Geometry.Conv.L2;
                this.ConvN  = Geometry.Conv.N;
            end
            
            this.polar = 'cart';            
            InitializationPts(this);            
        end                         
    end
    
    methods (Access = public)
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
        
        %*** Mapping function ***
        function [y1,dy1,ddy1] = PhysSpace1Tref(this,x1)                                    
            [y1,dy1,h1s,h2s,h3s,h4s,ddy1] = SqrtMap(x1,this.L1,inf);
            y1 = y1 + this.y10;
        end
        function xf = CompSpace1Tref(this,y1)            
            xf = InvSqrtMap(y1-this.y10,this.L1,inf);
        end
        function [y2,dy2,ddy2] = PhysSpace2Tref(this,x2)
            [y2,dy2,h1s,h2s,h3s,h4s,ddy2] = QuotientMap(x2,this.L2,this.y2Min,inf);            
        end
        function x2 = CompSpace2Tref(this,y2)
            x2  = InvQuotientMap(y2,this.L2,this.y2Min,inf);
        end                
    end
end