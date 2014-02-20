classdef HalfSpaceQuad < SpectralSpectral

    properties                
        y2Min = 0;
        L1,L2
        y10 = 0;
    end
    
    methods        
        function this = HalfSpaceQuad(Geometry,N)
            this@SpectralSpectral(N(1),N(2));
                        
            this.L1    = Geometry.L1;
            this.L2    = Geometry.L2;
            if(isfield(Geometry,'y10'))
                this.y10 = Geometry.y10;
            end
            if(isfield(Geometry,'y2Min'))
                this.y2Min = Geometry.y2Min;
            end
            
            InitializationPts(this);            
        end                         
    end
    
    methods (Access = public)
        function [y1,dy1,dx,ddx,dddx,ddddx] = PhysSpace1(this,x1)            
            [y1,dy1,dx,ddx,dddx,ddddx] = QuadMap(x1,this.L1,inf);            
            y1 = y1 + this.y10;
        end
        function xf = CompSpace1(this,y1)
            xf  = InvQuadMap(y1-this.y10,this.L1,inf);
        end
        function [th,dth,dx,ddx,dddx,ddddx] = PhysSpace2(this,x2)                
            [th,dth,dx,ddx,dddx,ddddx] = QuotientMap(x2,this.L2,...
                                                        this.y2Min,inf);
        end
        function x2 = CompSpace2(this,y2)                        
            x2  = InvQuotientMap(y2,this.L2,this.y2Min,inf);
        end        
        
    end
end