classdef InfWedge < Polar_SpectralSpectral

    properties   
        R_in
        LR
        th1,th2                       
    end
    
    methods        
        function this = InfWedge(Geometry)
            this@Polar_SpectralSpectral(Geometry.N(1),Geometry.N(2));
            
            this.R_in     = Geometry.R_in;            
            this.LR       = Geometry.LR;
            this.th1      = Geometry.th1;
            this.th2      = Geometry.th2;  
            
            InitializationPts(this);            
            if(isfield(Geometry,'Origin'))
                this.Origin = Geometry.Origin;                               
            end
        end                         
    end
    
    methods (Access = public)
        function [r,dr,dx,ddx,dddx,ddddx] = PhysSpace1(this,x1)            
            [r,dr,dx,ddx,dddx,ddddx] = QuotientMap(x1,this.LR,this.R_in);             
        end
        function x1 = CompSpace1(this,r)    
            x1  = InvQuotientMap(r,this.LR,this.R_in);            
        end
        function [th,dth,dx,ddx,dddx,ddddx] = PhysSpace2(this,x2)                
            [th,dth,dx,ddx,dddx,ddddx] = LinearMap(x2,this.th1,this.th2);
        end
        function x2 = CompSpace2(this,y2)            
            x2 = InvLinearMap(y2,this.th1,this.th2);
        end        
        
    end
end