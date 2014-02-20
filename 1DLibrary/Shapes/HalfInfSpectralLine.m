classdef HalfInfSpectralLine < Spectral
    
    properties
        L
        yMin=0;
    end
    
    methods
        function this = HalfInfSpectralLine(Geometry)
            this@Spectral(Geometry.N);
            
            this.L = Geometry.L;
            
            if(isfield(Geometry,'yMin'))
                this.yMin = Geometry.yMin;
            else
                this.yMin = 0;
            end

            this.polar = 'cart';
            
            InitializationPts(this);  
        end
        
    end
       
   methods (Access = public)
       
        function [th,dth,dx,ddx,dddx,ddddx] = PhysSpace(this,x)
            [th,dth,dx,ddx,dddx,ddddx] = QuotientMap(x,this.L,this.yMin,inf);
        end
        
        function x = CompSpace(this,y)                        
            x  = InvQuotientMap(y,this.L,this.yMin,inf);
        end          
        
   end
   
end