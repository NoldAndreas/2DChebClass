classdef InfSpectralLineQuad < InfSpectralLine
    
%     properties
%         L
%         y0=0;
%     end
%     
    methods
        function this = InfSpectralLineQuad(Geometry)
            %this@Spectral(Geometry.N);
            this@InfSpectralLine(Geometry);
            
            %this.L = Geometry.L;
            %this.polar = 'cart';
            
            %InitializationPts(this);  
        end
        
    end
       
   methods (Access = public)
        
       function [y,dy,dx,ddx,dddx,ddddx] = PhysSpace(this,x)
           [y,dy,dx,ddx,dddx,ddddx] = QuadMap(x,this.L,inf);            
           y = y + this.y0;
       end        
       function xf = CompSpace(this,y)            
            xf  = InvQuadMap(y-this.y0,this.L,inf);
        end               
        
   end
   
end
