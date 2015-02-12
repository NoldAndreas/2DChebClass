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
%        
%         function [th,dth] = PhysSpace(this,x)%,dx,ddx,dddx,ddddx
%             [th,dth] = ExpMap(x,this.L,this.yMin);%QuotientMap(x,this.L,this.yMin,inf);,dx,ddx,dddx,ddddx
%         end
%         
%         function x = CompSpace(this,y)                        
%             x  = InvExpMap(y,this.L,this.yMin);%InvQuotientMap(y,this.L,this.yMin,inf);
%         end          
        
   end
   
end