 classdef HalfSpaceMinusHalfDisk < Polar_M1SpectralSpectral
    properties 
        R
        Origin = [0,0]
        y2Wall
        L1                
        Rmax = Inf
    end
    
    methods
        function this = HalfSpaceMinusHalfDisk(Geometry)
            %Parameters:
            %L1,R,Origin,Rmax
            this@Polar_M1SpectralSpectral(Geometry.N(1),Geometry.N(2));
            
            this.L1         = Geometry.L1;
            this.R          = Geometry.R;
            this.Origin = Geometry.Origin;
            if(isfield(Geometry,'Rmax'))
                this.Rmax       = Geometry.Rmax;
            end
            
            InitializationPts(this);                        
            this.polar = 'polar';
        end        
        %***************************************************************
        %   Mapping functions:
        %***************************************************************             
        function ptsCart = GetCartPts(this,pts_y1,pts_y2)
            
            if(nargin == 1)
                ptsCart = GetCartPts@Shape(this);
            else
                ptsCart = GetCartPts@Shape(this,pts_y1,pts_y2);
            end
            ptsCart.y1_kv = ptsCart.y1_kv + this.Origin(1);
            ptsCart.y2_kv = ptsCart.y2_kv + this.Origin(2);
            
        end
        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)        
            
            O = ones(size(x1));            
            
            [y2_kv,dy2] = LinearMap(x2,0,pi);
            [y1_kv,dy1] = QuotientMap(x1,this.L1,this.R,this.Rmax);
            
             n = this.N1*this.N2;
             if(nargout >= 3)
                 J        = zeros(n,2,2);
                 J(:,1,1) = dy1;
                 %J(:,2,1) = ((1+x2)/2).*dgdx1;
                 J(:,2,2) = dy2;
             end
 
             if(nargout >= 4)
                 dH1        = zeros(n,2,2);                 
             end
 
             if(nargout >= 4)
                 dH2        = zeros(n,2,2);            
             end
        end
        function [x1,x2] = CompSpace(this,y1,y2)
            exc = MException('Segment:CompSpace','not yet implemented');
            throw(exc);
        end
        
        function [int] = ComputeIntegrationVector(this)
            int = ComputeIntegrationVector@Polar_M1SpectralSpectral(this);
        end
            
    end
end