 classdef StripMinusDisk < Polar_M1SpectralSpectral
    properties 
        R        
        y2Wall
        L1        
    end
    
    methods
        function this = StripMinusDisk(Geometry)
            %Parameters:
            %L1,R,Origin,y2Wall
            this@Polar_M1SpectralSpectral(Geometry.N(1),Geometry.N(2));
            
            this.L1         = Geometry.L1;
            this.R          = Geometry.R;
            this.Origin     = Geometry.Origin;
            this.y2Wall     = Geometry.y2Wall;
            
            InitializationPts(this);                        
            this.polar = 'polar';
        end        
        %***************************************************************
        %   Mapping functions:
        %***************************************************************             
        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)        
            
            O = ones(size(x1));
            
            [y2_kv,dy2]   = LinearMap(x2,pi,2*pi);
            
            rd            = abs((this.Origin(2)-this.y2Wall)./sin(y2_kv));
            rd((y2_kv==pi)|(y2_kv==2*pi)) =  inf;
            
            L1_r          = this.L1*(rd-this.R)./(this.L1*3+rd-this.R);
            L1_r(rd==inf) = this.L1;
            [y1_kv,dy1]   = QuotientMap(x1,L1_r,O*this.R,rd);                        
            
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