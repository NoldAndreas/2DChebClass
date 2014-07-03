 classdef HalfStripMinusDisk < Polar_M1SpectralSpectral
    properties 
        R
        y2Wall
        L1        
        LeftRight
        Rmax = Inf
    end
    
    methods
        function this = HalfStripMinusDisk(Geometry)
            this@Polar_M1SpectralSpectral(Geometry.N(1),Geometry.N(2));
            
            this.L1         = Geometry.L1;
            this.R          = Geometry.R;
            this.Origin     = Geometry.Origin;
            this.y2Wall     = Geometry.y2Wall;
            this.LeftRight  = Geometry.LeftRight;
            %this.Rmax       = Geometry.Rmax;
            
            InitializationPts(this);                        
            this.polar = 'polar';
        end        
        %***************************************************************
        %   Mapping functions:
        %***************************************************************             
        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)        
            
            O = ones(size(x1));
            
            if(strcmp(this.LeftRight,'Left'))
                thMin = pi;
                thMax = pi + asin((this.Origin(2)-this.y2Wall)/this.R);
            elseif(strcmp(this.LeftRight,'Right'))
                thMin = 2*pi - asin((this.Origin(2)-this.y2Wall)/this.R);
                thMax = 2*pi;
            else
                error('HalfStripMinusDisk:Choose Left or Right.');
            end
            
            [y2_kv,dy2] = LinearMap(x2,thMin,thMax);
            rMax        = max(abs((this.Origin(2)-this.y2Wall)./sin(y2_kv)),this.R);
            rMax        = min(rMax,this.Rmax);
            if(strcmp(this.LeftRight,'Left'))
                rMax(y2_kv == pi) = this.Rmax;
            else
                rMax(y2_kv == 2*pi) = this.Rmax;
            end
            rd                = rMax-this.R;
            %L1_r              = min(O*this.L1,rd/3);
            L1_r              = this.L1*rd./(this.L1*3+rd);
            L1_r(rd==inf)     = this.L1;
            [y1_kv,dy1] = QuotientMap(x1,L1_r,O*this.R,max(rMax,this.R));
            
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