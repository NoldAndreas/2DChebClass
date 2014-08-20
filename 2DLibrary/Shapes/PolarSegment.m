 classdef PolarSegment < Polar_M1SpectralSpectral
    properties 
        h,R,S        
        sphere
    end
    
    methods
        function this = PolarSegment(Geometry)
            this@Polar_M1SpectralSpectral(Geometry.N(1),Geometry.N(2));
            
            this.R = Geometry.R;
            if(isfield(Geometry,'Origin'))              
                this.Origin = Geometry.Origin;
            end
            if(isfield(Geometry,'S'))
                this.S = Geometry.S;
            elseif(~isfield(Geometry,'S') && Geometry.h<0)               
                this.S = -1;
            elseif(~isfield(Geometry,'S') && Geometry.h >= 0)
                this.S = 1;            
            end
            this.h = abs(Geometry.h);
            InitializationPts(this);            
            
            if(isfield(Geometry,'sphere'))
                this.sphere = Geometry.sphere;                
            end            
        end        
        %***************************************************************
        %   Mapping functions:
        %***************************************************************             
        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)        
            R = this.R;
            S = this.S;
            h = this.h;                        
            
            n  = length(x1);
            
            Phi     = acos(h/R);
            y2_kv   = pi/2 + Phi*x2;
            y2d     = y2_kv - pi/2;
            Rmin    = h./cos(y2d);
            y1_kv   = Rmin + (1+x1).*(R-Rmin)/2;

            if(nargout >= 3)
                J        = zeros(n,2,2);
                J(:,1,1) = (R-Rmin)/2;
                J(:,1,2) = h*sin(y2d)./((cos(y2d)).^2)*Phi.*(1-x1)/2;
                J(:,2,2) = Phi;
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
        
        function [int,area] = ComputeIntegrationVector(this)
            int = ComputeIntegrationVector@Polar_M1SpectralSpectral(this);
            %Check Accuracy            
            if(this.sphere)
                int  = 2*int.*real(sqrt(this.R^2 - (this.Pts.y1_kv).^2))';  
                this.Int = int;
                
                ht   = this.R - this.h;
                area = pi*ht^2/3*(3*this.R - ht);
            else
                th   = 2*acos(this.h/this.R);
                area = this.R^2/2*(th-sin(th));
            end
            if(area == 0)
                disp(['Segment: Area is zero, Absolute error is: ',...
                                    num2str(area-sum(this.Int))]);
            else
                disp(['Segment: Error of integration of area (ratio): ',...
                                    num2str(1-sum(this.Int)/area)]);
            end
        end
            
    end
end