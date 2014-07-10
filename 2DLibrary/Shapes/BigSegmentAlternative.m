classdef BigSegmentAlternative < M1SpectralSpectral
    properties
        R,h,sphere
    end
    
    methods
        function this = BigSegmentAlternative(Geometry)
            this@M1SpectralSpectral(Geometry.N(1),Geometry.N(2));
            
            this.Origin = Geometry.Origin;
            this.R      = Geometry.R;
            this.h      = Geometry.h;
            this.polar  = 'cart';
            InitializationPts(this);
            if(isfield(Geometry,'sphere'))
                this.sphere = Geometry.sphere;
            end
            
        end
        
        function [int,area] = ComputeIntegrationVector(this)           
            int = ComputeIntegrationVector@M1SpectralSpectral(this);
            
            if(this.sphere)
                y1s  = this.Pts.y1_kv;
                y2s  = this.Pts.y2_kv;
                int  = 2*int.*real(sqrt(this.R^2-y1s.^2-y2s.^2))';  
                this.Int = int;
                
                ht   = this.R - abs(this.h);
                area = 4/3*pi*this.R^3 - pi*ht^2/3*(3*this.R - ht);
            else
                %area      = area1 + area2;
                hh   = abs(this.h);
                th   = pi - acos(hh/this.R);
                area = this.R^2*(th-0) + ...
                        this.R^2/2*(sin(2*0) - sin(2*th));                
                this.Int  = int;
            end
            if(area == 0)
                disp(['BigSegmentAlternative: Error of integration of area (absolute,area = 0): ',...
                                    num2str(area-sum(this.Int))]);   
            else
                disp(['BigSegmentAlternative: Error of integration of area (%): ',...
                                    num2str(1-sum(this.Int)/area)]);   
            end
        end
        
                
        %***************************************************************
        %   Mapping functions:
        %***************************************************************             
        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)        
            n  = length(x1);

            C1      = (this.R+this.h)/2;

            y2H_kv  = (1+x2)*C1 - this.h;
            y2_kv   = y2H_kv;

            g       = sqrt(this.R^2-y2H_kv.^2);

            dgdx2   = -C1*y2H_kv./(sqrt(this.R^2-y2H_kv.^2));
            dgdx2(isinf(dgdx2)) = 0;
            ddgddx2 = -C1^2*this.R^2./((this.R^2-y2H_kv.^2).^(3/2));

            y1_kv   = x1.*g;

            if(nargout >= 3)
                J        = zeros(n,2,2);
                J(:,1,1) = g;
                J(:,1,2) = x1.*dgdx2;              
                J(:,2,2) = C1*ones(n,1);
            end

            if(nargout >= 4)
                dH1        = zeros(n,2,2);                 
                dH1(:,1,2) = dgdx2;
                dH1(:,2,1) = dgdx2;
                dH1(:,2,2) = ddgddx2;
            end

            if(nargout >= 4)
                dH2        = zeros(n,2,2);
            end

        end
        function [x1,x2] = CompSpace(this,y1,y2)
            exc = MException('BigSegmentAlternative:CompSpace','not yet implemented');
            throw(exc);
        end
    end          
        
end
