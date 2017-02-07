 classdef Segment < M1SpectralSpectral
    properties 
        h,R
        sphere
    end
    
    methods
        function this = Segment(Geometry)
            this@M1SpectralSpectral(Geometry.N(1),Geometry.N(2));
            
            this.R = Geometry.R;
            if(isfield(Geometry,'Origin'))              
                this.Origin = Geometry.Origin;
            end
            
            this.h = abs(Geometry.h);
            InitializationPts(this);            
            
            if(isfield(Geometry,'sphere'))
                this.sphere = Geometry.sphere;                
            end
            this.polar = 'cart';
        end        
        %***************************************************************
        %   Mapping functions:
        %***************************************************************             
        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)        
            R      = this.R;            
            h      = this.h;                        
            
            n  = length(x1);

            C1      = sqrt(R^2-h^2);
            y1_kv   = C1*x1;

            g       = sqrt(R^2-y1_kv.^2)-h;
            dgdx1   = -C1^2*x1./(sqrt(R^2-y1_kv.^2));
            ddgddx1 = -C1^2*R^2./((R^2-y1_kv.^2).^(3/2));
            
            y2_kv   = ((1+x2)/2).*g + h;

            if(nargout >= 3)
                J        = zeros(n,2,2);
                J(:,1,1) = C1*ones(n,1);
                J(:,2,1) = ((1+x2)/2).*dgdx1;
                J(:,2,2) = g/2;
            end

            if(nargout >= 4)
                dH1        = zeros(n,2,2);                 
            end

            if(nargout >= 4)
                dH2        = zeros(n,2,2);            
                dH2(:,1,1) = ddgddx1.*(1+x2)/2;
                dH2(:,1,2) = dgdx1/2;
                dH2(:,2,1) = dgdx1/2;            
            end

        end
        function [x1,x2] = CompSpace(this,y1,y2)
            exc = MException('Segment:CompSpace','not yet implemented');
            throw(exc);
        end        
        function [int,area] = ComputeIntegrationVector(this)
            int = ComputeIntegrationVector@M1SpectralSpectral(this);
            %Check Accuracy            
            if(this.sphere)
                y1s  = this.Pts.y1_kv;
                y2s  = this.Pts.y2_kv;
                int  = 2*int.*real(sqrt(this.R^2-y1s.^2-y2s.^2))';  
                this.Int = int;
                
                ht   = this.R - this.h;
                area = pi*ht^2/3*(3*this.R - ht);
            else
                th   = 2*acos(this.h/this.R);
                area = this.R^2/2*(th-sin(th));
            end
            if(nargout < 2)
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
end