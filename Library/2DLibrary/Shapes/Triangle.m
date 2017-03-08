classdef Triangle < M1SpectralSpectral
    properties
        Y,A
    end
    
    methods
        function this = Triangle(Geometry,N)
            this@M1SpectralSpectral(N(1),N(2));
            this.Y   = Geometry.Y;
            this.A   = GetTriangleParameters(this.Y);
            this.polar = 'cart';
            InitializationPts(this);
        end
                
        %***************************************************************
        %   Mapping functions:
        %***************************************************************             
        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)        
            n  = length(x1);
            O  = ones(n,1);
            A = this.A;

            h = A*[O';x1';x2';(x1.*x2)'];
            y1_kv = h(1,:)';  y2_kv = h(2,:)';

            if(nargout >= 3)
                J        = zeros(n,2,2);
                J(:,1,1) = A(1,2)*O + A(1,4)*x2;
                J(:,1,2) = A(1,3)*O + A(1,4)*x1;
                J(:,2,1) = A(2,2)*O + A(2,4)*x2;
                J(:,2,2) = A(2,3)*O + A(2,4)*x1;
            end

            if(nargout >= 4)
                dH1        = zeros(n,2,2);     
                dH1(:,1,2) = A(1,4)*O;
                dH1(:,2,1) = A(1,4)*O;
            end

            if(nargout >= 4)
                dH2        = zeros(n,2,2);            
                dH2(:,1,2) = A(1,4)*O;
                dH2(:,2,1) = A(1,4)*O;
            end

        end
        function [x1,x2] = CompSpace(y1,y2)
             exc = MException('Triangle:CompSpace','not yet implemented');
            throw(exc);
        end
        
        function [int,area] = ComputeIntegrationVector(this)
            int  = ComputeIntegrationVector@M1SpectralSpectral(this);
            %Check Accuracy
            y = this.Y;
            area = 1/2*abs((y(1,1)-y(3,1))*(y(2,2)-y(1,2))-...
                           (y(1,1)-y(2,1))*(y(3,2)-y(1,2)));       
            if(nargout < 2)
                if(area == 0)
                    disp(['Triangle: Error of integration of area(=0): ',...
                                        num2str(area-sum(this.Int))]);
                else
                    disp(['Triangle: Error of integration of area(ratio): ',...
                                        num2str(1-sum(this.Int)/area)]);
                end
            end
        end
    end
end

