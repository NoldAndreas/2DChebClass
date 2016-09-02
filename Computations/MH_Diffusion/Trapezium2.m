 classdef Trapezium2 < M1SpectralSpectral
    properties 
        L1
        L2
        h
        shift
        
    end
    
    methods
        function this = Trapezium2(Geometry)
            this@M1SpectralSpectral(Geometry.N(1),Geometry.N(2));
            
            this.L1= Geometry.L1;
            this.L2 = Geometry.L2;
            this.h = Geometry.h;
            if(isfield(Geometry,'shift'))
                this.shift = Geometry.shift;
            else
                this.shift = [0;0];
            end
            InitializationPts(this);            
            
            this.polar = 'cart';
        end        
        %***************************************************************
        %   Mapping functions:
        %***************************************************************             
        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)        
            L1      = this.L1;            
            L2      = this.L2;
            h      = this.h;
            shift = this.shift;
            
            n  = length(x1);
            
            y2Min = -h/2 + shift(2);
            y2Max =  h/2 + shift(2);
            
            %[y2_kv,dy2dx2] =  LinearMap(x2,-h/2,h/2);
            
            [y2_kv,dy2dx2] =  LinearMap(x2,y2Min,y2Max);
            
%            Ly = L2 - (L2-L1)/h*(h/2-y2_kv);
            Ly = L2 - (L2-L1)/h*(y2Max-y2_kv);
            y1_kv   = Ly.*x1/2 + shift(1);

            dLydx2 = (L2-L1)/2*ones(n,1);
            %ddLyddx2 = 0;
            
            % these may all (still) be wrong
            if(nargout >= 3)
                J        = zeros(n,2,2);
                J(:,1,1) = Ly/2;
                J(:,1,2) = x1/2.*dLydx2;
                J(:,2,2) = dy2dx2;
                %J(:,2,1) = zeros(n,1);
            end

            if(nargout >= 4)
                dH2        = zeros(n,2,2);
            end

           if(nargout >= 4)
               dH1        = zeros(n,2,2);            
               %dH1(:,1,1) = 0;
               %dH1(:,2,2) = 0;
               dH1(:,1,2) = dLydx2/2;
               dH1(:,2,1) = dLydx2/2;            
           end

        end
        function [x1,x2] = CompSpace(this,y1,y2)
            exc = MException('Trapezium:CompSpace','not yet implemented');
            throw(exc);
        end    
        
        function Ind    = ComputeIndices(this)
            Ind      = GetIndicesTrapezium(this);
            this.Ind = Ind;
        end    
        
        
        
        function [int,area] = ComputeIntegrationVector(this)
            int = ComputeIntegrationVector@M1SpectralSpectral(this);
            %Check Accuracy            
%             if(this.sphere)
%                 y1s  = this.Pts.y1_kv;
%                 y2s  = this.Pts.y2_kv;
%                 int  = 2*int.*real(sqrt(this.R^2-y1s.^2-y2s.^2))';  
%                 this.Int = int;
% 
%                 ht   = this.R - this.h;
%                 area = pi*ht^2/3*(3*this.R - ht);
%             else
%                 th   = 2*acos(this.h/this.R);
%                 area = this.R^2/2*(th-sin(th));
%             end
%             if(nargout < 2)
%                 if(area == 0)
%                     disp(['Segment: Area is zero, Absolute error is: ',...
%                                     num2str(area-sum(this.Int))]);
%                 else
%                     disp(['Segment: Error of integration of area (ratio): ',...
%                                         num2str(1-sum(this.Int)/area)]);
%                 end
%             end
        end
            
    end
end