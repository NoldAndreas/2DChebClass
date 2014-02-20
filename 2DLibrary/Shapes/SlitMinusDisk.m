 classdef SlitMinusDisk < M1SpectralSpectral
     
    properties 
        R
        OriginDisk = [0,0]
        y2Min,y2Max        
        L1
        LeftRight = 'Right';
    end
    
    methods
        function this = SlitMinusDisk(Geometry)
            this@M1SpectralSpectral(Geometry.N(1),Geometry.N(2));
            
            this.L1         = Geometry.L1;
            this.R          = Geometry.R;
            this.OriginDisk = Geometry.OriginDisk;
            this.y2Min      = Geometry.y2Min;
            this.y2Max      = Geometry.y2Max;
            this.LeftRight  = Geometry.LeftRight;
            
            InitializationPts(this);                        
            this.polar = 'cart';
        end        
        %***************************************************************
        %   Mapping functions:
        %***************************************************************             
        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)        
            
            [y1,dy1]    = QuotientMap(x1,this.L1,0,inf);
            [y2_kv,dy2] = LinearMap(x2,this.y2Min,this.y2Max);
            
            y1Shift    = sqrt(this.R^2 - (y2_kv - this.OriginDisk(2)).^2);        
            
            if(strcmp(this.LeftRight,'Right'))
                y1_kv  = y1 + this.OriginDisk(1) + y1Shift;
            elseif(strcmp(this.LeftRight,'Left'))
                y1_kv  = -y1 + this.OriginDisk(1) - y1Shift;
                dy1    = - dy1;
            else
                error('SlitMinusDisk:choose "Left" or "Right" for this.LeftRight.');
            end
            
%             C1      = sqrt(R^2-h^2);
%             y1_kv   = C1*x1 + Origin(1);
% 
%             g       = S*(sqrt(R^2-y1_kv.^2)-h);
%             dgdx1   = -S*C1^2*x1./(sqrt(R^2-y1_kv.^2));
%             ddgddx1 = -S*C1^2*R^2./((R^2-y1_kv.^2).^(3/2));
% 
%             y2_kv   = ((1+x2)/2).*g + Origin(2) + S*h;
% 
             n = this.N1*this.N2;
             if(nargout >= 3)
                 J        = zeros(n,2,2);
                 J(:,1,1) = dy1;
                 %J(:,2,1) = ((1+x2)/2).*dgdx1;
                 J(:,2,2) = dy2;
             end
% 
             if(nargout >= 4)
                 dH1        = zeros(n,2,2);                 
             end
% 
             if(nargout >= 4)
                 dH2        = zeros(n,2,2);            
%                 dH2(:,1,1) = ddgddx1.*(1+x2)/2;
%                 dH2(:,1,2) = dgdx1/2;
%                 dH2(:,2,1) = dgdx1/2;            
             end
        end
        function [x1,x2] = CompSpace(this,y1,y2)
            exc = MException('Segment:CompSpace','not yet implemented');
            throw(exc);
        end
        
        function [int,area] = ComputeIntegrationVector(this)
            int = ComputeIntegrationVector@M1SpectralSpectral(this);
            %Check Accuracy            
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