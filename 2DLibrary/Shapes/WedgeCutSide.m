 classdef WedgeCutSide < Polar_M1SpectralSpectral
    properties 
        R_in = 0;
        R_out
        h        
        leftRight = 'right'
    end
    
    methods
        function this = WedgeCutSide(Geometry)
            this@Polar_M1SpectralSpectral(Geometry.N(1),Geometry.N(2));
            
            if(Geometry.h > Geometry.R_in)
                exc = MException('Segment:CompSpace','not yet implemented');
                throw(exc);
            end
            
            if(isfield(Geometry,'R'))
                this.R_out    = Geometry.R;
            else
                this.R_out    = Geometry.R_out;
            end
            if(isfield(Geometry,'R_in'))
                    this.R_in     = Geometry.R_in;
            end
            this.h        = Geometry.h;
            this.leftRight = Geometry.leftRight;
            %this.th1      = Geometry.th1;
            %this.th2      = Geometry.th1 + mod(Geometry.th2-Geometry.th1,2*pi);  


            InitializationPts(this);  
            if(isfield(Geometry,'Origin'))
                this.Origin = Geometry.Origin;                    
            end       
            
        end        
        %***************************************************************
        %   Mapping functions:
        %***************************************************************             
        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)        
            
            thO     = acos(this.h/this.R_out);
            thI     = acos(this.h/this.R_in);
            
            if(strcmp(this.leftRight,'right'))
               th1 =  3/2*pi+thI;
               th2 =  3/2*pi+thO;
            else
               th1 =  3/2*pi-thO;
               th2 =  3/2*pi-thI;                
            end
            
            y2_kv   = LinearMap(x2,th1,th2);
            %r_out   = (this.h)./abs(cos(y2_kv-3/2*pi));
            r_out   = (this.h)./abs(sin(y2_kv));
            y1_kv   = LinearMap(x1,this.R_in*ones(size(x1)),r_out);
            
            
            if(nargout >= 3)
                J        = zeros(this.M,2,2);
                J(:,1,1) = (r_out-this.R_in)/2;                
                J(:,2,2) = (th2 - th1)/2;
            end

            if(nargout >= 4)
                dH1        = zeros(this.M,2,2);                 
            end

            if(nargout >= 4)
                dH2        = zeros(this.M,2,2);            
            end

        end
        function [x1,x2] = CompSpace(this,y1,y2)
            exc = MException('Segment:CompSpace','not yet implemented');
            throw(exc);
        end
        function [int,area] = ComputeIntegrationVector(this)
            int = ComputeIntegrationVector@Polar_M1SpectralSpectral(this);
            %Check Accuracy            
            h     = this.h;
            Rin   = this.R_in;
            Rout  = this.R_out;
            thO   = 2*acos(h/Rout);
            thI   = 2*acos(h/Rin);
            area  = 0.5*(h*sqrt(Rout^2-h^2) - (Rin^2*thO/2 - Rin^2/2*(thI-sin(thI))));
            %area  = this.R_out^2*sin(Dth)/2 - this.R_in^2*Dth/2;
            
            if(nargout < 2)
                if(area == 0)
                    disp(['WedgeCut: Area is zero, Absolute error is: ',...
                                        num2str(area-sum(this.Int))]);
                else
                    disp(['WedgeCut: Error of integration of area (ratio): ',...
                                        num2str(1-sum(this.Int)/area)]);
                end
            end
        end           
    end
end