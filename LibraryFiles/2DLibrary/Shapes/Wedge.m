classdef Wedge < Polar_SpectralSpectral

    properties   
        R_in = 0;
        R_out
        th1,th2                
    end
    
    methods        
        function this = Wedge(Geometry)
            this@Polar_SpectralSpectral(Geometry.N(1),Geometry.N(2));

            if(isfield(Geometry,'y1Min'))
                this.R_in  = Geometry.y1Min;
                this.R_out = Geometry.y1Max;
                this.th1   = Geometry.y2Min;
                this.th2   = Geometry.y2Max;                
            else                
                if(isfield(Geometry,'R'))
                    this.R_out    = Geometry.R;
                else
                    this.R_out    = Geometry.R_out;
                end
                if(isfield(Geometry,'R_in'))
                        this.R_in     = Geometry.R_in;
                end
                this.th1      = Geometry.th1;
                this.th2      = Geometry.th1 + mod(Geometry.th2-Geometry.th1,2*pi);  
            end

            InitializationPts(this);  
            if(isfield(Geometry,'Origin'))
                this.Origin = Geometry.Origin;                    
            end       
            
        end                         
    end
    
    methods (Access = public)
        function [r,dr,dx,ddx,dddx,ddddx] = PhysSpace1(this,x1)            
            [r,dr,dx,ddx,dddx,ddddx] = LinearMap(x1,this.R_in,this.R_out);
        end
        function x1 = CompSpace1(this,r)    
            x1  = InvLinearMap(r,this.R_in,this.R_out);
        end
        function [th,dth,dx,ddx,dddx,ddddx] = PhysSpace2(this,x2)                
            [th,dth,dx,ddx,dddx,ddddx] = LinearMap(x2,this.th1,this.th2);
        end
        function x2 = CompSpace2(this,y2)            
            x2 = InvLinearMap(y2,this.th1,this.th2);
        end                
        function [int,area] = ComputeIntegrationVector(this,t1Odd,t2Odd)
            if(nargin == 1)
                int  = ComputeIntegrationVector@Polar_SpectralSpectral(this);
            elseif(nargin == 2)
                int  = ComputeIntegrationVector@Polar_SpectralSpectral(this,t1Odd);
            elseif(nargin == 3)
                int  = ComputeIntegrationVector@Polar_SpectralSpectral(this,t1Odd,t2Odd);
            end            
            area = (this.th2-this.th1)/2*(this.R_out^2 - this.R_in^2);                
            if(nargout == 1)
                if(area == 0)
                    disp(['Wedge: Area is zero. Error of integration of area(abs): ', num2str(sum(this.Int))]);                                        
                else
                    disp(['Wedge: Error of integration of area(ratio): ', num2str(1-sum(this.Int)/area)]);                                        
                end            
%            elseif(t1Odd)
%                int  = 2*int.*real(sqrt(this.R^2-y1s.^2-y2s.^2))';  
%            elseif(t2Odd)
            end
        end        
        
    end
end