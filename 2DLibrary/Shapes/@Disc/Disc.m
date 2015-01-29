classdef Disc < Polar_SpectralFourier

    properties        
        R        
        Ch 
        sphere = false;        
    end
    
    methods        
        function this = Disc(Geometry)
            this@Polar_SpectralFourier(Geometry.N(1),Geometry.N(2));
            
            this.R      = Geometry.R;             
            if(isfield(Geometry,'sphere'))
                this.sphere = Geometry.sphere;
            end
            if(isfield(Geometry,'Origin'))
                this.Origin = Geometry.Origin;
            end
            InitializationPts(this);                                    
        end                         
    end
    
    methods (Access = public)
        function [r,dr] = PhysSpace1(this,x)%,dx,ddx,dddx,ddddx]
            %[r,dr,dx,ddx,dddx,ddddx] = LinearMap(x,-this.R,this.R);
            [r,dr] = QuadraticMap(x);
            r      = this.R*r;
            dr     = this.R*dr;
        end
        function xf = CompSpace1(this,r)
            r   = r/this.R;            
            xf = bisection(@QuadraticMap,-1,1,r);%,options)            
            %xf  = fzero(fun,r);
            %xf = InvLinearMap(r,-this.R,this.R);
        end
        function [th,dth,dx,ddx,dddx,ddddx] = PhysSpace2(this,xT)    
            [th,dth,dx,ddx,dddx,ddddx] = LinearMap01(xT,0,2*pi);
        end
        function xf = CompSpace2(this,th)
            xf = th/(2*pi);
        end                        
        function [int,area] = ComputeIntegrationVector(this)            
            %Check Accuracy
            
            if(this.sphere)                
                int = ComputeIntegrationVector@Polar_SpectralFourier(this,true);
                int      = 2*int.*sqrt(this.R^2 - this.Pts.y1_kv.^2)';
                if(sum(imag(int)) > 0)
                    error('Disc:ComputeIntegrationVector:Imaginary part in integration vector');
                end
                this.Int = int;
                area     = 4/3*pi*this.R^3;
            else
                int = ComputeIntegrationVector@Polar_SpectralFourier(this);
                area = pi*this.R^2;                
            end
            if(nargout < 2)
                disp(['Disc: Error of integration of area (ratio): ',...
                                        num2str(1-sum(this.Int)/area)]);                
            end
        end        
        
        acc = testIntegration(this,toCheck);
        c = testCount0(this,N) 
        c = testCount(this,N);      
        c = testCount2(this,N)        
        x = this.Ch(i);
        
    end
end