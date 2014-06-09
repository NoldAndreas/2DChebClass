classdef ArcNSEW < SpectralPath
    properties
        R,th1,th2
        Origin = [0,0]
    end
    
    methods
        function this = ArcNSEW(Geometry)
            this@SpectralPath(Geometry.N,'polar');
            
            this.R      = Geometry.R;
            this.th1    = Geometry.th1;
            this.th2    = Geometry.th2;                       
            
            InitializationPts(this);            
        end
        
        function [y1 ,y2 , dy1_dt , dy2_dt ]  = f_path(this,t)
                        
            y1  = ones(size(t))*this.R;
            y2  = this.th1 + (this.th2-this.th1)*(1+t)/2;
            
            dy1_dt = zeros(size(t));
            dy2_dt = (this.th2-this.th1)/2;
                        
        end
        
        function [int,length] = ComputeIntegrationVector(this)            
            int    = ComputeIntegrationVector@SpectralPath(this);
            length = this.R*(this.th2-this.th1);
            if(length == 0)
                disp(['Arc: Length is zero, Absolute error is: ',...
                                    num2str(length-sum(this.IntSc))]);
            else
                disp(['Arc: Error of integration of length(ratio): ',...
                                        num2str(1-sum(this.IntSc)/length)]);
            end
        end

    end
    
end
    
