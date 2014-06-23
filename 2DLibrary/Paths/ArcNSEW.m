classdef ArcNSEW < SpectralPath
    properties
        R,h,th1,th2
        % h>0 => inside box
        Origin = [0,0]
        WallPos = 'S';
    end
    
    methods
        function this = ArcNSEW(Geometry)
            this@SpectralPath(Geometry.N,'polar');
           
            if(isfield(Geometry,'Origin'))
                this.Origin = Geometry.Origin;
            end
            if(isfield(Geometry,'WallPos'))
                this.WallPos = Geometry.WallPos;
            end
            
            this.R      = Geometry.R;
            this.h      = Geometry.h;
            
            th = acos(-this.h/this.R);
            
            switch this.WallPos
                case 'N'
                    this.th1   = 3/2*pi-th;            
                    this.th2   = 3/2*pi+th;
                case 'S'
                    this.th1   = pi/2+th;
                    this.th2   = pi/2-th;
                case 'E'
                    this.th1   = pi+th;
                    this.th2   = pi-th;
                case 'W'
                    this.th1   = -th;
                    this.th2   = th;
            end

            InitializationPts(this);   
            
            ShiftArc(this);
            
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
        
        function ShiftArc(this)
            this.Pts = Pol2CartPts(this.Pts);
            this.Pts.y1_kv = this.Pts.y1_kv + this.Origin(1);
            this.Pts.y2_kv = this.Pts.y2_kv + this.Origin(2);
        end

    end
    
end
    
