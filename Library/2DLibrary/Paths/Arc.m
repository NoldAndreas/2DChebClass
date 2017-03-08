classdef Arc < SpectralPath
    properties
        R,th1,th2
        % h>0 => inside box
        Origin = [0,0]        
    end
    
    methods
        function this = Arc(Geometry)
            this@SpectralPath(Geometry.N,'polar');
           
            if(isfield(Geometry,'Origin'))
                this.Origin = Geometry.Origin;
            end
            
            this.R      = Geometry.R;                                    
            
            if(isfield(Geometry,'th1') && isfield(Geometry,'th2'))
                this.th1 = Geometry.th1;
                this.th2 = Geometry.th2;
            elseif(isfield(Geometry,'Corner'))
                CornerPos = Geometry.CornerPos;                    
                Corner    = Geometry.Corner;
                
                dw1 = sqrt(this.R^2 - (this.Origin(2)-CornerPos(2))^2);
                dw2 = sqrt(this.R^2 - (this.Origin(1)-CornerPos(1))^2);

                switch(Corner)                    
                    case 'SW'
                        yi1 = this.Origin(1) + dw1;
                        yi2 = this.Origin(2) + dw2;
                        this.th1 = - atan((this.Origin(2)-CornerPos(2))/(yi1-this.Origin(1)));
                        this.th2 = pi/2 + atan((this.Origin(1)-CornerPos(1))/(yi2-this.Origin(2)));
                    case 'NW'
                        yi1 = this.Origin(1) + dw1;
                        yi2 = this.Origin(2) - dw2;
                        this.th1 = atan( (CornerPos(2) - this.Origin(2))/(yi1 - this.Origin(1)) );
                        this.th2 = - pi/2 - atan( (this.Origin(1) - CornerPos(1))/(this.Origin(2) - yi2) );
                    case 'NE'
                        yi1 = this.Origin(1) - dw1;
                        yi2 = this.Origin(2) - dw2;
                        this.th1 = pi - atan( (CornerPos(2) - this.Origin(2))/(this.Origin(1) - yi1) );
                        this.th2 = 3*pi/2 + atan( (CornerPos(1) - this.Origin(1))/(this.Origin(2) - yi2) );
                    case 'SE'
                        yi1 = this.Origin(1) - dw1;
                        yi2 = this.Origin(2) + dw2;
                        this.th1 = pi/2 - atan( (CornerPos(1) - this.Origin(1))/(yi2 - this.Origin(2)));
                        this.th2 = pi + atan( (this.Origin(2) - CornerPos(2))/(this.Origin(1) - yi1) );
                end                
            elseif(isfield(Geometry,'WallPos'))
                th = acos(-Geometry.h/this.R);
                switch Geometry.WallPos
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
            end
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
            length = this.R*abs(this.th2-this.th1);
            if(length == 0)
                disp(['Arc: Length is zero, Absolute error is: ',...
                                    num2str(length-sum(this.IntSc))]);
            else
                disp(['Arc: Error of integration of length(ratio): ',...
                                        num2str(1-sum(this.IntSc)/length)]);
            end
        end
                
        function ptsCart = GetCartPts(this,pts_y1,pts_y2)
            if(nargin < 3)
                ptsCart  = GetCartPts@SpectralPath(this);
            else
                ptsCart  = GetCartPts@SpectralPath(this,pts_y1,pts_y2);                
            end
                        
            ptsCart.y1_kv = ptsCart.y1_kv + this.Origin(1);
            ptsCart.y2_kv = ptsCart.y2_kv + this.Origin(2);
        end
    end
    
end
    
