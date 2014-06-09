classdef CornerCircle < SpectralPath
    properties
        R,th1,th2
        Origin = [0,0];
        Corner    = 'SW';
        CornerPos = [0;0];
    end
    
    methods
        function this = CornerCircle(Geometry)
            this@SpectralPath(Geometry.N,'polar');
            
            if(isfield(Geometry,'Origin'))
                this.Origin = Geometry.Origin;
            end
            if(isfield(Geometry,'CornerPos'))
                this.CornerPos = Geometry.CornerPos;
            end
            if(isfield(Geometry,'Corner'))
                this.Corner = Geometry.Corner;
            end
            
            this.R      = Geometry.R;
            
            dw1 = sqrt(this.R^2 - (this.Origin(2)-this.CornerPos(2))^2);
            dw2 = sqrt(this.R^2 - (this.Origin(1)-this.CornerPos(1))^2);
            
            switch(this.Corner)
                case 'SW'
                    yi1 = this.Origin(1) + dw1;
                    yi2 = this.Origin(2) + dw2;
                    this.th1 = - atan((this.Origin(2)-this.CornerPos(2))/(yi1-this.Origin(1)));
                    this.th2 = pi/2 + atan((this.Origin(1)-this.CornerPos(1))/(yi2-this.Origin(2)));

                case 'NW'

                case 'NE'

                case 'SE'

            end
            
            
            InitializationPts(this);            

            this.Pts =  Pol2CartPts(this.Pts);
            this.polar = 'cart';
            this.Pts.y1_kv = this.Pts.y1_kv + this.Origin(1);
            this.Pts.y2_kv = this.Pts.y2_kv + this.Origin(2);
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
    
