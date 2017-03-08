classdef DoubleCutCircle < handle
   
    properties 
        Origin = [0;0];
        R
        Corner    = 'SW';
        CornerPos = [0;0];
        A1,A2
        Int
        Pts     
    end
    
    methods 
        function this = DoubleCutCircle(Geometry)
            if(isfield(Geometry,'Origin'))
                this.Origin = Geometry.Origin;
            end
            if(isfield(Geometry,'CornerPos'))
                this.CornerPos = Geometry.CornerPos;
            end
            if(isfield(Geometry,'Corner'))
                this.Corner = Geometry.Corner;
            end

            % Compute geometry for segment:
            this.R      = Geometry.R;
                               
            dw1 = sqrt(this.R^2 - (this.Origin(2)-this.CornerPos(2))^2);
            dw2 = sqrt(this.R^2 - (this.Origin(1)-this.CornerPos(1))^2);

            %xm = this.Origin(1) - dw1;
            xp = this.Origin(1) + dw1;
            %ym = this.Origin(2) - dw2;
            yp = this.Origin(2) + dw2;

            switch(this.Corner)
                case 'SW'
                    theta1 = atan( (yp - this.Origin(2))/(this.Origin(1) - this.CornerPos(1)));
                    theta2 = atan( (xp - this.Origin(1))/(this.Origin(2) - this.CornerPos(2)));
                    arcSh1.th1 = pi + theta1;
                    arcSh1.th2 = 3*pi/2 - theta2;
                    arcSh2.th1 = pi - theta1;
                    arcSh2.th2 = -pi/2 + theta2; 
                case 'NW'
                    theta1 = atan( (yp - this.Origin(2))/(this.Origin(1) - this.CornerPos(1)));
                    theta2 = atan( (xp - this.Origin(1))/(this.CornerPos(2) - this.Origin(2)));
                    arcSh1.th1 = pi/2 + theta2;
                    arcSh1.th2 = pi - theta1;
                    arcSh2.th1 = pi/2 - theta2;
                    arcSh2.th2 = -pi + theta1; 
                case 'NE'                    
                    theta1 = atan( (yp - this.Origin(2))/(this.CornerPos(1) - this.Origin(1)));
                    theta2 = atan( (xp - this.Origin(1))/(this.CornerPos(2) - this.Origin(2)));
                    arcSh1.th1 = pi/2 - theta2;
                    arcSh1.th2 = theta1;
                    arcSh2.th1 = 2*pi-theta1;
                    arcSh2.th2 = pi/2 + theta2; 
                case 'SE'
                    theta1 = atan( (yp - this.Origin(2))/(this.CornerPos(1) - this.Origin(1)));
                    theta2 = atan( (xp - this.Origin(1))/(this.Origin(2) - this.CornerPos(2)));
                    arcSh1.th1 = 2*pi - theta1;
                    arcSh1.th2 = 3*pi/2 + theta2;
                    arcSh2.th1 = theta1;
                    arcSh2.th2 = 3*pi/2 - theta2; 
            end
            
            arcSh1.R  = this.R;
            arcSh1.N  = Geometry.N;
              
            this.A1 = Arc(arcSh1);
            
            arcSh2.R  = this.R;
            arcSh2.N  = Geometry.N;
              
            this.A2 = Arc(arcSh2);
            
            ShiftArcs(this);
            
            this.Pts  = struct('y1_kv',[this.A1.Pts.y1_kv;this.A2.Pts.y1_kv],...
                               'y2_kv',[this.A1.Pts.y2_kv;this.A2.Pts.y2_kv]);

        end        
        
        function ptsCart = GetCartPts(this)            
            ptsCart = this.Pts;
        end         
        function [int,length] = ComputeIntegrationVector(this)           
            
            [A1int,lengthA1]  = this.A1.ComputeIntegrationVector(); 
            [A2int,lengthA2]  = this.A2.ComputeIntegrationVector(); 
            
            int      = [A1int,A2int];
            length   = lengthA1 + lengthA2;
            this.Int = int;

        end

        function ShiftArcs(this)

            this.A1.Pts =  Pol2CartPts(this.A1.Pts);
            this.A2.Pts =  Pol2CartPts(this.A2.Pts);

            this.A1.Pts.y1_kv = this.A1.Pts.y1_kv + this.Origin(1);
            this.A1.Pts.y2_kv = this.A1.Pts.y2_kv + this.Origin(2);
            this.A2.Pts.y1_kv = this.A2.Pts.y1_kv + this.Origin(1);
            this.A2.Pts.y2_kv = this.A2.Pts.y2_kv + this.Origin(2);


        end
    end
    
end    

