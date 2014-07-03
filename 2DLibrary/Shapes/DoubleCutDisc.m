classdef DoubleCutDisc < handle
   
    properties 
        Origin = [0;0];
        R
        Corner    = 'SW';
        CornerPos = [0;0];
        W1,W2,T1,T2
        Int
        Pts     
    end
    
    methods 
        function this = DoubleCutDisc(Geometry)
            if(isfield(Geometry,'Origin'))
                this.Origin     = Geometry.Origin;                
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

            theta1 = atan( dw2/(this.Origin(1) - this.CornerPos(1)));
            theta2 = atan( dw1/(this.Origin(2) - this.CornerPos(2)));
                    
            switch(this.Corner)
                case 'SW'
                    wedgeSh1.th1 = pi + theta1;
                    wedgeSh1.th2 = 3*pi/2 - theta2;
                    wedgeSh2.th1 = -pi/2 + theta2;
                    wedgeSh2.th2 =  pi - theta1;
                case 'NW'                    
                    theta2 = -theta2;
                    wedgeSh1.th1 = pi/2 + theta2;
                    wedgeSh1.th2 = pi - theta1;
                    wedgeSh2.th1 = -pi + theta1; 
                    wedgeSh2.th2 = pi/2 - theta2;
                case 'NE'            
                    theta1 = -theta1;
                    theta2 = -theta2;                    
                    wedgeSh1.th1 = theta1;
                    wedgeSh1.th2 = pi/2 - theta2;
                    wedgeSh2.th1 = pi/2 + theta2;
                    wedgeSh2.th2 = 2*pi-theta1; 
                case 'SE'
                    theta1 = -theta1;
                    wedgeSh1.th1 = 3*pi/2 + theta2;
                    wedgeSh1.th2 = 2*pi - theta1;
                    wedgeSh2.th1 = theta1;
                    wedgeSh2.th2 = 3*pi/2 - theta2; 
            end
            % Compute points for triangles:
            triSh1.Y = [0, 0;
                        this.CornerPos(1) - this.Origin(1), dw2;
                        this.CornerPos(1) - this.Origin(1), -dw2];
                   
            this.T1 = Triangle(triSh1,Geometry.NT);

            triSh2.Y = [0, 0;
                        dw1 , this.CornerPos(2) - this.Origin(2);
                        -dw1, this.CornerPos(2) - this.Origin(2)];
                   
            this.T2 = Triangle(triSh2,Geometry.NT);
            
            wedgeSh1.R  = this.R;
            wedgeSh1.N  = Geometry.NW;
              
            this.W1 = Wedge(wedgeSh1);
            
            wedgeSh2.R  = this.R;
            wedgeSh2.N  = Geometry.NW;
              
            this.W2 = Wedge(wedgeSh2);           
            
            this.Pts  = struct('y1_kv',[this.T1.GetCartPts.y1_kv;this.T2.GetCartPts.y1_kv; ...
                                        this.W1.GetCartPts.y1_kv;this.W2.GetCartPts.y1_kv],...
                               'y2_kv',[this.T1.GetCartPts.y2_kv;this.T2.GetCartPts.y2_kv; ...
                                        this.W1.GetCartPts.y2_kv;this.W2.GetCartPts.y2_kv]);
        end        
        
        function ptsCart = GetCartPts(this)            
            ptsCart       = this.Pts;
            ptsCart.y1_kv = ptsCart.y1_kv + this.Origin(1);
            ptsCart.y2_kv = ptsCart.y2_kv + this.Origin(2);
        end         
        function [int,area] = ComputeIntegrationVector(this)           
            
            [T1int,area1]  = this.T1.ComputeIntegrationVector(); 
            [T2int,area2]  = this.T2.ComputeIntegrationVector(); 
            [W1int,area3]  = this.W1.ComputeIntegrationVector();
            [W2int,area4]  = this.W2.ComputeIntegrationVector();
            
            int  = [T1int,T2int,W1int,W2int];
            area = area1 + area2 + area3 + area4;
            
            this.Int = int;

        end

    end
    
end    

