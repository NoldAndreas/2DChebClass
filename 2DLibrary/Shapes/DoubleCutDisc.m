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

            xm = this.Origin(1) - dw1;
            xp = this.Origin(1) + dw1;
            ym = this.Origin(2) - dw2;
            yp = this.Origin(2) + dw2;

            
            switch(this.Corner)
                case 'SW'
                    theta1 = atan( (yp - this.Origin(2))/(this.Origin(1) - this.CornerPos(1)));
                    theta2 = atan( (xp - this.Origin(1))/(this.Origin(2) - this.CornerPos(2)));
                    wedgeSh1.th1 = pi + theta1;
                    wedgeSh1.th2 = 3*pi/2 - theta2;
                    wedgeSh2.th1 = pi - theta1;
                    wedgeSh2.th2 = -pi/2 + theta2; 
                case 'NW'
                    theta1 = atan( (yp - this.Origin(2))/(this.Origin(1) - this.CornerPos(1)));
                    theta2 = atan( (xp - this.Origin(1))/(this.CornerPos(2) - this.Origin(2)));
                    wedgeSh1.th1 = pi/2 + theta2;
                    wedgeSh1.th2 = pi - theta1;
                    wedgeSh2.th1 = pi/2 - theta2;
                    wedgeSh2.th2 = -pi + theta1; 
                case 'NE'                    
                    theta1 = atan( (yp - this.Origin(2))/(this.CornerPos(1) - this.Origin(1)));
                    theta2 = atan( (xp - this.Origin(1))/(this.CornerPos(2) - this.Origin(2)));
                    wedgeSh1.th1 = pi/2 - theta2;
                    wedgeSh1.th2 = theta1;
                    wedgeSh2.th1 = 2*pi-theta1;
                    wedgeSh2.th2 = pi/2 + theta2; 
                case 'SE'
                    theta1 = atan( (yp - this.Origin(2))/(this.CornerPos(1) - this.Origin(1)));
                    theta2 = atan( (xp - this.Origin(1))/(this.Origin(2) - this.CornerPos(2)));
                    wedgeSh1.th1 = 2*pi - theta1;
                    wedgeSh1.th2 = 3*pi/2 + theta2;
                    wedgeSh2.th1 = theta1;
                    wedgeSh2.th2 = 3*pi/2 - theta2; 
            end

            % Compute points for triangles:
            triSh1.Y = [this.Origin(1), this.Origin(2);
                        this.CornerPos(1), yp;
                        this.CornerPos(1), ym];
                   
            this.T1 = Triangle(triSh1,Geometry.NT);

            triSh2.Y = [this.Origin(1), this.Origin(2);
                        xp, this.CornerPos(2);
                        xm, this.CornerPos(2)];
                   
            this.T2 = Triangle(triSh2,Geometry.NT);
            
            wedgeSh1.R  = this.R;
            wedgeSh1.N  = Geometry.NW;
              
            this.W1 = Wedge(wedgeSh1);
            
            wedgeSh2.R  = this.R;
            wedgeSh2.N  = Geometry.NW;
              
            this.W2 = Wedge(wedgeSh2);
            
            ShiftWedges(this);
            
            this.Pts  = struct('y1_kv',[this.T1.Pts.y1_kv;this.T2.Pts.y1_kv; ...
                                        this.W1.Pts.y1_kv;this.W2.Pts.y1_kv],...
                               'y2_kv',[this.T1.Pts.y2_kv;this.T2.Pts.y2_kv; ...
                                        this.W1.Pts.y2_kv;this.W2.Pts.y2_kv]);

%             this.Pts  = struct('y1_kv',[this.T1.Pts.y1_kv;this.T2.Pts.y1_kv],...
%                    'y2_kv',[this.T1.Pts.y2_kv;this.T2.Pts.y2_kv]);

%             this.Pts  = struct('y1_kv',[this.T1.Pts.y1_kv;this.T2.Pts.y1_kv; ...
%                                         this.W1.Pts.y1_kv;],...
%                                'y2_kv',[this.T1.Pts.y2_kv;this.T2.Pts.y2_kv; ...
%                                         this.W1.Pts.y2_kv;]);

%             this.Pts  = struct('y1_kv',[this.W1.Pts.y1_kv;this.W2.Pts.y1_kv],...
%                    'y2_kv',[this.W1.Pts.y2_kv;this.W2.Pts.y2_kv]);


        end        
        
        function ptsCart = GetCartPts(this)            
            ptsCart = this.Pts;
        end         
        function int = ComputeIntegrationVector(this)           
            
            [T1int,~]  = this.T1.ComputeIntegrationVector(); 
            [T2int,~]  = this.T2.ComputeIntegrationVector(); 
            [W1int,~]  = this.W1.ComputeIntegrationVector();
            [W2int,~]  = this.W2.ComputeIntegrationVector();
            
            int  = [T1int,T2int,W1int,W2int];
            this.Int = int;

        end

        function ShiftWedges(this)
            [x1,y1] = pol2cart(this.W1.Pts.y2_kv,this.W1.Pts.y1_kv);
            this.W1.Pts.y1_kv = x1 + this.Origin(1);
            this.W1.Pts.y2_kv = y1 + this.Origin(2);
            
            [x2,y2] = pol2cart(this.W2.Pts.y2_kv,this.W2.Pts.y1_kv);
            this.W2.Pts.y1_kv = x2 + this.Origin(1);
            this.W2.Pts.y2_kv = y2 + this.Origin(2);

        end
    end
    
end    

