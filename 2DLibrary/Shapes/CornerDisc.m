classdef CornerDisc < handle
   
    properties 
        Origin = [0;0];
        R
        Corner    = 'SW';
        CornerPos = [0;0];
        S,T        
        Int
        Pts     
    end
    
    methods 
        function this = CornerDisc(Geometry)
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
            
            switch(this.Corner)
                case 'SW'
                    yi1 = this.Origin(1) + dw1;
                    yi2 = this.Origin(2) + dw2;
                    theta = -atan((yi2-this.CornerPos(2))/(yi1-this.CornerPos(1)));
                case 'NW'
                    yi1 = this.Origin(1) + dw1;
                    yi2 = this.Origin(2) - dw2;
                    theta = pi-atan((yi2-this.CornerPos(2))/(yi1-this.CornerPos(1)));
                case 'NE'
                    yi1 = this.Origin(1) - dw1;
                    yi2 = this.Origin(2) - dw2;
                    theta = pi -atan((yi2-this.CornerPos(2))/(yi1-this.CornerPos(1)));
                case 'SE'
                    yi1 = this.Origin(1) - dw1;
                    yi2 = this.Origin(2) + dw2;
                    theta = -atan((yi2-this.CornerPos(2))/(yi1-this.CornerPos(1)));
            end

            y1p = 1/2*(yi1 + this.CornerPos(1));
            y2p = 1/2*(yi2 + this.CornerPos(2));

            
            segSh.R = this.R;
            segSh.h = sqrt( (this.Origin(1) - y1p)^2 + (this.Origin(2) - y2p)^2 );
            segSh.N = Geometry.NS;
            
            this.S  = Segment(segSh);
            ShiftRotateSegment(this,segSh.h,theta,[y1p;y2p]);
            
            % Compute points for triangle:
            triSh.Y = [this.CornerPos(1), this.CornerPos(2);
                       this.CornerPos(1), yi2;
                       yi1, this.CornerPos(2)];
                   
            this.T = Triangle(triSh,Geometry.NT);
            
            this.Pts  = struct('y1_kv',[this.S.Pts.y1_kv;this.T.Pts.y1_kv],...
                               'y2_kv',[this.S.Pts.y2_kv;this.T.Pts.y2_kv]);                            
        end        
        function ptsCart = GetCartPts(this)            
            ptsCart = this.Pts;
        end         
        function [int,area] = ComputeIntegrationVector(this)           
            
            [Tint,areaT]  = this.T.ComputeIntegrationVector(); 
            [Sint,areaS]  = this.S.ComputeIntegrationVector();
            
            int  = [Sint,Tint];
            area = areaT + areaS;
            this.Int = int;            

        end
        function ShiftRotateSegment(this,h,th,yp)
            
            y1_kv = this.S.Pts.y1_kv;
            y2_kv = this.S.Pts.y2_kv - h;

            y1_kvs = cos(th)*y1_kv - sin(th)*y2_kv + yp(1);
            y2_kvs = sin(th)*y1_kv + cos(th)*y2_kv + yp(2);

            this.S.Pts.y1_kv = y1_kvs;
            this.S.Pts.y2_kv = y2_kvs;

        end
    end
    
end    

