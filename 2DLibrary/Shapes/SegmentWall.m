classdef SegmentWall < Segment
    
    properties
        alpha = 0
    end
    
    methods
        function this = SegmentWall(Geometry)
            Geometry.h = strcmp(Geometry.Wall_VertHor,'horizontal')*...
                            abs(Geometry.Wall_Y - Geometry.Origin(2)) + ... 
                         strcmp(Geometry.Wall_VertHor,'vertical')*...
                            abs(Geometry.Wall_Y - Geometry.Origin(1));
            this@Segment(Geometry);

            if(strcmp(Geometry.Wall_VertHor,'horizontal') && (Geometry.Wall_Y > Geometry.Origin(2)))
                this.alpha = 0;
            elseif(strcmp(Geometry.Wall_VertHor,'horizontal') && (Geometry.Wall_Y <= Geometry.Origin(2)))
                this.alpha = -pi;
            elseif(strcmp(Geometry.Wall_VertHor,'vertical') && (Geometry.Wall_Y > Geometry.Origin(1)))
                this.alpha = -pi/2;
            elseif(strcmp(Geometry.Wall_VertHor,'vertical') && (Geometry.Wall_Y <= Geometry.Origin(1)))
                this.alpha = pi/2;                
            end
            
            if(isfield(Geometry,'alpha'))
                this.alpha = Geometry.alpha;
            end
            
            InitializationPts(this);
            
        end
        function [y1_kv,y2_kv,J] = PhysSpace(this,x1,x2)        
            
            y10 = this.Origin(1);
            y20 = this.Origin(2);
            
            [y1_kv,y2_kv,J] = PhysSpace@Segment(this,x1,x2);

            [y1R,y2R,JR] = RotateAroundOrigin(y1_kv - y10,...
                                           y2_kv - y20,this.alpha);
            y1_kv     = y1R + y10;
            y2_kv     = y2R + y20;
            J         = MultMatr(JR,J);
        end
    end
    
    
end