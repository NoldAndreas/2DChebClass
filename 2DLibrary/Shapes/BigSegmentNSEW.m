classdef BigSegmentNSEW < handle
   
    properties 
        Origin = [0;0];
        R,h
        WallPos = 'S';
        % h positive = inside wall at position WallPos
        % h negative = outside wall at position WallPos
        W,T        
        Int
        Pts     
        sphere
    end
    
    methods 
        function this = BigSegmentNSEW(Geometry)
            if(isfield(Geometry,'Origin'))
                this.Origin = Geometry.Origin;
            end
            if(isfield(Geometry,'WallPos'))
                this.WallPos = Geometry.WallPos;
            end
            if(isfield(Geometry,'sphere'))
                this.sphere    = Geometry.sphere;
                wedgeSh.sphere = Geometry.sphere;
            end
            this.R      = Geometry.R;            
            this.h      = Geometry.h;
                        
            if(this.h<0)
                switch this.WallPos
                    case 'N'
                        this.WallPos = 'S';
                    case 'S'
                        this.WallPos = 'N';
                    case 'E'
                        this.WallPos = 'W';
                    case 'W'
                        this.WallPos = 'E';
                end
                this.h = -this.h;
            end
  
            th = acos(-this.h/this.R);
            
            switch this.WallPos
                case 'N'
                    wedgeSh.th1   = 3/2*pi-th;            
                    wedgeSh.th2   = 3/2*pi+th;
                    c1 = [this.Origin(1)-sin(th)*this.R, ...
                          this.Origin(2)+this.h];
                    c2 = [this.Origin(1)+sin(th)*this.R, ...
                          this.Origin(2)+this.h];
                case 'S'
                    wedgeSh.th1   = pi/2+th;
                    wedgeSh.th2   = pi/2-th;
                    c1 = [this.Origin(1)-sin(th)*this.R, ...
                          this.Origin(2)-this.h];
                    c2 = [this.Origin(1)+sin(th)*this.R, ...
                          this.Origin(2)-this.h];
                case 'E'
                    wedgeSh.th1   = pi+th;
                    wedgeSh.th2   = pi-th;
                    c1 = [this.Origin(1)+this.h, ...
                          this.Origin(2)-sin(th)*this.R];
                    c2 = [this.Origin(1)+this.h, ...
                          this.Origin(2)+sin(th)*this.R];
                case 'W'
                    wedgeSh.th1   = -th;
                    wedgeSh.th2   = th;
                    c1 = [this.Origin(1)-this.h, ...
                          this.Origin(2)-sin(th)*this.R];
                    c2 = [this.Origin(1)-this.h, ...
                          this.Origin(2)+sin(th)*this.R];
            end
            
            wedgeSh.R_out = this.R;
            wedgeSh.N     = Geometry.NW;

            this.W       = Wedge(wedgeSh);                   
            [y1c,y2c]    = WedgeToCart(this);

            Y = [this.Origin(1), this.Origin(2);
                c1;
                c2];

            this.T = Triangle(v2struct(Y),Geometry.NT);
            
            this.Pts  = struct('y1_kv',[y1c;this.T.Pts.y1_kv],...
                               'y2_kv',[y2c;this.T.Pts.y2_kv]);                            
        end        
        function ptsCart = GetCartPts(this)            
            ptsCart = this.Pts;
        end         
        function [int,area] = ComputeIntegrationVector(this)           
            
            [Tint,h1s] = this.T.ComputeIntegrationVector();                        
            
            if(this.sphere)                                
                [wint,h]  = this.W.ComputeIntegrationVector(true,false);
                int       = [wint,Tint];            
                
                y1s  = this.Pts.y1_kv - this.Origin(1);                
                y2s  = this.Pts.y2_kv - this.Origin(2);
                int  = 2*int.*real(sqrt(this.R^2-y1s.^2-y2s.^2))';  
                this.Int = int;
                
                ht   = this.R - abs(this.h);
                
                th = acos(abs(this.h)/this.R);
                area1 = 4/3*(pi-th)*this.R^3;
                intW  = int(1:length(this.W.Int));
                
                area = 4/3*pi*this.R^3 - pi*ht^2/3*(3*this.R - ht);                
                
                area2 = area - area1;
                intT = int(length(this.W.Int)+1:end);
                
                if(nargout < 2)
                    disp(['Error of integration of (spherical) Wedge (ratio):',num2str(1-abs(sum(intW)/area1))]);                
                    disp(['Error of integration of (spherical) Wedge (abs):',num2str(abs(sum(intW)-area1))]);                

                    disp(['Error of integration of (spherical) Triangle (ratio):',num2str(1-abs(sum(intT))/area2)]);                
                    disp(['Error of integration of (spherical) Triangle (abs):',num2str(abs(sum(intT)-area2))]);                
                end
            else
                [wint,h]  = this.W.ComputeIntegrationVector();
                int       = [wint,Tint];
                %area      = area1 + area2;
                hh   = abs(this.h);
                th   = pi - acos(hh/this.R);
                area = this.R^2*(th-0) + ...
                        this.R^2/2*(sin(2*0) - sin(2*th));                
                this.Int  = int;
            end
            if(nargout < 2)
                if(area ~= 0)
                    disp(['BigSegment: Error of integration of area (ratio): ',...
                                        num2str(1-sum(this.Int)/area)]);   
                end
                disp(['BigSegment: Error of integration of area (absolute,area = 0): ',...
                                        num2str(area-sum(this.Int))]);   
            end
        end
        function [y1c,y2c]= WedgeToCart(this)
            [x,y] = pol2cart(this.W.Pts.y2_kv,this.W.Pts.y1_kv);
            
            y1c = x + this.Origin(1);
            y2c = y + this.Origin(2);    
        end
    end
    
end    

