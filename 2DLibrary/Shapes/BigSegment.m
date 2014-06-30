classdef BigSegment < handle
   
    properties 
        Origin = [0;0];        
        R,h
        W,T        
        Int
        Pts     
        sphere
    end
    
    methods 
        function this = BigSegment(Geometry)
            if(isfield(Geometry,'Origin'))
                this.Origin = Geometry.Origin;
            end
            if(isfield(Geometry,'sphere'))
                this.sphere    = Geometry.sphere;
                wedgeSh.sphere = Geometry.sphere;
            end
            this.R      = Geometry.R;   
            
            %input th1, th2 
            if(isfield(Geometry,'Wall_Y'))        
                if(strcmp(Geometry.Wall_VertHor,'horizontal'))
                    h = Geometry.Wall_Y - Geometry.Origin(2);
                    
                    Geometry.th1 = pi - asin(h/this.R);              
                    Geometry.th2 = asin(h/this.R);                    
                else
                    h = Geometry.Wall_Y - Geometry.Origin(1);
                    
                    Geometry.th1 = acos(h/this.R);              
                    Geometry.th2 = -acos(h/this.R);                          
                end                                                
            end                        
            
            
            r   = Geometry.R;                                        
            
            Geometry.th1 = mod(Geometry.th1,2*pi);
            Geometry.th2 = mod(Geometry.th2,2*pi);
            
            if(mod(Geometry.th1 - Geometry.th2,2*pi) < pi)
                th2 = Geometry.th2;
                th1 = Geometry.th1;
            else
                th1 = Geometry.th2;
                th2 = Geometry.th1;
            end
            
            
            wedgeSh.th1   = th1;
            wedgeSh.th2   = th2;                          
            wedgeSh.R_out = this.R;
            wedgeSh.N     = Geometry.NW;

            this.W       = Wedge(wedgeSh);                   
            [y1c,y2c]    = WedgeToCart(this);
            
            %Compute points for triangle:
            Y = [this.Origin(1),this.Origin(2);
                 this.Origin(1)+cos(th1)*r,this.Origin(2)+sin(th1)*r;
                 this.Origin(1)+cos(th2)*r,this.Origin(2)+sin(th2)*r];
            
            this.T = Triangle(v2struct(Y),Geometry.NT);
            
            this.Pts  = struct('y1_kv',[y1c;this.T.Pts.y1_kv],...
                               'y2_kv',[y2c;this.T.Pts.y2_kv]);                            
                           
                           
            this.h = r*abs(cos((th1-th2)/2));
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

