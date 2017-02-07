classdef Sphere < SpectralSpectral%SpectralEvenFourier%
    %Computes the integral over the surface of a sphere, where it is
    %assumed that the x-y-plane is a axis of symmetry for the function to
    %evaluate. (otherwise inherit from SpectralFourier, not
    %SpectralEvenFourier)
    %
    %here, we use spherical coordinates and fix the radial variable to a
    %cosntant. We then consider the projection onto the following plane:
    %
    %here, x = r sin(theta)*cos(phi)
    %      y = r cos(theta)
    % 
    % we consider invariance in direction z = r sin(theta)*sin(phi)
    % where phi in [0,pi] and 
    
    properties 
        theta1 = 0
        theta2 = pi
        R
        PtsCart    
        volume = true
    end
    
    %**********************************************
    %************   Constructor   *****************
    %**********************************************
    methods 
        function this = Sphere(Geometry)
             this@SpectralSpectral(Geometry.N(1),Geometry.N(2)); %EvenFourier
             
             if(isfield(Geometry,'theta1'))
                  this.theta1 = Geometry.theta1;
             end
             if(isfield(Geometry,'theta2'))
                this.theta2 = Geometry.theta2;
             end
             this.R      = Geometry.R;
             if(isfield(Geometry,'Origin'))
                 this.Origin = Geometry.Origin;
             end
             if(isfield(Geometry,'volume'))
                 this.volume = Geometry.volume;
             end
             
             InitializationPts(this);             
             this.polar = 'sphSurf';             
        end      
    end
    
    %**********************************************
    %************   Initializations    ************
    %**********************************************
    methods (Access = public) 
        function [int,area] = ComputeIntegrationVector(this)
            R   = this.R;
            y1  = this.Pts.y1_kv;
            y2  = this.Pts.y2_kv;
            th1 = this.theta1;  
            th2 = this.theta2;
            
            int = ComputeIntegrationVector@SpectralSpectral(this);            %EvenFourier
            %int = int.*(2*this.R^2*sin(this.Pts.y1_kv)'); %weight from spherical coordinates
            int = int.*(R^2*sin(y2)'.*sin(y1').^2); %weight from spherical coordinates            
                           
           
            if(this.volume)
                int  = int.*(2*R*sin(y1).*sin(y2))';                                
                %area = 4/3*pi*R^3*( 1 - ((1-cos(th1)).^2).*(2+cos(th1))/4);
                area = pi*R^3/6*( (cos(th2)*(cos(2*th2)-5)) - ...
                                  (cos(th1)*(cos(2*th1)-5)));
            else
                %area = R^2*(pi-(2*th1-sin(2*th1))/2);%this.R^2*pi*(cos(this.theta1)-cos(this.theta2));
                area = R^2/2*(2*(th2-th1)-(sin(2*th2)-sin(2*th1))); %this.R^2*pi*(cos(this.theta1)-cos(this.theta2));
                %area = R^2*(3/4*pi+1/2);
            end
            
            if(nargout < 2)
                if(area == 0)
                    disp(['Ball: Area = 0, Error in computation (absolute):',num2str(area-sum(int))]);
                else
                    disp(['Ball: Error in computation of area (ratio):',num2str(1-sum(int)/area)]);
                end
            end
            
            
            this.Int = int;
            
        end                  
        function ComputeConvolutionMatrix(this)
            exc = MException('Segment2:ComputeConvolutionMatrix','case not implemented');
            throw(exc);           
        end       
        function Ind    = ComputeIndices(this)
            exc = MException('Segment2:ComputeIndices','case not implemented');
            throw(exc);           
        end           
        function Diff   = ComputeDifferentiationMatrix(this)
            exc = MException('Segment2:ComputeDifferentiationMatrix','case not implemented');
            throw(exc);           
        end               
           
        %**************
        %**** Maps ****
        %**************
        function [y,dy,dx,ddx,dddx,ddddx] = PhysSpace1(this,x)            
            [y,dy,dx,ddx,dddx,ddddx] = LinearMap(x,this.theta1,this.theta2);            
        end
        function xf = CompSpace1(this,y)
            xf  = InvLinearMap(y,this.theta1,this.theta2);           
        end
        function [y,dy,dx,ddx,dddx,ddddx] = PhysSpace2(this,x)    
            [y,dy,dx,ddx,dddx,ddddx] = LinearMap(x,0,pi);%CCC
        end
        function xf = CompSpace2(this,y)
            xf = InvLinearMap(y,0,pi);%y/(2*pi);
        end                                        
        
        acc = testIntegration(this,toCheck);
    end
    
end


