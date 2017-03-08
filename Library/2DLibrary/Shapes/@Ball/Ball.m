classdef Ball < SpectralEvenFourier
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
        theta1,theta2,R
        PtsCart    
    end
    
    %**********************************************
    %************   Constructor   *****************
    %**********************************************
    methods 
        function this = Ball(Geometry)
             this@SpectralEvenFourier(Geometry.N(1),Geometry.N(2));
             
             this.theta1 = Geometry.theta1;
             this.theta2 = Geometry.theta2;
             this.R      = Geometry.R;
             if(isfield(Geometry,'Origin'))
                 this.Origin = Geometry.Origin;
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
            int = ComputeIntegrationVector@SpectralEvenFourier(this);            
            int = int.*(2*this.R^2*sin(this.Pts.y1_kv)'); %weight from spherical coordinates
            %int = int.*(this.R^2*sin(this.Pts.y1_kv)'.*sin(this.Pts.y2_kv').^2); %weight from spherical coordinates
            
            this.Int = int;
            
            area = 2*this.R^2*pi*(cos(this.theta1)-cos(this.theta2));
            if(nargout < 2)
                if(area == 0)
                    disp(['Ball: Area = 0, Error in computation (absolute):',num2str(area-sum(int))]);
                else
                    disp(['Ball: Error in computation of area (%):',num2str(1-sum(int)/area)]);
                end
            end
            
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
        function [y1,dy1,dx,ddx,dddx,ddddx] = PhysSpace1(this,x1)            
            [y1,dy1,dx,ddx,dddx,ddddx] = LinearMap(x1,this.theta1,this.theta2);
        end
        function xf = CompSpace1(this,y1)
            xf  = InvLinearMap(y1,this.theta1,this.theta2);
        end
        function [th,dth,dx,ddx,dddx,ddddx] = PhysSpace2(this,xT)    
            [th,dth,dx,ddx,dddx,ddddx] = LinearMap01(xT,0,2*pi);%CCC
        end
        function xf = CompSpace2(this,th)
            xf = th/(2*pi);
        end                                        
        
        acc = testIntegration(this,toCheck);
    end
    
end


