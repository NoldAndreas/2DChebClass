classdef Segment2 < Ball
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
    %and
    %      z = r sin(theta)*sin(phi)
    % 
    % we consider invariance in direction z
    % where phi in [0,pi] and     
    
    %**********************************************
    %************   Constructor   *****************
    %**********************************************
    methods 
        function this = Segment2(Geometry)
             this@Ball(Geometry);                          
        end      
    end
    
    %**********************************************
    %************   Initializations    ************
    %**********************************************
    methods (Access = public)
        function [int,area] = ComputeIntegrationVector(this)
            int = ComputeIntegrationVector@Ball(this);
            int = int.*(sin(this.Pts.y1_kv)'.*sin(this.Pts.y2_kv)'/2);
                                    
            area = this.R^2*(this.theta2-this.theta1) + ...
                   this.R^2/2*(sin(2*this.theta1) - sin(2*this.theta2));
            if(area == 0)
                disp(['Segment2: Area = 0, Error in computation (absolute):',num2str(area-sum(int))]);
            else
                disp(['Segment2: Error in computation of area (%):',num2str(1-sum(int)/area)]);
            end
            
        end                                
    end
    
end


