classdef InfCapillary_FMT < InfCapillary 
    properties
        R
        AD        
    end
    
    methods
        
        %% Exemple for grids such as defined by HalfSpace_FMT
        % 
        % <<GridsDensities.PNG>>
        % 
        % *Caption:* the grid depicted on the left is defined by _this_ object and its parents, and the
        % grid on the right by the property object _AD_ of this class         
        
        function this = InfCapillary_FMT(Geometry,R)   
            
            if(nargin < 2)
                R    = Geometry.R;
            end
            
            % multiply with factor to take into account skewed grid by
            % angle alpha
            
            alpha_rad = pi/2;
            L1        = Geometry.L1/sin(alpha_rad);
%            L2        = Geometry.L2/sin(alpha_rad);           
%            L2_AD     = Geometry.L2_AD/sin(alpha_rad);            
            
            shapeIC       = Geometry;            
            shapeIC.L1    = L1;
%            shapeIC.L2    = L2;
            
            this@InfCapillary(shapeIC);
                        
            this.R        = R;            
            
            %Geometry includes:
           % y2Min,h,L1,L2,N
            %Check how to handle multiple Species!
            %this.Boundary(1,length(R)) = InfCapillary();
            %for i = 1:length(R)            
            
            %                             'L2',L2_AD,...
            shapeAD = struct('L1',L1,...
                             'y2Min',Geometry.y2Min,...
                             'y2Max',Geometry.y2Max,...
                             'R',R,...
                             'N',Geometry.N,'N2bound',Geometry.N2bound);

            this.AD = InfCapillary_Composed(shapeAD);                                                 
        end                        
                
        function [AD,AAD] = GetAverageDensities(this,area,weights)            
            %%
            % 
            % This function is designed to compute the averaged densities and the
            % functional derivates of the free energy in for Fundamental Measure Theory (FMT).
            % With AD, it computes the convolution of a function defined on a HalfSpace,
            % with a list of weight functions of finite support. The support of the resulting function 
            % is a half space, which goes beyond the half space of the original
            % function. With AAD, the convolution of a function defined on this
            % extended halfspace is computed with the same list of weight functions defined in the input.
            %
            %% Input
            %
            % * area - structure with two methods (1) [int,A] = ComputeIntegrationVector() and
            %                   (2) pts = GetCartPts(), with pts = struct(y1_kv,y2_kv)
            % * weights - a cell with a list of functions which get a structure with 
            %         two arrays [y1_kv,y2_kv] as input and returns one
            %         array of same size. [y1_kv,y2_kv] represents a point in polar
            %         coordinates, representing the radial and angular
            %         component, respectively.
            %
            %% Output
            %
            % Both AD and AAD are operators for the following convolution:
            %
            % $$X_{i,:,k}\cdot \rho = \int_{A}  \rho({\bf r}_i+{\bf r}')w_k({\bf r}') d{\bf r}'$$
            %
            %
            % * AD  : $A = area.GetCartPts() \cap this.GetCartPts()$, and ${\bf r}_i \in this.AD.GetCartPts()$
            %    (average density)
            %
            % * AAD  : $A = area.GetCartPts() \cap this.AD.GetCartPts()$, and ${\bf r}_i \in this.GetCartPts()$
            %       (average the average densities to compute free energy)
            %            
            AAD  = this.AD.ComputeConvolutionFiniteSupport(area,weights,this.Pts);
            AD   = ComputeConvolutionFiniteSupport(this,area,weights,this.AD.Pts);  
        end
    end
    
    
end