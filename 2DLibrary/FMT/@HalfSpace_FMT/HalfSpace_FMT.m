classdef HalfSpace_FMT < HalfSpaceSkewed & ConvolutionFiniteSupport
    properties
        R
        AD        
        y2wall
    end
    
    methods
        
        %% Exemple for grids such as defined by HalfSpace_FMT
        % 
        % <<GridsDensities.PNG>>
        % 
        % *Caption:* the grid depicted on the left is defined by _this_ object and its parents, and the
        % grid on the right by the property object _AD_ of this class         
        
        function this = HalfSpace_FMT(Geometry,R)   
            
            if(~isfield(Geometry,'alpha') && ~isfield(Geometry,'alpha_deg'))
                Geometry.alpha = pi/2; %default: 90 deg
                Geometry.alpha_deg = 90; %default: 90 deg
            elseif(~isfield(Geometry,'alpha') && isfield(Geometry,'alpha_deg'))                
                Geometry.alpha = Geometry.alpha_deg*pi/180;                            
            elseif(isfield(Geometry,'alpha') && ~isfield(Geometry,'alpha_deg'))                
                Geometry.alpha_deg = Geometry.alpha*180/pi;
            end                       
            if(nargin < 2)
                R    = Geometry.R;
            end
            
            % multiply with factor to take into account skewed grid by
            % angle alpha
            
            alpha_rad = Geometry.alpha_deg*pi/180;
            L1        = Geometry.L1/sin(alpha_rad);
            L2        = Geometry.L2/sin(alpha_rad);
            L2_AD     = Geometry.L2_AD/sin(alpha_rad);
            %config.optsNum.PhysArea.L2        = 2/sin(config.optsNum.PhysArea.alpha_deg*pi/180);
            
            shapeHS       = Geometry;
            shapeHS.y2Min = Geometry.y2wall+R;
            shapeHS.L1    = L1;
            shapeHS.L2    = L2;
            
            this@HalfSpaceSkewed(shapeHS);
                        
            this.R        = R;
            this.y2wall   = Geometry.y2wall;            
            
            %Geometry includes:
           % y2Min,h,L1,L2,N
            %Check how to handle multiple Species!
            %this.Boundary(1,length(R)) = InfCapillary();
            %for i = 1:length(R)            
            shapeAD = struct('L1',L1,...
                             'L2',L2_AD,...
                             'y2Min',Geometry.y2wall,...
                             'h',Geometry.h,...
                             'N',Geometry.N,'N2bound',Geometry.N2bound,...
                             'alpha',Geometry.alpha);

            this.AD = HalfSpace_Composed(shapeAD);                        
            
             %There are four numberings: 
             % (here, there is only one domain,simplifying the numbering, as
             % all identifiers = 1)
             %(1) the numberings of the domains. The points with identifyer
             %     X are given by mark_id(:,X)       
             this.mark_id = true(this.M,1);
             %(2) the numbering of the domains in the first dimension. Points
             %in the first dimension belonging to the part with identifyer Y
             %are given by mark_id_1(:,Y)
             this.mark_id_1 = true(this.N1,1);
             %(3) the numbering of the domains in the second dimension. Points
             %in the second dimension belonging to the part with identifyer Z
             %are given by mark_id_2(:,Z)                      
             this.mark_id_2 = true(this.N2,1);
             %mark12 makes a link between identifyers of dimensions 1,2 and 
             % the global identifyer. Points which are in the first dimension
             % located in a domain with identifyer Y and in the second
             % direction in a domain with identifyer Z, belong to the global
             % domain with the identifyer X = mark_12(Y,Z).                 
             this.mark_12   = ones(1,1); 
        end                  
        function do1Dplot(this,val,y1)  
            global PersonalUserOutput
            if(~PersonalUserOutput)
                return;
            end
            
            %TODO: Do proper interpolation
            if(nargin == 2)
                y1 = inf;
            end
            
            %subPtsAAD.y2_kv           = (0.5:0.01:3.5)';
            %subPtsAAD.y1_kv           = inf*ones(size(subPtsAAD.y2_kv));           
            %IP                        = SubShapePts(subPtsAAD);
            
            y2_kv_1D                  = this.Pts.y2_kv(this.Pts.y1_kv == y1);
            y2_1D                     = this.Interp.pts2(Interp.pts1 == y1);
            Interp1D_IP               = this.Interp.InterPol(Interp.pts1 == y1,this.Pts.y1_kv == y1);
            
            %Interp1D_AAD.InterPol     = IP(:,this.Pts.y1_kv == inf);
            %Interp1D_AAD.pts1         = subPtsAAD.y1_kv; 
            %Interp1D_AAD.pts2         = subPtsAAD.y2_kv;
            
            plot(y2_kv_1D,val,'o','MarkerEdgeColor','k','MarkerFaceColor','g'); hold on
            plot(y2_1D,Interp1D_IP.InterPol*val,'linewidth',1.5);
            xlim([min(Interp1D_IP.pts2) max(Interp1D_IP.pts2)]);
            h = xlabel('$y$');     set(h,'Interpreter','Latex'); set(h,'fontsize',25);
        	%h = ylabel('$\rho$');  set(h,'Interpreter','Latex'); set(h,'fontsize',25);                        
            %title(sel(i));    
            set(gca,'fontsize',20);                        
            set(gca,'linewidth',1.5);                                    
        end
        
        AD = ComputeConvolutionFiniteSupport(this,area,weights,pts,params)
        
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

        M_int  = ComputeAreaIntegrationFiniteSupport(this,areaGeom,f,params,convolution)

    end
    
    
end