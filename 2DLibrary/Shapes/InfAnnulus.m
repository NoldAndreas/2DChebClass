classdef InfAnnulus < Polar_SpectralFourierNoOrigin

    properties        
        RMin
        L        
    end
    methods        
        function this = InfAnnulus(Geometry)
            this@Polar_SpectralFourierNoOrigin(Geometry.N(1),Geometry.N(2));
            
            this.RMin = Geometry.RMin; 
            
            if(isfield(Geometry,'L'))
                this.L    = Geometry.L; 
            elseif(isfield(Geometry,'f'))
                this.L = FindOptimalL(@RadMap,this.Pts.x1,@Geometry.f,10);
            else
                error('InfAnnulus: No L or f given in input.');
            end
            
            if(isfield(Geometry,'Origin'))
                this.Origin = Geometry.Origin;
            end
            
            InitializationPts(this);   
            
            function r = RadMap(x,L)
                r = QuotientMap(x,L,this.RMin,inf);
            end
        end                         
    end
    
    methods (Access = public)     
        function [r,dr,dx,ddx,dddx,ddddx] = PhysSpace1(this,x)
            [r,dr,dx,ddx,dddx,ddddx] = QuotientMap(x,this.L,this.RMin,inf);
        end 
        function xf = CompSpace1(this,r)
            xf = InvQuotientMap(r,this.L,this.RMin,inf);
        end
        function [th,dth,dx,ddx,dddx,ddddx] = PhysSpace2(this,xT)    
            [th,dth,dx,ddx,dddx,ddddx] = LinearMap01(xT,0,2*pi);
        end
        function xf = CompSpace2(this,th)
            xf = th/(2*pi);
        end                
        function [int,area] = ComputeIntegrationVector(this)
            int  = ComputeIntegrationVector@Polar_SpectralFourierNoOrigin(this);            
            area = Inf;
            %Check Accuracy
%             r = this.Pts.y1_kv;
%             t = this.Pts.y2_kv;
%             
%             % f = exp(-r^2); 
%             intF1 = pi*exp(-this.RMin^2);
%             f1 = exp(-r.^2);
%             disp(['InfAnnulus: Error of integration of exp(-r^2) (ratio): ',...
%                                     num2str(1-sum(this.Int*f1)/intF1)]); 
%                                 
%             % f = exp(-r^2)*sin(t); 
%             intF2 = 0;
%             f2 = exp(-r.^2).*sin(t);
%             disp(['InfAnnulus: Error of integration of exp(-r^2) sin(t) (ratio): ',...
%                                     num2str(intF2-sum(this.Int*f2))]); 
%                                 
%             % f = exp(-r^2)*sin(t); 
%             intF3 = pi/2*exp(-this.RMin^2);
%             f3 = exp(-r.^2).*sin(t).^2;
%             disp(['InfAnnulus: Error of integration of exp(-r^2) sin^2(t)(ratio): ',...
%                                     num2str(1-sum(this.Int*f3)/intF3)]);
            
        end               
        function Diff = ComputeDifferentiationMatrix(this)
            Diff = ComputeDifferentiationMatrix@Polar_SpectralFourierNoOrigin(this);
            
            %Check Accuracy

            % grad f = rhat \partial_r f + that 1/r \partial_t f
            % div V  = 1/r \partial_r (r V_r) + 1\r \partial_t V_t
            % Lap    = \partial_r^2 + 1/r \partial_r + 1/r^2 \partial_t^2 
            
            r = this.Pts.y1_kv;
            t = this.Pts.y2_kv;
                             
            % f = exp(-r^2)*sin(t); 
            f1 = exp(-r.^2).*sin(t);
            gradf1 = [ -2*r.*exp(-r.^2).*sin(t) ; r.^(-1).*exp(-r.^2).*cos(t) ];
            Lapf1  = (-2 + 4*r.^2).*exp(-r.^2).*sin(t) - 2*exp(-r.^2).*sin(t) ...
                         - r.^(-2).*exp(-r.^2).*sin(t);
            disp(['InfAnnulus: Max error of gradient of exp(-r^2)*sin(t): ',...
                                    num2str( max(Diff.grad*f1 - gradf1) ) ]); 
            disp(['InfAnnulus: Max error of Laplacian of exp(-r^2)*sin(t): ',...
                                    num2str( max(Diff.Lap*f1 - Lapf1) ) ]); 
            
        end        
        function Interp = ComputeInterpolationMatrix(this,interp1,interp2,fullInterpGrid,saveBool)

            if(nargin<2)
               interp1 = (-1+0.05:0.05:1-0.05)';
               interp2 = (0:0.02:1)';
               fullInterpGrid = true;
               saveBool = true;
            end
            
            Interp = ComputeInterpolationMatrix@Polar_SpectralFourierNoOrigin ...
                       (this,interp1,interp2,fullInterpGrid,saveBool);
            
%             %Check Accuracy
%        
%             r = this.Pts.y1_kv;
%             t = this.Pts.y2_kv;
%               
%             rPlot = Interp.pts1;
%             tPlot = Interp.pts2;
%             
%             % f = exp(-(r-Rmin)^2 cos(3t)
%             f1 = exp(-(r-this.RMin).^2).*cos(3*t);
%             f1Plot = exp(-(rPlot-this.RMin).^2).*cos(3*tPlot);
%             disp(['InfAnnulus: Max error of interpolation of exp(-(r-Rmin)^2 cos(3t): ',...
%                       num2str( max(abs(Interp.InterPol*f1 - f1Plot)) ) ])
%             
%             figure
%             doPlots_SC_Polar(Interp, this.Pts, f1)
                   
                   
        end          
    end
end