classdef Annulus < Polar_SpectralFourierNoOrigin

    properties        
        RMin,RMax        
    end
    
    methods        
        function this = Annulus(Geometry)
            this@Polar_SpectralFourierNoOrigin(Geometry.N(1),Geometry.N(2));
            
            if(isfield(Geometry,'R_in'))
                this.RMin   = Geometry.R_in; 
                this.RMax   = Geometry.R_out; 
            else
                this.RMin   = Geometry.RMin; 
                this.RMax   = Geometry.RMax; 
            end
            
             if(isfield(Geometry,'Origin'))
                 this.Origin = Geometry.Origin;
             end
            
            InitializationPts(this);                                    
        end                         
    end
    
    methods (Access = public)        
        
        function [r,dr,dx,ddx,dddx,ddddx] = PhysSpace1(this,x)
            [r,dr,dx,ddx,dddx,ddddx] = LinearMap(x,this.RMin,this.RMax);
        end
        function xf = CompSpace1(this,r)
            xf = InvLinearMap(r,this.RMin,this.RMax);
        end
        function [th,dth,dx,ddx,dddx,ddddx] = PhysSpace2(this,xT)    
            [th,dth,dx,ddx,dddx,ddddx] = LinearMap01(xT,0,2*pi);
        end
        function xf = CompSpace2(this,th)
            xf = th/(2*pi);
        end                      
        function [int,area] = ComputeIntegrationVector(this)
            int = ComputeIntegrationVector@Polar_SpectralFourierNoOrigin(this);
            area = pi*(this.RMax^2-this.RMin^2);
            
            %Check Accuracy            
            if(nargout == 1)                
                disp(['Annulus: Error of integration of area (ratio): ',...
                                        num2str(1-sum(this.Int)/area)]);      

                % f = r; 
                intF1 = 2*pi/3*(this.RMax^3-this.RMin^3);
                f1 = this.Pts.y1_kv;
                disp(['Annulus: Error of integration of r (ratio): ',...
                                        num2str(1-sum(this.Int*f1)/intF1)]); 

                % f = r*sin(theta); 
                intF2 = 0;
                f2 = this.Pts.y1_kv.*sin(this.Pts.y2_kv);
                disp(['Annulus: Error of integration of r*sin(theta) (ratio): ',...
                                        num2str(intF2-sum(this.Int*f2))]); 

                % f = r*sin^2(theta); 
                intF3 = pi/3*(this.RMax^3-this.RMin^3);
                f3 = this.Pts.y1_kv.*sin(this.Pts.y2_kv).^2;
                disp(['Annulus: Error of integration of r*sin^2(theta) (ratio): ',...
                                        num2str(1-sum(this.Int*f3)/intF3)]);

                % f = exp(-r^2); 
                intF4 = -pi*(exp(-this.RMax^2)-exp(-this.RMin^2));
                f4 = exp(-this.Pts.y1_kv.^2);
                disp(['Annulus: Error of integration of exp(-r^2) (ratio): ',...
                                        num2str(1-sum(this.Int*f4)/intF4)]); 
            end
            
        end        
        function Diff = ComputeDifferentiationMatrix(this)
            Diff = ComputeDifferentiationMatrix@Polar_SpectralFourierNoOrigin(this);
            
            %Check Accuracy

            % grad f = rhat \partial_r f + that 1/r \partial_t f
            % div V  = 1/r \partial_r (r V_r) + 1\r \partial_t V_t
            % Lap    = \partial_r^2 + 1/r \partial_r + 1/r^2 \partial_t^2 
            
            size1 = size(this.Pts.y1_kv);
            size2 = size(this.Pts.y2_kv);
            r = this.Pts.y1_kv;
            t = this.Pts.y2_kv;
            
            % f = r;        
            f1     = r;
            gradf1 = [ ones(size1) ; zeros(size2) ];
            Lapf1  = r.^(-1);
            disp(['Annulus: Max error of gradient of r: ',...
                                    num2str( max(Diff.grad*f1 - gradf1) ) ]); 
            disp(['Annulus: Max error of Laplacian of r: ',...
                                    num2str( max(Diff.Lap*f1 - Lapf1) ) ]); 
                                                
            % f = r*sin(theta); 
            f2     = r.*sin(t);
            gradf2 = [ sin(t) ; cos(t) ];
            Lapf2  = r.^(-1).*sin(t) - r.^(-1).*sin(t);
            disp(['Annulus: Max error of gradient of r*sin(t): ',...
                                    num2str( max(Diff.grad*f2 - gradf2) ) ]); 
            disp(['Annulus: Max error of Laplacian of r*sin(t): ',...
                                    num2str( max(Diff.Lap*f2 - Lapf2) ) ]); 
                               
            % f = exp(-r^2)*sin(t); 
            f3 = exp(-r.^2).*sin(t);
            gradf3 = [ -2*r.*exp(-r.^2).*sin(t) ; r.^(-1).*exp(-r.^2).*cos(t) ];
            Lapf3  = (-2 + 4*r.^2).*exp(-r.^2).*sin(t) - 2*exp(-r.^2).*sin(t) ...
                         - r.^(-2).*exp(-r.^2).*sin(t);
            disp(['Annulus: Max error of gradient of exp(-r^2)*sin(t): ',...
                                    num2str( max(Diff.grad*f3 - gradf3) ) ]); 
            disp(['Annulus: Max error of Laplacian of exp(-r^2)*sin(t): ',...
                                    num2str( max(Diff.Lap*f3 - Lapf3) ) ]); 
            
        end
        function Interp = ComputeInterpolationMatrix(this,interp1,interp2,fullInterpGrid,saveBool)

            if(nargin<2)
               interp1 = (-1:0.05:1)';
               interp2 = (0:0.02:1)';
               fullInterpGrid = true;
               saveBool = true;
            end
            
            Interp = ComputeInterpolationMatrix@Polar_SpectralFourierNoOrigin ...
                       (this,interp1,interp2,fullInterpGrid,saveBool);
            
            %Check Accuracy
       
%             r = this.Pts.y1_kv;
%             t = this.Pts.y2_kv;
%               
%             rPlot = Interp.pts1;
%             tPlot = Interp.pts2;
            
%             % f = r
%             f1 = r;
%             f1Plot = rPlot;
%             disp(['Annulus: Max error of interpolation of r: ',...
%                       num2str( max(abs(Interp.InterPol*f1 - f1Plot)) ) ]); 
%             
%             figure
%             doPlots_SC_Polar(Interp, this.Pts, f1)
%             
%             % f = r sin^2 t
%             f2 = r.*sin(t).^2;
%             f2Plot = rPlot.*sin(tPlot).^2;
%             disp(['Annulus: Max error of interpolation of r sin(t)^2: ',...
%                       num2str( max(abs(Interp.InterPol*f2 - f2Plot)) ) ]); 
%             
%             figure
%             doPlots_SC_Polar(Interp, this.Pts, f2)
%             
%             
%             % f = exp(-(r-Rmin)^2 cos^3(t)
%             f3 = exp(-(r-this.RMin).^2).*cos(3*t);
%             f3Plot = exp(-(rPlot-this.RMin).^2).*cos(3*tPlot);
%             disp(['Annulus: Max error of interpolation of exp(-(r-Rmin)^2 cos^3(t): ',...
%                       num2str( max(abs(Interp.InterPol*f3 - f3Plot)) ) ])
%             
%             figure
%             doPlots_SC_Polar(Interp, this.Pts, f3)
                   
                   
        end  
        function PlotBorders(this)
            ComputeIndices(this);
             
            rIn = this.Pts.y1_kv(this.Ind.inR);
            geomIn.R = rIn(1);
            geomIn.N = this.Pts.N2;
            
            cIn = Circle(geomIn);
                       
            rOut = this.Pts.y1_kv(this.Ind.outR);
            geomOut.R = rOut(1);
            geomOut.N = this.Pts.N2;
            
            cOut = Circle(geomOut);
            
            cIn.PlotFunction;
            cOut.PlotFunction;
           
        end 
    end
end