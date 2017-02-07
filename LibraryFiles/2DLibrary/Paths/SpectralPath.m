classdef (Abstract) SpectralPath < handle
   
    properties
        Pts,Diff,Interp,Conv             
        N        
        normal,tang
        IntSc,IntNormal,IntTang
        polar = 'undefined';    
        ds
    end
    
    methods (Abstract)
        [y1,y2] = f_path(this,t)
    end
     
    methods
        function this = SpectralPath(N,polar)
            this.Pts.t   = ClenCurtFlip(N-1);                        
            this.Pts.N   = N;
            this.N       = N;            
            this.polar   = polar;            
        end
        
        function InitializationPts(this)
            [this.Pts.y1_kv , this.Pts.y2_kv,dy1_dt,dy2_dt]  = f_path(this,this.Pts.t);    
            dy1_dt(dy1_dt == inf | dy1_dt == -inf) = 0;
                        
            if(strcmp(this.polar,'polar'))
                r_dy2_dt = (this.Pts.y1_kv).*dy2_dt;
            else
                r_dy2_dt = dy2_dt;
            end
            this.ds   = sqrt(dy1_dt.^2 + r_dy2_dt.^2);            
            dsInv     = 1./this.ds;
            dsInv(this.ds == 0) = 0;
            
            this.normal  =  diag(dsInv)*[diag(r_dy2_dt) diag(-dy1_dt)];                
            this.tang    =  diag(dsInv)*[diag(-dy1_dt)  diag(r_dy2_dt)];    
            
            
        end        
        function int = ComputeIntegrationVector(this)            
            [h1,wInt]    = ClenCurtFlip(this.N-1);
            this.IntSc   = this.ds'.*wInt;
            this.IntSc(this.IntSc == inf)   = 0; %assuming that the function to integrate
                                         %converges to zero fast enough, s.t. the
                                         %value at inf can be ignored                                                                                
            this.IntSc(isnan(this.IntSc)) = 0;                                         
            
            this.IntNormal =  this.IntSc*this.normal;
            this.IntTang   =  this.IntSc*this.tang;             

            int = this.IntSc;
        end        
        function Diff = ComputeDifferentiationMatrix(this)
            %Differentiation  
            Diff_t   = barychebdiff(t,4);
            this.Diff     = struct('Dy1',diag(dy1_dt)*Diff_t.Dx,...
                              'Dy2',diag(dy2_dt)*Diff_t.Dx);    
            Diff = this.Diff;
        end        
        function interp = ComputeInterpolationMatrix(this,t_interpPts)
            %if(exist('interpPts','var'))        
            [y1_IP,y2_IP]    = f_path(this,t_interpPts);
            IP               = barychebevalMatrix(this.Pts.t,t_interpPts);         
            interp           = struct('InterPol',IP,...
                                      'pts_y1',y1_IP,'pts_y2',y2_IP,...                        
                                      'Nplot',length(t_interpPts),...
                                      'N',this.N,'t',t_interpPts);
            %end
            this.Interp = interp;
        end
        function Conv = ComputeConvolutionMatrix(this,f_conv)
            this.Conv = this.Int.*f_conv(t);
            Conv      = this.Conv;
        end
        
        function PlotPath(this,polar)
            if((nargin > 1) && strcmp(polar,'polar'))
                pts = Pol2CartPts(this.Pts);                
            else
                pts = this.Pts;
            end
            plot(pts.y1_kv,pts.y2_kv,'m','linewidth',1.5);            
        end
        
        
        function ptsCart = GetCartPts(this,pts_y1,pts_y2)
            
            if(nargin == 1)
                pts_y1 = this.Pts.y1_kv;
                pts_y2 = this.Pts.y2_kv;
            end
            
            if(strcmp(this.polar,'polar'))
                [ptsCart.y1_kv,ptsCart.y2_kv] = pol2cart(pts_y2,pts_y1);
            elseif(strcmp(this.polar,'cart'))                
                ptsCart.y1_kv = pts_y1;  ptsCart.y2_kv = pts_y2;
            else
                exc = MException('Shape:GetCartPts','select {polar,cart}');
                throw(exc);                
            end
            
        end 
                
    end
        
    
end

