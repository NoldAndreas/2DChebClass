classdef HalfSpace < SpectralSpectral & ConvolutionPointwise

    properties                
        y2Min = 0;
        L1,L2
        y10 = 0;
        LConv=[];
        L2Conv=[];
    end
    
    methods        
        function this = HalfSpace(Geometry)
            this@SpectralSpectral(Geometry.N(1),Geometry.N(2));

            this.L1    = Geometry.L1;
            this.L2    = Geometry.L2;
            if(isfield(Geometry,'y10'))
                this.y10 = Geometry.y10;
            end
            if(isfield(Geometry,'y2Min'))
                this.y2Min = Geometry.y2Min;
            end
            if(isfield(Geometry,'Conv') && ~isempty(Geometry.Conv))
                this.ConvN = Geometry.Conv.N;
                if(isfield(Geometry.Conv,'L'))
                    this.LConv = Geometry.Conv.L;                
                end
                if(isfield(Geometry.Conv,'L2'))
                    this.L2Conv = Geometry.Conv.L2;
                end
            end
            this.polar = 'cart';            
            InitializationPts(this);            
        end                         
    end
    
    methods (Access = public)
        function M_conv = ComputeConvolutionMatrix(this,f,saveBool,ptsCheck)
            if(nargin>=4)
                M_conv  = ComputeConvolutionMatrix_Pointwise(this,f,ptsCheck);
            else
                M_conv  = ComputeConvolutionMatrix_Pointwise(this,f);
            end

            if((nargin == 3) && islogical(saveBool) && saveBool)
                this.Conv = M_conv;
            end            
        end 
        
        %*** Mapping function ***
        function [y1,dy1,dx,ddx,dddx,ddddx] = PhysSpace1(this,x1)                        
            [y1,dy1,dx,ddx,dddx,ddddx] = SqrtMap(x1,this.L1,inf);
            y1 = y1 + this.y10;
        end
        function xf = CompSpace1(this,y1)            
            xf  = InvSqrtMap(y1-this.y10,this.L1,inf);
        end
        function [th,dth,dx,ddx,dddx,ddddx] = PhysSpace2(this,x2)                
            [th,dth,dx,ddx,dddx,ddddx] = QuotientMap(x2,this.L2,this.y2Min,inf);
        end
        function x2 = CompSpace2(this,y2)                        
            x2  = InvQuotientMap(y2,this.L2,this.y2Min,inf);
        end        
        %*** Mapping functions for the convolution ***
        function [y1D,dy1D] = PhysSpace1_Sub(this,y10,y20,x,xCheb,f)
            %return y1D = y1-y10
            %[z,dz] = QuadMapAB(x,LConv,-inf,inf);    
            if((nargin > 4)&& isempty(this.LConv))
                Lmax        = 10;
                this.LConv = FindOptimalL(@SqrtMap,xCheb,f,Lmax);            
            elseif((nargin == 4) && isempty(this.LConv))
                error('HalfSpace:PhysSpace1_Sub: L1Conv not initialized.');
            end
            [y1D,dy1D]  = SqrtMapAB(x,this.LConv,-inf,inf);                    
        end        
        function x = CompSpace1_Sub(this,y10,y20,y1D)
             %x = InvQuadMapAB(z-y10,LConv,-inf,inf);        
             x = InvSqrtMapAB(y1D,this.LConv,-inf,inf);        
        end        
        
        function [y2D,dy2D,d_Pade,ep_Pade] = PhysSpace2_Sub(this,y10,y20,x,xCheb,f,opts)
            if(nargin < 7)
                LConv  = this.LConv;
                L2Conv = this.L2Conv;
                y2Min  = this.y2Min;
            else
                L2Conv = opts.L2Conv;
                y2Min  = opts.y2Min;
            end
            
            %returns the difference y2D = y2 - y20 and its derivative
            d_Pade  = []; ep_Pade = [];
            if(y20==inf)
                Lmax       = 10;
                LConv      = FindOptimalL(@SqrtMap,xCheb,f,Lmax);
                [y2D,dy2D] = SqrtMap(x,LConv,inf);
            elseif(y20>L2Conv)            
                ep       = sqrt(1+(y20/L2Conv)^2) - (y20/L2Conv);                    
                d        = 0;                    
                if(nargin >= 6)
                    %1st Step: Analyze using intuitive guess
                    [x1,dx1] = Tref(xCheb,d,ep); %x1 = x; dx1 = ones(size(x1));                
                    [y2,dy2] = QuotientMap(x1,y20-y2Min,y2Min,inf);
                    y2D      = y2 - y20;
                    fx       = f(y2D);
                    
                    %2nd Step: Do Pade-Approximation to obtain ep,d-values                                        
                    [d_Pade,ep_Pade] = GetPolesPadeApproximation(fx,2);
                                        
                    %Compute values using new intermediate map
                    [x,dx2Tee] = Tref(x,d_Pade,ep_Pade);
                else                    
                    dx2Tee  = ones(size(x));
                end        
                [x1,dx1] = Tref(x,d,ep); %x1 = x; dx1 = ones(size(x1));                
                [y2,dy2] = QuotientMap(x1,y20-y2Min,y2Min,inf);
                y2D      = y2 - y20;
                dy2D     = dy2.*dx1.*dx2Tee;
            else
                [y2,dy2D] = QuotientMap(x,L2Conv,y2Min,inf);            
                y2D       = y2 - y20;
            end
        end
        function x = CompSpace2_Sub(this,y10,y20,y2D,d_Tee,ep_Tee)
             if(y20<=this.L2Conv)
                    x     = InvQuotientMap(y2D + y20,this.L2Conv,this.y2Min,inf);
             elseif(y20 == inf)                                  
                 x     = InvSqrtMap(y2D,this.LConv,inf);
             else   %if(y20>this.L2Conv)
                ep     = sqrt(1+(y20/this.L2Conv)^2) - (y20/this.L2Conv);
                %ep     = 0.1;
                x1     = InvQuotientMap(y2D + y20,y20-this.y2Min,this.y2Min,inf);
                x      = InvTref(x1,0,ep);            
                if((nargin == 6) && (~isempty(d_Tee)) && (~isempty(ep_Tee)))
                    %Compute values using new intermediate map
                    x = InvTref(x,d_Tee,ep_Tee);                
                end             
             end
        end
       
        function do1DPlotNormal(this,V,sym,col)
            if(nargin < 3)
                sym = 'o';
            end 
            if(nargin < 4)
                colEdge = 'k';
                colFace = 'g';                
            else
                colEdge = col;
                colFace = col;
            end
            
            V = V(:);
            
            %V: Vector of length N2
            ptsM        = GetInvCartPts(this,inf,15);
            y2Max       = ptsM.y2_kv;%10
            mark        = (this.Pts.y1_kv == inf);            
            IP          = ComputeInterpolationMatrix(this,1,(-1:0.01:CompSpace2(this,y2Max))',true);         
            
            if(length(V) == this.M)
                V = V(mark);
            end
            IP.InterPol = IP.InterPol(:,mark);
            
            PtsCart     = GetCartPts(this);
            IPPtsCart   = GetCartPts(this,IP.pts1,IP.pts2);          
            y2MinC      = min(PtsCart.y2_kv(mark));
            if(~isempty(sym))
                plot(PtsCart.y2_kv(mark)-y2MinC,V,sym,'MarkerEdgeColor',colEdge,'MarkerFaceColor',colFace); hold on;
            end
            plot(IPPtsCart.y2_kv-y2MinC,IP.InterPol*V,colEdge,'linewidth',1.5);
            xlim([0 max(IPPtsCart.y2_kv)]);
            xlabel('$y_{2,Cart}$','Interpreter','Latex','fontsize',25);            
            set(gca,'fontsize',20);                        
            set(gca,'linewidth',1.5);       
        end        
        function [pts_y1,V_pts] = do1DPlotParallel(this,V)  
            %V: Vector of length N2
            y1W         = 5;
            mark        = (this.Pts.y2_kv == inf);
            x1W         = CompSpace1(this,y1W);
            IP          = ComputeInterpolationMatrix(this,(-x1W:0.01:x1W)',1,true);
            IP.InterPol = IP.InterPol(:,mark);            
            
            plot(this.Pts.y1_kv(mark),V,'o','MarkerEdgeColor','k','MarkerFaceColor','g'); 
            hold on;
            pts_y1 = IP.pts1;
            V_pts  = IP.InterPol*V;
            plot(pts_y1,V_pts,'linewidth',1.5);
            xlim([min(IP.pts1) max(IP.pts1)]);
            xlabel('$y_{1}$','Interpreter','Latex','fontsize',25);        
            set(gca,'fontsize',20);                        
            set(gca,'linewidth',1.5);       
        end     
        
    end
end