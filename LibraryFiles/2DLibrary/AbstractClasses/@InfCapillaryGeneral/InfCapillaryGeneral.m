classdef (Abstract) InfCapillaryGeneral < SpectralSpectral & ConvolutionPointwise

    properties        
        L1       
        y2Min,y2Max        
        L1Conv,L2Conv
        y10 = 0
    end
    
    methods        
        function this = InfCapillaryGeneral(Geometry)

            this@SpectralSpectral(Geometry.N(1),Geometry.N(2));

            this.L1      = Geometry.L1; 
            this.y2Min   = Geometry.y2Min; 
            this.y2Max   = Geometry.y2Max; 
            if(isfield(Geometry,'y10'))
                this.y10 = Geometry.y10;
            end

            InitializationPts(this);            
            
            if(isfield(Geometry,'Conv') && ~isempty(Geometry.Conv))
                this.L1Conv = Geometry.Conv.L1;
                this.L2Conv = Geometry.Conv.L2;
                this.ConvN  = Geometry.Conv.N;
            end
            this.polar = 'cart';
        end
    end
    
    methods (Access = public)
%         function [y1,dy1,dx,ddx,dddx,ddddx] = PhysSpace1(this,x1)
%             [y1,dy1,dx,ddx,dddx,ddddx] = SqrtMap(x1,this.L1,inf);            
%             y1 =  this.y10 + y1; 
%         end
%         function xf = CompSpace1(this,y1)
%             xf  = InvSqrtMap(y1 - this.y10,this.L1,inf);
%         end
        function [th,dth,dx,ddx,dddx,ddddx] = PhysSpace2(this,x2)                
            [th,dth,dx,ddx,dddx,ddddx] = LinearMap(x2,this.y2Min,this.y2Max);
        end
        function x2 = CompSpace2(this,y2)                        
            x2  = InvLinearMap(y2,this.y2Min,this.y2Max);
        end              
        
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
                
        function [dataDisk] = AverageDiskPt(this,y20,r,N,sphere)
            shape.N  = N;
            shape.R  = r;     
            
            if((nargin == 5) && strcmp(sphere,'sphere'))
                shape.sphere = true;
            end
            
            if(y20 >= this.y2Max + r)
                dataDisk.pts        = 0;
                dataDisk.ptsPolLoc  = 0;
                dataDisk.int        = 0;
                dataDisk.area       = -1;                                
                return;
            end
                        
            if((y20 > this.y2Max) && (y20 < this.y2Max + r)) 
                %Segment
                shape.h            = y20 - this.y2Max;
                shape.S            = -1;
                                
                area               = Segment(shape);                
                dataDisk.pts       = area.Pts;                
                
            elseif(y20 == this.y2Max)                
                %Semicircle
                shape.th1       = pi;
                shape.th2       = 2*pi;        
                
                area            = Wedge(shape);                
                dataDisk.pts    = Pol2CartPts(area.Pts);                            
            else
                exc = MException('InfCapillary:AverageDiskPt','case not implemented');
                throw(exc);
            end                                    
                        
            %Shift in y2-direction
            dataDisk.ptsPolLoc   = Cart2PolPts(area.Pts);            
            dataDisk.pts.y2_kv   = dataDisk.pts.y2_kv + y20;
            
            [dataDisk.int,dataDisk.area]     = area.ComputeIntegrationVector();
                        
            scatter(dataDisk.pts.y1_kv,dataDisk.pts.y2_kv,'.');            
        end                
        function [dataCircle] = AverageCirclePt(this,y20,r,Ncircle)
            
            shapeLine.N = Ncircle;
            shapeLine.R = r;
            
            if(y20 >= this.y2Max + r)                
                dataCircle.pts        = 0;
                dataCircle.ptsPolLoc  = 0;
                dataCircle.int        = 0;
                dataCircle.area       = -1;
                return;
            end                          
            
            if((y20 >= this.y2Max) && (y20 < this.y2Max + r)) 
                th                 = acos((y20 - this.y2Max)/r);
                shapeLine.th1      = 3/2*pi - th;
                shapeLine.th2      = 3/2*pi + th;
                
                line               = Arc(shapeLine);
                dataCircle.pts     = Pol2CartPts(line.Pts);
            else
                exc = MException('InfCapillary:AverageDiskPt','case not implemented');
                throw(exc);
            end
                        
            dataCircle.pts       = Pol2CartPts(line.Pts);            
            dataCircle.pts.y2_kv = dataCircle.pts.y2_kv + y20;
            dataCircle.ptsPolLoc = line.Pts;            
                        
            [dataCircle.int,dataCircle.area] = line.ComputeIntegrationVector();                                                   
            
            %Plot data            
            scatter(dataCircle.pts.y1_kv,dataCircle.pts.y2_kv,'.');            
        end        
        function dataBall = AverageBallPt(this,y20,r,N)
            shape.N  = N;
            shape.R  = r;     
            
            if(y20 >= this.y2Max + r)                
                dataBall.pts        = 0;
                dataBall.ptsPolLoc  = 0;
                dataBall.int        = 0;
                dataBall.area       = -1;
                return;
            end    
                    
            if((y20 >= this.y2Max) && (y20 < this.y2Max + r)) 
                %1b. if part of disk is in HalfSpace  (>= half)
                %1b1. Integrate over segment in HalfSpace
                %shape.Origin    = [0,y20];
                th                = acos((y20 - this.y2Max)/r);
                shape.theta1      = pi-th;
                shape.theta2      = pi;                
                                
            %1c. if part of disk is in HalfSpace  (< half)    
            else
                exc = MException('HalfSpace_FMT:AverageDisk','case not implemented');
                throw(exc);                
            end            
            
            area               = Ball(shape); 
            dataBall.pts       = area.PtsCart;
            
            %Shift in y2-direction
            dataBall.ptsPolLoc = Cart2PolPts(area.Pts);            
            dataBall.pts.y2_kv = dataBall.pts.y2_kv + y20;
                        
            [dataBall.int,dataBall.area]     = area.ComputeIntegrationVector();
                                   
            scatter(dataBall.pts.y1_kv,dataBall.pts.y2_kv,'.');            
        end                        
                
        function [z,dz] = PhysSpace1_Sub(this,y10,y20,x,xCheb,fy1) 
            %[z,dz] = QuadMapAB(x,L1Conv,-inf,inf);        
            [z,dz] = SqrtMapAB(x,this.L1Conv,-inf,inf);        
            %z      = y10 + z;
        end
        function [y2d,dy2d,d_Pade,ep_Pade] = PhysSpace2_Sub(this,y10,y20,x,xCheb,fy2)
            [z,dy2d] = TrefInterval(x,this.y2Min,this.y2Max,y20,this.L2Conv);
            y2d      = z - y20;
            d_Pade   = 0;
            ep_Pade  = 0;
        end        
        function x = CompSpace1_Sub(this,y10,y20,y1D)
             %x = InvQuadMapAB(z-y10,L1Conv,-inf,inf);        
             x = InvSqrtMapAB(y1D,this.L1Conv,-inf,inf);        
         end
        function x = CompSpace2_Sub(this,y10,y20,y2D,d_Pade,ep_Pade)
             x      = InvTrefInterval(y2D + y20,this.y2Min,this.y2Max,y20,this.L2Conv);
         end
        
        [X,checkSum] = InterpolateAndIntegratePtsOrigin(this,ptsOr,data,weights);       
        
        function do1DPlotNormal(this,V)
            %plotLine(this,[inf inf],[this.y2Min this.y2Max],rho);            
            V = V(:);
            
            %V: Vector of length N2
            ptsM        = GetInvCartPts(this,inf,this.y2Max);
            y2Max       = ptsM.y2_kv;%10
            mark        = (this.Pts.y1_kv == inf);            
            IP          = ComputeInterpolationMatrix(this,1,(-1:0.01:CompSpace2(this,y2Max))',true);         
            
            if(length(V) == this.M)
                V = V(mark);
            end
            IP.InterPol = IP.InterPol(:,mark);
            
            PtsCart     = GetCartPts(this);
            IPPtsCart   = GetCartPts(this,IP.pts1,IP.pts2);            
            plot(PtsCart.y2_kv(mark),V,'o','MarkerEdgeColor','k','MarkerFaceColor','g'); 
            hold on;
            plot(IPPtsCart.y2_kv,IP.InterPol*V,'linewidth',1.5);
            xlim([min(IPPtsCart.y2_kv) max(IPPtsCart.y2_kv)]);
            xlabel('$y_{2,Cart}$','Interpreter','Latex','fontsize',25);            
            set(gca,'fontsize',20);                        
            set(gca,'linewidth',1.5);       
        end
       
             
        
    end
end