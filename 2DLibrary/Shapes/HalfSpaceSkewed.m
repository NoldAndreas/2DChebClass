classdef HalfSpaceSkewed < HalfSpace 

    properties  
        alpha = pi/2;       
    end
    
    methods        
        function this = HalfSpaceSkewed(Geometry)
            Geometry.y2Min = Geometry.y2Min/sin(Geometry.alpha);
            this@HalfSpace(Geometry);
                        
            if(isfield(Geometry,'alpha'))                
                this.alpha  = Geometry.alpha;
            elseif(isfield(Geometry,'alpha_deg'))
               this.alpha  = Geometry.alpha_deg*pi/180; 
            end
            
            this.polar = 'cartSkewed';
            InitializationPts(this);      
        end                         
    end
    
    methods (Access = public)
        
        %*** Mapping function ***
        function ptsCart = GetCartPts(this,pts_y1,pts_y2)            
            if(nargin == 1)
                pts_y1 = this.Pts.y1_kv;
                pts_y2 = this.Pts.y2_kv;
            end
                                  
            if(strcmp(this.polar,'cartSkewed'))
                [ptsCart.y1_kv,ptsCart.y2_kv] = SkewedGrid(pts_y1,pts_y2,this.alpha);
            else
                exc = MException('HalfSpaceSkewed:GetCartPts','select cartSkewed');
                throw(exc);                
            end            
        end       
        function pts = GetInvCartPts(this,ptsCart_y1,ptsCart_y2)            
                                  
            if(strcmp(this.polar,'cartSkewed'))
                [pts.y1_kv,pts.y2_kv] = InvSkewedGrid(ptsCart_y1,ptsCart_y2,this.alpha);
            else
                exc = MException('HalfSpaceSkewed:GetInvCartPts','select cartSkewed');
                throw(exc);                
            end            
        end       
        function [Int,Int1,Int2]    = ComputeIntegrationVector(this)
            %This gives the integration in the real (cartesian) space!
            [Int,Int1,Int2] = ComputeIntegrationVector@SpectralSpectral(this);  
            Int             = Int*sin(this.alpha);
            Int2            = Int2*sin(this.alpha);

            this.Int  = Int;
        end
        function Diff = ComputeDifferentiationMatrix(this)
            DiffSk = ComputeDifferentiationMatrix@SpectralSpectral(this);
            
            ATInv = zeros(2,2);    
            ATInv(1,1) = 1;
            ATInv(2,1) = -1/tan(this.alpha);
            ATInv(2,2) = 1/sin(this.alpha);
            %A(1,2) = cos(this.alpha);
            %A(2,2) = sin(this.alpha);
            
            % y = A*yt
            %nabla_y f = (A^(T))^(-1)* nabla_y'*f
            
            h         = zeros(this.M,this.M,2);
            h(:,:,1)  = DiffSk.Dy1;
            h(:,:,2)  = DiffSk.Dy2;
            
            hh        = MultScalarMatrVecOp(ATInv,h);
            Diff.Dy1  = hh(:,:,1);
            Diff.Dy2  = hh(:,:,2);
            
            %Diff.Dy1  = ATInv(1,1)*DiffSk.Dy1 + ATInv(1,2)*DiffSk.Dy2;
            %Diff.Dy2  = ATInv(2,1)*DiffSk.Dy1 + ATInv(2,2)*DiffSk.Dy2;
            Diff.div   = [Diff.Dy1 Diff.Dy2];
            Diff.grad  = sparse([Diff.Dy1;Diff.Dy2]);
            
            H          = zeros(this.M,this.M,2,2);
            H(:,:,1,1) = DiffSk.DDy1;
            H(:,:,1,2) = DiffSk.Dy1Dy2;
            H(:,:,2,1) = DiffSk.Dy1Dy2;
            H(:,:,2,2) = DiffSk.DDy2;
            
            w          = MultScalarMatrOp(ATInv,H);   
                        
            ht         = w(:,:,1,2);
            w(:,:,1,2) = w(:,:,2,1);
            w(:,:,2,1) = ht;
            
       %     w          = permute(w,[1 2 4 3]);
            w          = MultScalarMatrOp(ATInv,w);            
            
            ht         = w(:,:,1,2);
            w(:,:,1,2) = w(:,:,2,1);
            w(:,:,2,1) = ht;            
            
            Diff.DDy1    = w(:,:,1,1);
            Diff.Dy1Dy2  = w(:,:,1,2);
            Diff.DDy2    = w(:,:,2,2);%DiffSk.DDy2;%
            
            Diff.Lap     = Diff.DDy1 + Diff.DDy2;                        
            %Diff.DDy1 = ;
            
        end        
        function M_conv = ComputeConvolutionMatrix(this,f,saveBool)
            
            %M_conv  = ComputeConvolutionMatrix_Pointwise(this,f)*sin(this.alpha);
            M_conv  = ComputeConvolutionMatrix_PointwiseX(this,f);            
            
            if((nargin == 3) && islogical(saveBool) && saveBool)
                this.Conv = M_conv;
            end
        end
        
         function do1DPlotParallel(this,V)    
             global PersonalUserOutput
            if(~PersonalUserOutput)
                return;
            end            
            %V: Vector of length N2
            y1N         = 5;
            y1W         = y1N/sin(this.alpha);
            mark        = (this.Pts.y2_kv == inf);
            x1W         = CompSpace1(this,y1W);
            IP          = ComputeInterpolationMatrix(this,(-x1W:0.01:x1W)',1,true);
            IP.InterPol = IP.InterPol(:,mark);            
            
            plot(sin(this.alpha)*this.Pts.y1_kv(mark),V,'o','MarkerEdgeColor','k','MarkerFaceColor','g'); 
            hold on;
            plot(sin(this.alpha)*IP.pts1,IP.InterPol*V,'linewidth',1.5);
            xlim([-y1N y1N]);
            xlabel('$\sin(\alpha) y_{1}$','Interpreter','Latex','fontsize',25);        
            set(gca,'fontsize',20);                        
            set(gca,'linewidth',1.5);       
        end     
    end
end
        