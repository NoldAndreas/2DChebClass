classdef InfCapillarySkewed < InfCapillary 

    properties  
        alpha = pi/2;       
    end
    
    methods        
        function this = InfCapillarySkewed(Geometry)
            Geometry.y2Min = Geometry.y2Min/sin(Geometry.alpha);
            Geometry.y2Max = Geometry.y2Max/sin(Geometry.alpha);
            
            this@InfCapillary(Geometry);
                        
            this.alpha  = Geometry.alpha;
            
            this.polar  = 'cartSkewed';
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
            [ptsCart.y1_kv,ptsCart.y2_kv] = SkewedGrid(pts_y1,pts_y2,this.alpha);
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
            
            Diff.grad = sparse([Diff.Dy1;Diff.Dy2]);            
            
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
            % 
            %Diff.DDy1 = ;
            
        end        
        function M_conv = ComputeConvolutionMatrix(this,f,saveBool,ptsCheck)
            if(nargin>=4)
                M_conv  = ComputeConvolutionMatrix@InfCapillary(this,f,ptsCheck);
            else
                M_conv  = ComputeConvolutionMatrix@InfCapillary(this,f);
            end            
            
            M_conv = M_conv*sin(this.alpha);
            
            if((nargin == 3) && islogical(saveBool) && saveBool)
                this.Conv = M_conv;
            end            
            
        end
    end
end
        