classdef InfCapillary_Interface < M1SpectralSpectral

    properties
        L1       
        y2Min,y2Max        
        dy1
    end
    
    methods        
        function this = InfCapillary_Interface(Geometry)

            this@M1SpectralSpectral(Geometry.N(1),Geometry.N(2));

            this.L1      = Geometry.L1; 
            this.y2Min   = Geometry.y2Min; 
            this.y2Max   = Geometry.y2Max; 
            this.polar   = 'cart';
            this.dy1     = zeros(Geometry.N(2),1);                        

            InitializationPts(this);            
            
            %this.dy1 = 2*(this.Pts.x2*(this.y2Min - this.y2Max)/2).^2;
            InitializationPts(this);            
        end
    end
    
    methods (Access = public)

        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)
            %[y1_kv,y2_kv,J,DDy1,DDy2] = Comp_to_Phys(x1,x2)
            %[y1_kv,dy1dx1] = LinearMap(x1,-L1,L1); %,dx1,ddx1,dddx,ddddx]                
            n  = length(x1);           

            [y2_kv,dy2dx2] =  LinearMap(x2,this.y2Min,this.y2Max);
            ddy2ddx2 = zeros(size(y2_kv));
            [h2,dh2,ddh2]  =  h(this,x2);
            % c2 = 0.3; [y1_kv,Diffy1] =  M1SqrtMap(x1,c2,inf);% M1Tref(x1,0,0.07);                         
            [y1_kv,Diffy1] =  M1QuadMap(x1,this.L1,inf);% M1Tref(x1,0,0.07);                                     
            y1_kv          =  y1_kv + h2;

            if(nargout >= 3)
                J        = zeros(n,2,2);
                J(:,1,1) = Diffy1.dydx; 
                J(:,1,2) = dh2.*dy2dx2;
                J(:,2,1) = zeros(n,1);
                J(:,2,2) = dy2dx2;            
            end

            if(nargout >= 4)
                dH1        = zeros(n,2,2);
                dH1(:,1,1) = Diffy1.dyddx; %d2y1dx1dx1; 
                dH1(:,1,2) = zeros(this.N1*this.N2,1);%d2y1dx1dx2;             
                dH1(:,2,1) = zeros(this.N1*this.N2,1);%d2y1dx1dx2;     
                dH1(:,2,2) = ddh2.*(dy2dx2).^2+dy2dx2.*ddy2ddx2;
            end

            if(nargout >= 4)
                dH2        = zeros(n,2,2);
                dH2(:,1,1) = zeros(n,1); %d2y2dx1dx1; 
                dH2(:,1,2) = zeros(n,1); %d2y2dx1dx2;             
                dH2(:,2,1) = zeros(n,1); %d2y2dx1dx2;     
                dH2(:,2,2) = zeros(n,1); %d2y2dx2dx2;             
            end

            %[z,dz,dx,ddx,dddx,ddddx] = QuadMapAB(xR,L1,-L1,L1);
        end
        function [x1,x2] = CompSpace(this,y1,y2)
            x2  = InvLinearMap(y2,this.y2Min,this.y2Max);
            h2  =  h(this,x2);
            x1  = InvQuadMap(y1-h2,this.L1,inf);
            
             %xf  = InvQuadMapAB(z,L1,-L1,L1);
            %error = MException('InfCapillary_Interface',...
%                                'inverse function not yet implemented');
%             throw(error);
        end     
        function [h2,dh2,ddh2] = h(this,x2)%,y)

            X2Diff = barychebdiff(this.Pts.x2);
            IP     = barychebevalMatrix(this.Pts.x2,x2);
            L      = (this.y2Max - this.y2Min)/2;
            
            h2     = IP*this.dy1;            
            dh2    = IP*X2Diff.Dx*this.dy1/L; 
            ddh2   = IP*X2Diff.DDx*this.dy1/(L^2); 
%
            %h2       = kron(h2,ones(this.N2,1));
            %dh2      = kron(dh2,ones(this.N2,1));
            %ddh2     = kron(ddh2,ones(this.N2,1));
            
%             y     = y - (this.y2Min + this.y2Max)/2;
%             c     = 2;
%             h2    = c*y.^2;
%             dh2   = 2*c*y;
%             ddh2  = 2*c*ones(size(y));
% %
            %h2   = c*y;
            %dh2  = c*ones(size(y));
            %ddh2 = zeros(size(y));
        end   
        
    end
end