classdef Box_Interface < M1SpectralSpectral

    properties
        y1Min = 0
        y2Min = 0
        y1Max,y2Max            
    end
    
    methods        
        function this = Box_Interface(Geometry)
            
            this@M1SpectralSpectral(Geometry.N(1),Geometry.N(2));
            
            if(isfield(Geometry,'L1'))
                if(isfield(Geometry,'Origin'))
                    this.y1Min = Geometry.Origin(1);
                    this.y2Min = Geometry.Origin(2);
                end
                this.y1Max      = this.y1Min + Geometry.L1; 
                this.y2Max      = this.y2Min + Geometry.L2;            
            else
                this.y1Min = Geometry.y1Min;
                this.y1Max = Geometry.y1Max;                
                this.y2Min = Geometry.y2Min;
                this.y2Max = Geometry.y2Max;
            end
            
            this.polar = 'cart';
            
            InitializationPts(this);                
        end
    end
    
    methods (Access = public)


        function [y1_kv,y2_kv,J,dH1,dH2] = PhysSpace(this,x1,x2)
            %[y1_kv,y2_kv,J,DDy1,DDy2] = Comp_to_Phys(x1,x2)
            %[y1_kv,dy1dx1] = LinearMap(x1,-L1,L1); %,dx1,ddx1,dddx,ddddx]                
            n  = length(x1);
            c1 = 0.4;
            c3 = 0.05;

            [y1_kv,Diffy1] =  M1Tref(x1,x2*c1+c3,0.07);                 
            [y2_kv,dy2dx2] =  LinearMap(x2,this.y2Min,this.y2Max);

            if(nargout >= 3)
                J        = zeros(n,2,2);
                J(:,1,1) = Diffy1.dydx; 
                J(:,1,2) = Diffy1.dyd_d*c1;           
                J(:,2,1) = zeros(n,1);
                J(:,2,2) = dy2dx2;            
            end

            if(nargout >= 4)
                dH1      = zeros(n,2,2);
                dH1(:,1,1) = Diffy1.dyddx; %d2y1dx1dx1; 
                dH1(:,1,2) = c1*Diffy1.dydxd_d;%d2y1dx1dx2;             
                dH1(:,2,1) = c1*Diffy1.dydxd_d;%d2y1dx1dx2;     
                dH1(:,2,2) = (c1^2).*Diffy1.dydd_d;
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
             %xf  = InvQuadMapAB(z,L1,-L1,L1);
             error = MException('InfCapillary_Interface',...
                                'inverse function not yet implemented');
             throw(error);
        end             
        
    end
end