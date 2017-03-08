classdef M1Box < M1SpectralSpectral

    properties        
        L1,L2                
    end
    
    methods        
        function this = M1Box(Geometry,N)
            this@M1SpectralSpectral(N(1),N(2));
            
            this.L1      = Geometry.L1; 
            this.L2      = Geometry.L2; 
            
            if(isfield(Geometry,'Origin'))
                this.Origin = Geometry.Origin;
            end
            
            InitializationPts(this);            
        end                         
    end
    
    methods (Access = public)
    
        function [y1_kv,y2_kv,J,dH1,dH2] = Comp_to_Phys(x1,x2)
            %[y1_kv,y2_kv,J,DDy1,DDy2] = Comp_to_Phys(x1,x2)
            %[y1_kv,dy1dx1] = LinearMap(x1,-L1,L1); %,dx1,ddx1,dddx,ddddx]                
            n  = length(x1);
            c1 = 0.4;
            c3 = 0.05;

            [y1_kv,Diffy1] =  M1Tref(x1,x2*c1+c3,0.07);                 
            [y2_kv,dy2dx2] =  LinearMap(x2,-L2,L2);

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
      
        
    end
end