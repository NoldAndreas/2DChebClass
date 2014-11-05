classdef DiffuseInterfaceBinaryFluid < DiffuseInterface
    
   properties (Access = public)
       p  % pressure                   
   end
       
	methods (Access = public) 
       function this = DiffuseInterfaceBinaryFluid(config)           
             this@DiffuseInterface(config);
       end
       
       IterationStepFullProblem(this,noIterations)       
            
       function vec = GetInitialCondition(this) 
            vec = GetInitialCondition@DiffuseInterface(this);           
            %if(length(vec) == this.IC.M*4+3)
%                vec = [0;vec];
%            end
            if(isempty(this.p))
                this.p   = zeros(this.IC.M,1);                
            end           
           vec      = [vec;this.p];
       end        
       function [v_mom,A_mom] = GetVelBC(this,uv,a,deltaX,theta)
           M = this.IC.M;           
           if(nargin == 2)
                [v_mom,A_mom] = GetVelBC@DiffuseInterface(this,uv);
           else
                [v_mom,A_mom] = GetVelBC@DiffuseInterface(this,uv,a,deltaX,theta);
           end
           A_mom         = [A_mom,zeros(size(A_mom,1),M)]; 
       end
       function [v_mu_T,A_mu_T] = GetPhiBC(this,phi,theta)
           M = this.IC.M;           
           [v_mu_T,A_mu_T] = GetPhiBC@DiffuseInterface(this,phi,theta);                      
           A_mu_T    = [A_mu_T,zeros(size(A_mu_T,1),M)];
           if(IsSeppecher(this))
               A_mu_T = [zeros(size(A_mu_T,1),1),A_mu_T];
           end
       end       
       function [v_G_IBB,A_G_IBB] = GetChemPotBC(this,uv,phi,G,theta)
           
            m              = this.optsPhys.mobility;
            Ind            = this.IC.Ind;
            Diff           = this.IC.Diff;
            M    = this.IC.M;
            EYM  = eye(M); 
            Z    = zeros(M);
            F       = false(M,1);   T       = true(M,1);    
           
            % (BC4) nu*grad(G) = 0        
            A_G                         = zeros(M,5*M);
            v_G                         = zeros(M,1);
            
            A_G(Ind.bottom|Ind.top,[F;F;F;T;F]) = Diff.Dy2(Ind.bottom|Ind.top,:);    
            v_G(Ind.bottom|Ind.top)             = Diff.Dy2(Ind.bottom|Ind.top,:)*G;       
            
            if(IsSeppecher(this))
                 a_direction               = [-sin(theta)*EYM,cos(theta)*EYM];
                 a_direction_theta         = -[cos(theta)*EYM,sin(theta)*EYM];            

                 A_G(Ind.top,[F;F;F;T;F])  = a_direction(Ind.top,:)*Diff.grad;
                 
                 A_G_four                                    = zeros(M,4);                 
                 A_G_four(Ind.top,[false;false;false;true])  = a_direction_theta(Ind.top,:)*(Diff.grad*G);                 
                 v_G(Ind.top)                                = a_direction(Ind.top,:)*(Diff.grad*G);                                                                                                                    
                %A_G(Ind.top,:)           = 0;
                %A_G(Ind.top,[F;F;F;T;F]) = Diff.Dy1(Ind.top,:);  
                %v_G(Ind.top)             = Diff.Dy1(Ind.top,:)*G;                
            end

            %(BC7) 
            indBC7 = Ind.left | Ind.right;
            A_G(indBC7,:)           = 0;
            A_G(indBC7,[F;F;F;T;F]) = Diff.Dy2(indBC7,:);
            v_G(indBC7,:)           = Diff.Dy2(indBC7,:)*G;

            %(BC8) [uv;phi;G;p]                    
            A_G(Ind.left & Ind.bottom,:)           = 0;
            A_G(Ind.left & Ind.bottom,[F;T;F;F;F]) = this.IC.borderBottom.IntSc*(diag(phi)*Diff.Dy2);
            A_G(Ind.left & Ind.bottom,[F;F;T;F;F]) = + EYM(Ind.right & Ind.bottom,:) ...
                                                      - EYM(Ind.left & Ind.bottom,:) ...                                              
                                                      + this.IC.borderBottom.IntSc*(diag(Diff.Dy2*(uv(1+end/2:end))));
            A_G(Ind.left & Ind.bottom,[F;F;F;T;F]) = - m*this.IC.borderBottom.IntSc*Diff.DDy2;
            v_G(Ind.left & Ind.bottom)             =    phi(Ind.right & Ind.bottom) ...
                                                      - phi(Ind.left & Ind.bottom) ...
                                                      + this.IC.borderBottom.IntSc*(phi.*(Diff.Dy2*(uv(1+end/2:end))))...
                                                      - m*this.IC.borderBottom.IntSc*Diff.DDy2*G;
                              
            % ****************
            IP                      = this.IC.SubShapePtsCart(this.RightCapillary.GetCartPts);           
            IP0                     = this.IC.SubShapePtsCart(struct('y1_kv',0,'y2_kv',0));
            IntPath_Half_Low        = this.RightCapillary.borderBottom.IntSc*IP;            
            
            A_G(Ind.right & Ind.bottom,:)           = 0;
            A_G(Ind.right & Ind.bottom,[F;T;F;F;F]) = IntPath_Half_Low*(diag(phi)*Diff.Dy2);
            A_G(Ind.right & Ind.bottom,[F;F;T;F;F]) = + EYM(Ind.right & Ind.bottom,:) - IP0 ...                                              
                                                      + IntPath_Half_Low*(Diff.Dy2*(uv(1+end/2:end)));
            A_G(Ind.right & Ind.bottom,[F;F;F;T;F]) = - m*IntPath_Half_Low*Diff.DDy2;
            
            v_G(Ind.right & Ind.bottom) = phi(Ind.right & Ind.bottom) - IP0*phi ...
                                          + IntPath_Half_Low*(phi.*(Diff.Dy2*(uv(1+end/2:end))))...
                                          - m*IntPath_Half_Low*Diff.DDy2*G;
                                      
            if(IsSeppecher(this))
                A_G  = [A_G_four,A_G];
            end                                      
            
            A_G_IBB = A_G(Ind.bound,:);
            v_G_IBB = v_G(Ind.bound);  
            
        end    
       function [v_SeppAdd,A_SeppAdd] = GetSeppecherConditions(this,uv,phi,G,p,deltaX,theta)

            M              = this.IC.M;
            Ind            = this.IC.Ind;
            Diff           = this.IC.Diff;
            IntNormalUp    = this.IC.borderTop.IntNormal;
            Cn             = this.optsPhys.Cn;
            m              = this.optsPhys.mobility;
            y2Max          = this.optsNum.PhysArea.y2Max;
            lbC            = Ind.left & Ind.bottom;
            rbC            = Ind.right & Ind.bottom;
            F              = false(M,1);
            T              = true(M,1);   
 
            [fWP,fW,fWPP]  = DoublewellPotential(phi,Cn);

            % Three extra conditions [uv;phi;G,p]
            % (EX 1) int((phi+rho_m)*u_y|_y2Max,y1=-infty..infty) = 2*y2Max        
            % Conservation of phase 1 matter
            A_a                = zeros(1,5*M);        
            A_a([T;T;F;F;F])   = IntNormalUp*diag(repmat(phi,2,1));
            A_a([F;F;T;F;F])   = IntNormalUp*[diag(uv(1:end/2));diag(uv(1+end/2:end))];
            A_a([F;F;rbC;F;F]) = A_a([F;F;rbC;F;F]) + y2Max;
            A_a([F;F;lbC;F;F]) = A_a([F;F;lbC;F;F]) - y2Max;            
            A_a([F;F;F;T;F])   = -m*IntNormalUp*Diff.grad;
            A_a                = [0,0,0,0,A_a];        
            
            v_a                = IntNormalUp*(repmat(phi,2,1).*uv - m*Diff.grad*G) +(phi(rbC)-phi(lbC))*y2Max;

            % Conservation of mass
            A_b                = zeros(1,5*M);        
            A_b([T;T;F;F;F])   = IntNormalUp;
            A_b                = [0,0,0,0,A_b];        
            v_b                = IntNormalUp*uv;

            % (EX 2) phi(y2Max/tan(theta) + deltaX,y2Max) = 0
            InterpMatchPos       = this.IC.SubShapePtsCart(...
                                    struct('y1_kv',deltaX + y2Max/tan(theta),...
                                           'y2_kv',y2Max));
            A_deltaX               = zeros(1,5*M);
            A_deltaX([F;F;T;F;F])  = InterpMatchPos;
            A_deltaX_deltaX        = InterpMatchPos*(Diff.Dy1*phi);
            A_deltaX_theta         = -y2Max*(1/sin(theta))^2*InterpMatchPos*(Diff.Dy1*phi);
            A_deltaX               = [0,0,A_deltaX_deltaX,A_deltaX_theta,...
                                      A_deltaX];    
            v_deltaX               = InterpMatchPos*phi;

            % (EX 3) mu(y1=-infty) = 0
             %A_theta              = zeros(1,4*M);
             %A_theta([F;F;F;lbC]) = 1;
             %A_theta              = [0,0,0,A_theta];
             %v_theta              = G(lbC);
              A_theta              = zeros(1,5*M);
              A_theta([F;F;F;F;rbC]) = 1;
              A_theta              = [0,0,0,0,A_theta];
              v_theta              = p(rbC);
                   
             %A_theta = [0,0,0,1,zeros(1,5*M)];    
             %v_theta = theta - pi/2;

             v_SeppAdd = [v_b;v_a;v_deltaX;v_theta];
             A_SeppAdd = [A_b;A_a;A_deltaX;A_theta];
             %A_SeppAdd =  [A_SeppAdd,zeros(size(A_SeppAdd,1),M)];
        end
       
       function SavePlotResults(this)           
           SavePlotResults@DiffuseInterface(this);
           
           figure('Position',[0 0 800 600],'color','white');
           PlotU(this); hold on;            
           this.IC.doPlots(this.p,'contour');     
           SaveCurrentFigure(this,[this.filename '_Pressure']);                
       end       
       function PlotU(this,uv,opts) 
           
            if(nargin < 2)
                opts = [];
            end
           
            if(nargin == 1)
                uv = [];
            end
            
            if(isempty(this.StagnationPoint))
               FindStagnationPoint(this,[-2,2],[-2,this.IC.y2Max-2]); 
            end            
            
            y1Pts = [0;...
                     this.optsNum.PlotArea.y1Min+0.1;...
                     this.optsNum.PlotArea.y1Max-0.1];
            y2Pts = [this.optsNum.PlotArea.y2Max/2;...
                     0.5;...
                     this.optsNum.PlotArea.y2Max-0.5];
            
            PlotU@DiffuseInterface(this,uv,y1Pts,y2Pts,opts);
        end       
       function CheckResultResolution(this)
           figure('Position',[0 0 1000 1000],'color','white');
           
           y2M = this.IC.y2Max;
           
           subplot(2,2,1);           
           this.IC.doPlotFLine([-20 20],[0 0],this.IC.Diff.Lap*this.mu);           
           this.IC.doPlotFLine([-20 20],[y2M y2M]/2,this.IC.Diff.Lap*this.mu);
           this.IC.doPlotFLine([-20 20],[y2M y2M],this.IC.Diff.Lap*this.mu);
           title('$\Delta \mu$','Interpreter','Latex','fontsize',15);
           
           subplot(2,2,2);
           this.IC.doPlotFLine([-20 20],[0 0],this.mu);           
           this.IC.doPlotFLine([-20 20],[y2M y2M]/2,this.mu);
           this.IC.doPlotFLine([-20 20],[y2M y2M],this.mu);
           title('$\mu$','Interpreter','Latex','fontsize',15);
           
           subplot(2,2,3);
           this.IC.doPlotFLine([-20 20],[0 0],this.IC.Diff.Lap*this.uv(1:end/2));
           this.IC.doPlotFLine([-20 20],[y2M y2M]/2,this.IC.Diff.Lap*this.uv(1:end/2));
           this.IC.doPlotFLine([-20 20],[y2M y2M],this.IC.Diff.Lap*this.uv(1:end/2));
           title('$\Delta u$','Interpreter','Latex','fontsize',15);           
           
           subplot(2,2,4);
           this.IC.doPlotFLine([-20 20],[0 0],this.IC.Diff.Lap*this.uv(1+end/2:end));
           this.IC.doPlotFLine([-20 20],[y2M y2M]/2,this.IC.Diff.Lap*this.uv(1+end/2:end));
           this.IC.doPlotFLine([-20 20],[y2M y2M],this.IC.Diff.Lap*this.uv(1+end/2:end));      
           title('$\Delta v$','Interpreter','Latex','fontsize',15);
           
       end
       
       function [v_cont,A_cont] = Continuity(this,uv,phi,G,p)            
           
           Ind  = this.IC.Ind;
           Diff = this.IC.Diff;
           M    = this.IC.M;
           EYM  = eye(M);
           Z    = zeros(M);
           
           F    = false(M,1);   
           T    = true(M,1);    
           
           Cn          = this.optsPhys.Cn;
           Cak         = this.optsPhys.Cak;
           y2Max       = this.optsNum.PhysArea.y2Max;
                      
           IntPathUpLow   = this.IC.borderTop.IntSc - this.IC.borderBottom.IntSc;     
           [fWP,fW,fWPP]  = DoublewellPotential(phi,Cn);
           
           A_cont         = [Diff.div,Z,Z,Z];    
           v_cont         = Diff.div*uv;
           
           % (BC0) [uv;phi;G;p]
           indBC0 = Ind.top;
           A_cont(indBC0,[T;T;F;F;F]) = Diff.LapDiv(indBC0,:);
           A_cont(indBC0,[F;F;T;F;F]) = Diff.div(indBC0,:)*(diag(repmat(G,2,1))*Diff.grad)/Cak;
           A_cont(indBC0,[F;F;F;T;F]) = Diff.div(indBC0,:)*([diag(Diff.Dy1*phi);diag(Diff.Dy2*phi)])/Cak;
           A_cont(indBC0,[F;F;F;F;T]) = - Diff.Lap(indBC0,:);
           
           v_cont(indBC0) = - Diff.Lap(indBC0,:)*p...
                                          + Diff.LapDiv(indBC0,:)*uv ...
                                          + Diff.div(indBC0,:)*(repmat(G,2,1).*(Diff.grad*phi))/Cak;
            
           % (BC1) p = 0  at y1 = +/- infinity
           A_cont(Ind.left,:)           = 0;
           A_cont(Ind.left,[F;F;F;F;T]) = EYM(Ind.left,:);
           v_cont(Ind.left)             = p(Ind.left);

            % (BC2) [uv;phi;G;p]
           
           A_cont(Ind.right,:)           = 0;
           A_cont(Ind.right,[F;F;F;F;T]) = Diff.Dy2(Ind.right,:);
           v_cont(Ind.right,:)           = Diff.Dy2(Ind.right,:)*p;
           
           indBC2 = Ind.right & Ind.top;
           
           A_m11 = zeros(M,5*M);            
           A_m11(:,[F;F;T;F;F]) = diag(fWP);
           A_m11(:,[F;F;F;F;T]) = -Cak*eye(M);
           m11 =  -p*Cak + fW;                 
           
           A_m12 = zeros(M,5*M);
           A_m12(:,[T;T;F;F;F]) = Cak*[Diff.Dy2 , Diff.Dy1];
           A_m12(:,[F;F;T;F;F]) = -Cn*(diag(Diff.Dy1*phi)*Diff.Dy2 + diag(Diff.Dy2*phi)*Diff.Dy1); 
           m12 = Cak*[Diff.Dy2 , Diff.Dy1]*uv - Cn*((Diff.Dy1*phi).*(Diff.Dy2*phi));                       
           
           A_cont(indBC2,:) = IntPathUpLow*A_m12  ...
                                + y2Max*(A_m11(Ind.right & Ind.top,:) - ...
                                         A_m11(Ind.left & Ind.top,:));           
           v_cont(indBC2) = IntPathUpLow*m12...
                                +y2Max*(m11(Ind.right & Ind.top) ...
                                      - m11(Ind.left & Ind.top));                                

       end
       function [v_mu,A_mu] = ChemicalPotential(this,phi,G)
           [v_mu,A_mu] = ChemicalPotential@DiffuseInterface(this,phi,G);
           A_mu = [A_mu,zeros(this.IC.M)];
       end
       function [v_G,A_G] = Phasefield(this,uv,phi,G)
           
            m              = this.optsPhys.mobility;            
            Diff           = this.IC.Diff;
            M    = this.IC.M;            
            Z    = zeros(M);            
            
            uvTdiag        = [diag(uv(1:end/2)),diag(uv(end/2+1:end))];        
            gradPhiMT      = [diag(Diff.Dy1*phi),diag(Diff.Dy2*phi)];
              
            A_G            = [-gradPhiMT,-uvTdiag*Diff.grad,m*Diff.Lap,Z];
            v_G            = m*Diff.Lap*G - gradPhiMT*uv;           
       end
       function [v_mom,A_mom]   = Momentum(this,uv,phi,G,p)                
           
            %[uv;phi;G;p]         
            Diff = this.IC.Diff;
            Cak  = this.optsPhys.Cak;            
           
            A_mom          = [Diff.LapVec,...
                               diag(repmat(G,2,1))*Diff.grad/Cak,...                               
                               [diag(Diff.Dy1*phi);diag(Diff.Dy2*phi)]/Cak,...
                               - Diff.grad];
            v_mom          = Diff.LapVec*uv ...
                              + repmat(G,2,1).*(Diff.grad*phi)/Cak ...
                              - Diff.grad*p;
       end
       
       function [v_mom_IBB,A_mom_IBB] = GetSeppecher_Vel(this,uv,a,deltaX,theta)

            b = a(1);
            a = a(2);
            %b = 0;
            
            Ind            = this.IC.Ind;
            Diff           = this.IC.Diff;
            M              = this.IC.M;
            PtsCart        = this.IC.GetCartPts();
            y2Max          = this.IC.y2Max;
            ITT            = repmat(Ind.top,2,1);  
            IBB            = repmat(Ind.bound,2,1); 
        	Dy12           = blkdiag(Diff.Dy1,Diff.Dy1);
            d_theta        = 0.01;
            EYMM           = eye(2*M);
            F              = false(M,1);   
            T              = true(M,1);   
            EYM            = eye(M);
            
            u_flow         = GetSeppecherSolutionCart_Blurred([PtsCart.y1_kv - deltaX,...
                                                PtsCart.y2_kv],1,0,0,theta);                  
                                            
            y1_Interface  = PtsCart.y1_kv - (deltaX + y2Max/tan(theta));
            
%             b_corr        = b*y1_Interface./(1+y1_Interface.^2);
%             b_corr(Ind.left | Ind.right) = 0;
%             b_corr_b      = y1_Interface./(1+y1_Interface.^2);        
%             b_corr_b(Ind.left | Ind.right) = 0;
%             b_corr_deltaX = -b*(1-y1_Interface.^2)./(1+y1_Interface.^2).^2;
%             b_corr_deltaX(Ind.left | Ind.right) = 0;
%             b_corr_theta  = -b_corr_deltaX*y2Max/(sin(theta)).^2;            


            b_corr        = b./(1+2*y1_Interface.^2);
            b_corr_b      = 1./(1+2*y1_Interface.^2);        
            b_corr_deltaX = 4*b*y1_Interface./(1+2*y1_Interface.^2).^2;
            b_corr_deltaX(Ind.left | Ind.right) = 0;
            b_corr_theta  = -b_corr_deltaX*y2Max/(sin(theta)).^2;            
            
            
            a_corr        = a./(1+y1_Interface.^2);
            a_corr_a      = 1./(1+y1_Interface.^2);        
            a_corr_deltaX = 2*a*y1_Interface./(1+y1_Interface.^2).^2;
            a_corr_deltaX(Ind.left | Ind.right) = 0;
            a_corr_theta  = -a_corr_deltaX*y2Max/(sin(theta)).^2;            
                        
            corr = 1 + a_corr + b_corr;
    
            uvBound        = u_flow .*repmat(corr,2,1);                                        

            uvBound_a      = u_flow.*repmat(a_corr_a,2,1);               
            uvBound_b      = u_flow.*repmat(b_corr_b,2,1);               
            uvBound_deltaX = -(Dy12*u_flow).*repmat(corr,2,1)+...
                               u_flow .*repmat(a_corr_deltaX+b_corr_deltaX,2,1);

            u_flow_PTheta    = GetSeppecherSolutionCart([PtsCart.y1_kv - deltaX,PtsCart.y2_kv],1,0,0,theta+d_theta);
            u_flow_d         = (u_flow_PTheta - u_flow)/d_theta;
            uvBound_theta    = u_flow_d.*repmat(corr,2,1) +...
                               u_flow.*repmat(a_corr_theta + b_corr_theta,2,1); 

            A_mom_ITT              = zeros(sum(ITT),4*M);            
            A_mom_ITT(:,[T;T;F;F]) = EYMM(ITT,:);            
            
            A_mom_a_ITT            = -uvBound_a(ITT);
            A_mom_b_ITT            = -uvBound_b(ITT);
            A_mom_deltaX_ITT       = -uvBound_deltaX(ITT); 
            A_mom_theta_ITT        = -uvBound_theta(ITT);        

            A_mom(ITT,:)           = [A_mom_b_ITT,...
                                      A_mom_a_ITT,...
                                      A_mom_deltaX_ITT,...
                                      A_mom_theta_ITT,...
                                      A_mom_ITT];
                                 
            v_mom(ITT)              = uv(ITT) - uvBound(ITT);
            
            u_Wall = [ones(M,1);zeros(M,1)];
        
            A_mom(IBB & ~ITT,:)                = 0;
            A_mom(IBB & ~ITT,[false;false;false;false;IBB & ~ITT;F;F]) = eye(sum(IBB & ~ITT));
            v_mom(IBB & ~ITT)                  = uv(IBB & ~ITT) - u_Wall(IBB & ~ITT);
            
            A_mom_IBB = A_mom(IBB,:);
            v_mom_IBB = v_mom(IBB);
        end
              
   end
end