classdef DiffuseInterfaceBinaryFluid < DiffuseInterface
    
   properties (Access = public)
       p  % pressure    
       flux
	end
       
	methods (Access = public) 
       function this = DiffuseInterfaceBinaryFluid(config)     
            if(nargin == 0)
                 config = [];
            end
            if(~isempty(config) && ~isfield(config.optsPhys,'mobility') && ...
                    isfield(config.optsPhys,'l_diff'))
                config.optsPhys.mobility = (config.optsPhys.l_diff)^2/(config.optsPhys.Cak);
            end
            this@DiffuseInterface(config);            
       end
       
       IterationStepFullProblem(this)        
       
       function SetCak(this,Cak)
           if(isfield(this.optsPhys,'l_diff'))
               this.optsPhys.mobility = (this.optsPhys.l_diff)^2/Cak;
           end            
           this.optsPhys.Cak = Cak;           
           SaveConfig(this);   
       end
       function vec = GetInitialCondition(this,theta) 
            if(nargin == 1)
                theta = pi/2;
            end
            vec = GetInitialCondition@DiffuseInterface(this,theta);           
            if(length(vec) == this.IDC.M*4+3)
                vec = [0;vec];
            end
            if(isempty(this.p))
                this.p   = zeros(this.IDC.M,1);                
            end           
           vec      = [vec;this.p];
       end        
       function [v_mom,A_mom]         = GetVelBC(this,uv,a,deltaX,theta)
           M = this.IDC.M;           
           if(nargin == 2)
                [v_mom,A_mom] = GetVelBC@DiffuseInterface(this,uv);
           else
                [v_mom,A_mom] = GetVelBC@DiffuseInterface(this,uv,a,deltaX,theta);
           end
           A_mom         = [A_mom,zeros(size(A_mom,1),M)]; 
       end
       function [v_mu_T,A_mu_T]       = GetPhiBC(this,phi,theta)
           M = this.IDC.M;           
           [v_mu_T,A_mu_T] = GetPhiBC@DiffuseInterface(this,phi,theta);                      
           A_mu_T          = [A_mu_T,zeros(size(A_mu_T,1),M)];
           if(IsSeppecher(this))
               A_mu_T = [zeros(size(A_mu_T,1),1),A_mu_T];
           end
       end       
       function [v_G_IBB,A_G_IBB]     = GetChemPotBC(this,uv,phi,G,theta)
           
            m              = this.optsPhys.mobility;
            Ind            = this.IDC.Ind;
            Diff           = this.IDC.Diff;
            M    = this.IDC.M;
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
            %A_G(Ind.left & Ind.bottom,[F;T;F;F;F]) = this.IDC.borderBottom.IntSc*(diag(phi)*Diff.Dy2);
            A_G(Ind.left & Ind.bottom,[F;F;T;F;F]) = + EYM(Ind.right & Ind.bottom,:) ...
                                                     - EYM(Ind.left & Ind.bottom,:) ;%...
                                                     %+ this.IDC.borderBottom.IntSc*(diag(Diff.Dy2*(uv(1+end/2:end))));
            A_G(Ind.left & Ind.bottom,[F;F;F;T;F]) = - m*this.IDC.borderBottom.IntSc*Diff.DDy2;
            v_G(Ind.left & Ind.bottom)             =    phi(Ind.right & Ind.bottom) ...
                                                      - phi(Ind.left & Ind.bottom) ...                                                     
                                                      - m*this.IDC.borderBottom.IntSc*Diff.DDy2*G;%+ this.IDC.borderBottom.IntSc*(phi.*(Diff.Dy2*(uv(1+end/2:end))))...
                              
            % ****************
            IP                      = this.IDC.SubShapePtsCart(this.RightCapillary.GetCartPts);           
            IP0                     = this.IDC.SubShapePtsCart(struct('y1_kv',0,'y2_kv',0));
            IntPath_Half_Low        = this.RightCapillary.borderBottom.IntSc*IP;            
            
            A_G(Ind.right & Ind.bottom,:)           = 0;
            %A_G(Ind.right & Ind.bottom,[F;T;F;F;F]) = IntPath_Half_Low*(diag(phi)*Diff.Dy2);
            A_G(Ind.right & Ind.bottom,[F;F;T;F;F]) = + EYM(Ind.right & Ind.bottom,:) - IP0 ;%...
                                                      %+ IntPath_Half_Low*(Diff.Dy2*(uv(1+end/2:end)));
            A_G(Ind.right & Ind.bottom,[F;F;F;T;F]) = - m*IntPath_Half_Low*Diff.DDy2;
            
            v_G(Ind.right & Ind.bottom) = phi(Ind.right & Ind.bottom) - IP0*phi ...                                          
                                          - m*IntPath_Half_Low*Diff.DDy2*G;%+ IntPath_Half_Low*(phi.*(Diff.Dy2*(uv(1+end/2:end))))...
                                      
            if(IsSeppecher(this))
                A_G  = [A_G_four,A_G];
            end                                      
            
            A_G_IBB = A_G(Ind.bound,:);
            v_G_IBB = v_G(Ind.bound);  
            
        end    
       function [v_SeppAdd,A_SeppAdd] = GetSeppecherConditions(this,uv,phi,G,p,deltaX,theta)

            M              = this.IDC.M;
            Ind            = this.IDC.Ind;
            Diff           = this.IDC.Diff;
            IntNormalUp    = this.IDC.borderTop.IntNormal;
            Cn             = this.optsPhys.Cn;
            m              = this.optsPhys.mobility;
            y2Max          = this.optsNum.PhysArea.y2Max;
            lbC            = Ind.left & Ind.bottom;
            rbC            = Ind.right & Ind.bottom;
            F              = false(M,1);
            T              = true(M,1);   
 
            [fWP,fW,fWPP]  = DoublewellPotential(phi,Cn);

            % Three extra conditions 
            % [uv;phi;G,p]
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
            InterpMatchPos       = this.IDC.SubShapePtsCart(...
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
                   
          %   A_theta = [0,0,0,1,zeros(1,5*M)];    
          %   v_theta = theta - pi/2;

             v_SeppAdd = [v_a;v_b;v_deltaX;v_theta];
             A_SeppAdd = [A_a;A_b;A_deltaX;A_theta];
             %A_SeppAdd =  [A_SeppAdd,zeros(size(A_SeppAdd,1),M)];
       end 
       function [v_mom_IBB,A_mom_IBB] = GetSeppecher_Vel(this,uv,a,deltaX,theta)

            a1 = a(1);
            a2 = a(2);            
            
            Ind            = this.IDC.Ind;
            Diff           = this.IDC.Diff;
            M              = this.IDC.M;
            PtsCart        = this.IDC.GetCartPts();
            y2Max          = this.IDC.y2Max;
            ITT            = repmat(Ind.top,2,1);  
            IBB            = repmat(Ind.bound,2,1); 
        	Dy12           = blkdiag(Diff.Dy1,Diff.Dy1);
            d_theta        = 0.01;            
            F              = false(M,1);   
            T              = true(M,1);   
            Z              = zeros(M,1);
            
            %*****************************************
            %Get beta from
            % theta = beta + Ca*f(beta)           
            [beta,dthetadbeta] = GetBeta(this,theta);            
            
            %Get Delta Z
            deltaZ         = deltaX + y2Max*(1/tan(theta) - 1/tan(beta));
            
            u_flow         = GetSeppecherSolutionCart([PtsCart.y1_kv - deltaZ,PtsCart.y2_kv],1,0,0,beta);                  
            %u_flow         = GetSeppecherSolutionCart_Blurred([PtsCart.y1_kv - deltaX,PtsCart.y2_kv],1,0,0,theta);                  
            %*****************************************
                                            
            y1_Interface  = PtsCart.y1_kv - (deltaX + y2Max/tan(theta));
            
            d = 4;
            
            f_a1  = (y1_Interface/d).*exp(-(y1_Interface/d).^2);
            fP_a1 = ((1-2/d*(y1_Interface.^2))/d).*exp(-(y1_Interface/d).^2);            
            %f_a1  = y1_Interface./(1+y1_Interface.^2);
            %fP_a1 = (1-y1_Interface.^2)./(1+y1_Interface.^2).^2;            
            f_a1(Ind.left | Ind.right) = 0;            
            fP_a1(Ind.left | Ind.right) = 0;
           
            f_a2  = exp(-(y1_Interface/d).^2);
            fP_a2 = -2*y1_Interface/d.*f_a2;
            %f_a2  = 1./(1+y1_Interface.^2);
            %fP_a2 = -2*y1_Interface./(1+y1_Interface.^2).^2;
            fP_a2(Ind.left | Ind.right) = 0;
                       
            a1_corr        = a1*f_a1;   
            a1_corr_a1     = f_a1;
            a1_corr_deltaX = -a1*fP_a1;
            a1_corr_theta  = a1*fP_a1*y2Max/(sin(theta)).^2;                                
            
            a2_corr        = a2*f_a2;
            a2_corr_a2     = f_a2;  
            a2_corr_deltaX = -a2*fP_a2;
            a2_corr_deltaX(Ind.left | Ind.right) = 0;
            a2_corr_theta  = a2*fP_a2*y2Max/(sin(theta)).^2;            
                        
            corr = a1_corr + a2_corr; %1+
    
            uvBound        = u_flow + [Z;corr];

            uvBound_a1      = [Z;a1_corr_a1];     %u_flow.*
            uvBound_a2      = [Z;a2_corr_a2];     %u_flow.*
            uvBound_deltaX = -(Dy12*u_flow)+...%.*repmat(corr,2,1)+...
                               [Z;a1_corr_deltaX+a2_corr_deltaX]; %u_flow .*

            %u_flow_PTheta    = GetSeppecherSolutionCart([PtsCart.y1_kv - deltaX,PtsCart.y2_kv],1,0,0,theta+d_theta);
            dbeta            = d_theta/dthetadbeta;
            %deltaZ          = deltaX + y2Max*(1/tan(theta) - 1/tan(beta));
            ddeltaZ          = deltaX + y2Max*(1/tan(theta+d_theta) - 1/tan(beta+dbeta));

            u_flow_PTheta    = GetSeppecherSolutionCart([PtsCart.y1_kv - ddeltaZ,PtsCart.y2_kv],1,0,0,beta+dbeta);
            u_flow_d         = (u_flow_PTheta - u_flow)/d_theta;
            uvBound_theta    = u_flow_d +... %repmat(corr,2,1)
                               [Z;a1_corr_theta + a2_corr_theta];  %u_flow.*

            A_mom_ITT        = [-uvBound_a1,...
                                -uvBound_a2,...
                                -uvBound_deltaX,...
                                -uvBound_theta,...
                                eye(2*M),zeros(2*M,2*M)];
                                  
            A_mom(ITT,:)     = A_mom_ITT(ITT,:);                                 
            v_mom(ITT)       = uv(ITT) - uvBound(ITT);
            
            % Boundary condition at the wall and at +/- infinity
            u_Wall = [ones(M,1);zeros(M,1)];
        
            A_mom(IBB & ~ITT,:)                                        = 0;
            A_mom(IBB & ~ITT,[false;false;false;false;IBB & ~ITT;F;F]) = eye(sum(IBB & ~ITT));
            v_mom(IBB & ~ITT)                                          = uv(IBB & ~ITT) - u_Wall(IBB & ~ITT);
            
            A_mom_IBB = A_mom(IBB,:);
            v_mom_IBB = v_mom(IBB);
            
       end
       function [beta,dthetadbeta] = GetBeta(this,theta)
           if(nargin < 2)
               theta = this.theta;
           end

           Ca          = 3/4*this.optsPhys.Cak;
           
           fsolveOpts  = optimset('Display','off');
           beta        = fsolve(@func,theta,fsolveOpts);
           [~,dfdbeta] = f_stokes(beta,1);
           dthetadbeta = 1 + Ca*dfdbeta;
           
           function y = func(b)
                y = b + Ca*f_stokes(b,1) - theta;
            end
       end

       vec_a  = FindAB(this,phi,mu,deltaX,a_ig)              
       
       function [v_cont,A_cont] = Continuity(this,uv,phi,G,p)            
           
           Ind  = this.IDC.Ind;
           Diff = this.IDC.Diff;
           M    = this.IDC.M;
           EYM  = eye(M);
           Z    = zeros(M);
           
           F    = false(M,1);   
           T    = true(M,1);    
           
           Cn          = this.optsPhys.Cn;
           Cak         = this.optsPhys.Cak;
           y2Max       = this.optsNum.PhysArea.y2Max;
                      
           IntPathUpLow   = this.IDC.borderTop.IntSc - this.IDC.borderBottom.IntSc;     
           [fWP,fW,fWPP]  = DoublewellPotential(phi,Cn);
           
           A_cont         = [Diff.div,Z,Z,Z];    
           v_cont         = Diff.div*uv;
           
           % (BC0) [uv;phi;G;p]
%            indBC0 = Ind.top;
%            A_cont(indBC0,[T;T;F;F;F]) = Diff.LapDiv(indBC0,:);
%            A_cont(indBC0,[F;F;T;F;F]) = Diff.div(indBC0,:)*(diag(repmat(G,2,1))*Diff.grad)/Cak;
%            A_cont(indBC0,[F;F;F;T;F]) = Diff.div(indBC0,:)*([diag(Diff.Dy1*phi);diag(Diff.Dy2*phi)])/Cak;
%            A_cont(indBC0,[F;F;F;F;T]) = - Diff.Lap(indBC0,:);
%            
%            v_cont(indBC0) = - Diff.Lap(indBC0,:)*p...
%                                           + Diff.LapDiv(indBC0,:)*uv ...
%                                           + Diff.div(indBC0,:)*(repmat(G,2,1).*(Diff.grad*phi))/Cak;
            
           % (BC1) p = 0  at y1 = +/- infinity
           A_cont(Ind.left,:)           = 0;
           A_cont(Ind.left,[F;F;F;F;T]) = EYM(Ind.left,:);
           v_cont(Ind.left)             = p(Ind.left);

            % (BC2) [uv;phi;G;p]
           A_cont(Ind.right,:)           = 0;
           A_cont(Ind.right,[F;F;F;F;T]) = Diff.Dy2(Ind.right,:);
           v_cont(Ind.right,:)           = Diff.Dy2(Ind.right,:)*p;

           indBC2 = Ind.right & Ind.bottom;

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
       function [v_mu,A_mu]     = ChemicalPotential(this,phi,G)
           [v_mu,A_mu] = ChemicalPotential@DiffuseInterface(this,phi,G);
           A_mu = [A_mu,zeros(this.IDC.M)];
       end
       function [v_G,A_G]       = PhasefieldEq(this,uv,phi,G)
           
            m              = this.optsPhys.mobility;            
            Diff           = this.IDC.Diff;
            M    = this.IDC.M;            
            Z    = zeros(M);            
            
            uvTdiag        = [diag(uv(1:end/2)),diag(uv(end/2+1:end))];        
            gradPhiMT      = [diag(Diff.Dy1*phi),diag(Diff.Dy2*phi)];
              
            A_G            = [-gradPhiMT,-uvTdiag*Diff.grad,m*Diff.Lap,Z];
            v_G            = m*Diff.Lap*G - gradPhiMT*uv;           
       end
       function [v_mom,A_mom]   = Momentum(this,uv,phi,G,p)                
           
            %[uv;phi;G;p]         
            Diff = this.IDC.Diff;
            Cak  = this.optsPhys.Cak;            
           
            A_mom          = [Diff.LapVec,...
                               diag(repmat(G,2,1))*Diff.grad/Cak,...                               
                               [diag(Diff.Dy1*phi);diag(Diff.Dy2*phi)]/Cak,...
                               - Diff.grad];
            v_mom          = Diff.LapVec*uv ...
                              + repmat(G,2,1).*(Diff.grad*phi)/Cak ...
                              - Diff.grad*p;
       end       

       function [y1M,y2M,fl_y1,fl_y2,startMask1,startMask2] = PlotFlux(this)
            [y1M,y2M,fl_y1,fl_y2,startMask1,startMask2] = PlotU(this,this.flux); 
       end
       function PlotResults(this)           
           PlotResults@DiffuseInterface(this);
           
           figure('Position',[0 0 800 600],'color','white');
           PlotU(this); hold on;            
           this.IDC.plot(this.p,'contour');     
           SaveCurrentFigure(this,'Pressure');
       end       
       function [y1M,y2M,fl_y1,fl_y2,startMask1,startMask2] = PlotU(this,uv,opts) 
           
            if(nargin < 3)
                %opts = [];
                opts.linewidth = 1.4;
                opts.color = 'k';
            end
           
            if(nargin == 1)
                uv = [];
            end
            
            if(isempty(this.StagnationPoint))
               FindStagnationPoint(this,[-2,2],[-2,this.IDC.y2Max-2]); 
            end            
            
            y1Pts = [this.optsNum.PlotArea.y1Min+0.1;...
                     this.optsNum.PlotArea.y1Max-0.1];%0;...
            y2Pts = [0.5;...
                     this.optsNum.PlotArea.y2Max-0.5]; %this.optsNum.PlotArea.y2Max/2;...
            
            [y1M,y2M,fl_y1,fl_y2,startMask1,startMask2] = PlotU@DiffuseInterface(this,uv,y1Pts,y2Pts,opts);
       end       
       function PlotComponentVelocities(this)
           
           uv       = this.uv;
           mu       = this.mu;
           phi      = repmat(this.phi,2,1);
           
           m        = this.optsPhys.mobility;
           Diff     = this.IDC.Diff;
           
           flux     = m*(Diff.grad*mu);
           
           n      = 20;
           y2Min  = this.optsNum.PhysArea.y2Min;
            y2Max = this.optsNum.PhysArea.y2Max;
            
            y1Min = this.optsNum.PlotArea.y1Min;
            y1Max = this.optsNum.PlotArea.y1Max;
                        
            y2L = y2Min + (y2Max-y2Min)*(0:n-1)'/(n-1);            
            y1L = y1Min + (y1Max-y1Min)*(0:n-1)'/(n-1);            
            
            startPtsy1    = [y1Max*ones(size(y2L))-0.1;...
                         y1Min*ones(size(y2L))+0.1;...
                         y1L];
            startPtsy2    = [y2L;y2L;y2Max*ones(size(y1L))];
           
           uv_1 = uv - flux./(1+phi);
           uv_2 = uv + flux./(1-phi);
           
           PlotRightArea  = struct('N',[40,40],...
                                   'y1Min',0,'y1Max',10,...
                                   'y2Min',0,'y2Max',this.IDC.y2Max);   	
                       
           bxRight        = Box(PlotRightArea);    
           PlotRightArea.N1 = 90;
           PlotRightArea.N2 = 90;
           bxRight.ComputeIndices();
           bxRight.InterpolationPlotCart(PlotRightArea,true);     
           IPRight        = this.IDC.SubShapePtsCart(bxRight.GetCartPts());
           IPR2           = blkdiag(IPRight,IPRight);
           
           uv_1_Right     = IPR2*uv - (IPR2*flux)./(1+IPR2*phi);
           
           figure('color','white','Position',[0 0 1000 1000]);
           opts = struct('color','g','linewidth',2.5);
           bxRight.plotStreamlines(uv_1_Right,startPtsy1,startPtsy2,opts); hold on;
           bxRight.plotFlux(uv_1_Right); hold on;
           this.IDC.plot(this.phi,'contour');        hold on;                                      
 %          PlotU(this);
          % PlotU(this,flux,struct('color','m','linewidth',2.5));
           
           
           
           
          % subplot(2,1,2);
%           this.IDC.plot(this.phi,'contour');                
%%           AddStreamlines(this,uv_1);           
%           AddStreamlines(this,uv_2);
%           PlotU(this,uv_2,opts);                      
       end
       function CheckResultResolution(this)
           figure('Position',[0 0 1000 1000],'color','white');
           
           y2M = this.IDC.y2Max;
           
           subplot(2,2,1);           
           this.IDC.plotLine([-20 20],[0 0],this.IDC.Diff.Lap*this.mu);           
           this.IDC.plotLine([-20 20],[y2M y2M]/2,this.IDC.Diff.Lap*this.mu);
           this.IDC.plotLine([-20 20],[y2M y2M],this.IDC.Diff.Lap*this.mu);
           title('$\Delta \mu$','Interpreter','Latex','fontsize',15);
           
           subplot(2,2,2);
           this.IDC.plotLine([-20 20],[0 0],this.mu);           
           this.IDC.plotLine([-20 20],[y2M y2M]/2,this.mu);
           this.IDC.plotLine([-20 20],[y2M y2M],this.mu);
           title('$\mu$','Interpreter','Latex','fontsize',15);
           
           subplot(2,2,3);
           this.IDC.plotLine([-20 20],[0 0],this.IDC.Diff.Lap*this.uv(1:end/2));
           this.IDC.plotLine([-20 20],[y2M y2M]/2,this.IDC.Diff.Lap*this.uv(1:end/2));
           this.IDC.plotLine([-20 20],[y2M y2M],this.IDC.Diff.Lap*this.uv(1:end/2));
           title('$\Delta u$','Interpreter','Latex','fontsize',15);           
           
           subplot(2,2,4);
           this.IDC.plotLine([-20 20],[0 0],this.IDC.Diff.Lap*this.uv(1+end/2:end));
           this.IDC.plotLine([-20 20],[y2M y2M]/2,this.IDC.Diff.Lap*this.uv(1+end/2:end));
           this.IDC.plotLine([-20 20],[y2M y2M],this.IDC.Diff.Lap*this.uv(1+end/2:end));      
           title('$\Delta v$','Interpreter','Latex','fontsize',15);
           
       end
       function PostProcess_Flux(this)
           Diff      = this.IDC.Diff;
           m         = this.optsPhys.mobility;
           uv        = this.uv;           
           %this.flux = repmat(1+this.phi,2,1).*uv/2 - m*(Diff.grad*this.mu)/2; %Particles 1
           %this.flux = repmat(1-this.phi,2,1).*uv/2 + m*(Diff.grad*this.mu)/2; %Particles 2
           %this.flux = -m*(Diff.grad*this.mu);            
           this.flux = repmat(this.phi,2,1).*uv - m*(Diff.grad*this.mu); 
           
           disp(['Error of flux through subArea: ',num2str(this.Int_of_path*uv)]);
           disp(['Error of phasefield flux through subArea: ',num2str(this.Int_of_path*this.flux)]);
       end
       function PostProcessInterface(this)
            PostProcessInterface@DiffuseInterface(this);
            
            pts.y1_kv = this.IsoInterface.h;
            pts.y2_kv = this.IDC.Pts.y2;            
            IP        = this.IDC.SubShapePtsCart(pts);            
            
            this.IsoInterface.p     = IP*this.p;      
            
            if(~isempty(this.flux))               
                
                fl_1      = this.flux(1:end/2);
                fl_2      = this.flux(1+end/2:end);
                th        = this.IsoInterface.theta;
               % this.IsoInterface.flux_n   = cos(th).*(IP*fl_1)  + sin(th).*(IP*fl_2);
               % this.IsoInterface.flux_t   = -sin(th).*(IP*fl_1) + cos(th).*(IP*fl_2);                  
                this.IsoInterface.flux_n   = sin(th).*(IP*fl_1)  - cos(th).*(IP*fl_2);
                this.IsoInterface.flux_t   = cos(th).*(IP*fl_1) + sin(th).*(IP*fl_2);                                  
            end
       end
       function PostProcess(this)
           GetYueParameters(this);           
           FindStagnationPoint(this,[-2,2],[-2,this.IDC.y2Max-2]); 
           PostProcess_Flux(this);
           PostProcessInterface(this);
       end
       function GetYueParameters(this)
           if(IsSeppecher(this))
                Cak      = this.optsPhys.Cak;
                Cn       = this.optsPhys.Cn;
                mobility = this.optsPhys.mobility;
                S        = sqrt(Cak*mobility);
                ls_Yue   = 2.5*S*Cn;
                disp(['Diffusion length/L: ',num2str(S)]);
                disp(['Prediction of slip length / distance of stagnation point to wall by Yue et. al (2010): ',num2str(ls_Yue)]);
           end            
           
       end
       function Plot3DProfiles(this)
           figure('color','white','Position',[0 0 800 800]);
           this.IDC.plot(this.mu,{'SC'});
           PlotU(this);
           title('Chemical Potential');
           SaveCurrentFigure(this,'3DChemPot');
           
           figure('color','white','Position',[0 0 800 800]);
           this.IDC.plot(this.phi,{'SC'});
           PlotU(this);
           title('Phasefield');
           SaveCurrentFigure(this,'3DPhasefield');
           
           figure('color','white','Position',[0 0 800 800]);
           this.IDC.plot(this.p,{'SC'});
           PlotU(this);
           title('Pressure');
           SaveCurrentFigure(this,'3DPressure');
       end      
       function PlotValuesThroughStagnationPoint(this)
           
           figure('Position',[0 0 800 800],'color','white');
           y1Int = 3*[this.optsNum.PlotArea.y1Min,this.optsNum.PlotArea.y1Max];
           %y2Int =  this.StagnationPoint.y2_kv(1)*[1,1];
           y2Int =  (this.StagnationPoint.y2_kv(1) - 2)*[1,1];
           y2Int =  0*[1,1];
           this.IDC.plotLine(y1Int,y2Int,100*this.mu,struct('color','b'));
           this.IDC.plotLine(y1Int,y2Int,this.p,struct('color','r'));
           this.IDC.plotLine(y1Int,y2Int,this.uv(1:end/2),struct('color','m'));
           this.IDC.plotLine(y1Int,y2Int,this.uv(1+end/2:end),struct('color','k'));
           set(gca,'linewidth',1.5);
           set(gca,'fontsize',15);
           SaveCurrentFigure(this,'ValuesThroughStagnationPoint');
           
       end
   end
end