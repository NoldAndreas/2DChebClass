classdef DiffuseInterfaceSingleFluid < DiffuseInterface
         
    methods (Access = public) 
         function this = DiffuseInterfaceSingleFluid(config)           
             this@DiffuseInterface(config);
         end
                 
         IterationStepFullProblem(this,noIterations)                
             
         function [v_cont,A_cont] = Continuity(this,uv,phi,G)
            
            Cn             = this.optsPhys.Cn;
            Cak            = this.optsPhys.Cak; 
            Ind            = this.IDC.Ind;
            zeta           = this.optsPhys.zeta;
                        
            lbC            = Ind.left & Ind.bottom;
            rbC            = Ind.right & Ind.bottom;
            y2Max          = this.optsNum.PhysArea.y2Max;    
            IntPathUpLow   = this.IDC.borderTop.IntSc - this.IDC.borderBottom.IntSc; 
            
            IP                      = this.IDC.SubShapePtsCart(this.RightCapillary.GetCartPts);           
            IntPath_Half_UpLow      = (this.RightCapillary.borderTop.IntSc - this.RightCapillary.borderBottom.IntSc)*IP;
            IntPath_Half_RightLeft  = (this.RightCapillary.borderRight.IntSc - this.RightCapillary.borderLeft.IntSc)*IP;
            
                        
            Diff           = this.IDC.Diff;            
            phi_m          = this.optsPhys.phi_m; 
            M              = this.IDC.M;
            F              = false(M,1);   
            T              = true(M,1);
            uvTdiag        = [diag(uv(1:end/2)),diag(uv(end/2+1:end))];
            
            [fWP,fW,fWPP]  = DoublewellPotential(phi,Cn);
            
            
            % Continuity   %[uv;phi;G]
            A_cont         = [diag(phi+phi_m)*Diff.div + [diag(Diff.Dy1*phi),diag(Diff.Dy2*phi)]...
                              diag(Diff.div*uv)+uvTdiag*Diff.grad,...
                              zeros(M)];
            v_cont         = Diff.div*(uv.*repmat(phi+phi_m,2,1));
                        
            
            % ************[uv;phi;G]
            ys             = fWP - Cn*Diff.Lap*phi;    
            ysP            = diag(fWPP) - Cn*Diff.Lap;
            Cmu            = -Diff.Lap*diag(phi+phi_m);
            
            Cuv            = Cak*(zeta + 4/3)*Diff.LapDiv; 
            
            A_cont(Ind.top|Ind.bottom,[T;T;F;F])  = Cuv(Ind.top|Ind.bottom,:);
            A_cont_phi                            = -Diff.Lap*diag(G) ...
                                                    + Diff.div*(diag(repmat(ys,2,1))*Diff.grad)...
                                                    + Diff.div*(diag(Diff.grad*phi)*repmat(ysP,2,1));
                                                    
            A_cont(Ind.top|Ind.bottom,[F;F;T;F])  = A_cont_phi(Ind.top|Ind.bottom,:);
            A_cont(Ind.top|Ind.bottom,[F;F;F;T])  = Cmu(Ind.top|Ind.bottom,:);
            
            v_cont(Ind.top|Ind.bottom)    =   Cmu(Ind.top|Ind.bottom,:)*G ...
                                            + Cuv(Ind.top|Ind.bottom,:)*uv ...
                                            + Diff.div(Ind.top|Ind.bottom,:)*(repmat(ys,2,1).*(Diff.grad*phi));
                                        
            % (BC1) y1 = +/- infinity
            A_cont(Ind.left|Ind.right,:)         = 0;
            A_cont(Ind.left|Ind.right,[F;F;F;T]) = Diff.Dy2(Ind.left|Ind.right,:);
            v_cont(Ind.left|Ind.right,:)         = Diff.Dy2(Ind.left|Ind.right,:)*G;                                        
            
            % ************            
            
            A_m11 = zeros(M,4*M);
            A_m11(:,[T;T;F;F]) = Cak*((zeta -2/3)*Diff.div + 2*[Diff.Dy1,zeros(M,M)]);
            A_m11(:,[F;F;T;F]) = diag(fWP - G) + Cn*( -diag(Diff.Dy1*phi)*Diff.Dy1 + diag(Diff.Dy2*phi).*Diff.Dy2);
            A_m11(:,[F;F;F;T]) = - diag(phi+phi_m);
            
            m11 =  Cak*((zeta -2/3)*Diff.div*uv + 2*Diff.Dy1*uv(1:end/2))...
                 + (fW - G.*(phi+phi_m) + Cn/2*( -(Diff.Dy1*phi).^2 + (Diff.Dy2*phi).^2 ));
            
            
            A_m12 = zeros(M,4*M);
            A_m12(:,[T;T;F;F]) = Cak*[Diff.Dy2 , Diff.Dy1];
            A_m12(:,[F;F;T;F]) = -Cn*(diag(Diff.Dy1*phi)*Diff.Dy2 + diag(Diff.Dy2*phi)*Diff.Dy1);
            
            m12 = Cak*[Diff.Dy2 , Diff.Dy1]*uv - Cn*((Diff.Dy1*phi).*(Diff.Dy2*phi));
            
            
            A_cont(lbC,:) = IntPath_Half_UpLow*A_m12 +  IntPath_Half_RightLeft*A_m11;
            v_cont(lbC)   = IntPath_Half_UpLow*m12   +  IntPath_Half_RightLeft*m11;            
            
            % (BC2) [uv;phi;G]                   
            A_cont(rbC,:)   = IntPathUpLow*A_m12 + y2Max*(A_m11(rbC,:) - A_m11(lbC,:));
            v_cont(rbC)     = IntPathUpLow*m12 +y2Max*(m11(rbC) - m11(lbC));
                 
         end        
         function [v_mom,A_mom]   = Momentum(this,uv,phi,G)                
            %[uv;phi;G;p] 
            %
            %        A_mom          = [tauM,...
            %                          -[diag(Diff.Dy1*G);diag(Diff.Dy2*G)],...
            %                          - diag(repmat(phi+phi_m,2,1))*Diff.grad];
            % 
            %        v_mom         = tauM*uv - repmat(phi+phi_m,2,1).*(Diff.grad*G);

            Ind            = this.IDC.Ind;
            Diff           = this.IDC.Diff;            
            phi_m          = this.optsPhys.phi_m;             
            Cak            = this.optsPhys.Cak;    
            Cn             = this.optsPhys.Cn;
            zeta           = this.optsPhys.zeta;
            M              = this.IDC.M;
            
            [fWP,fW,fWPP]  = DoublewellPotential(phi,Cn);
            tauM           = Cak*(Diff.LapVec + (zeta + 1/3)*Diff.gradDiv);            
            phiM2          = repmat(phi+phi_m,2,1);
           
            A_mom_phi      = -[diag(Diff.Dy1*G);diag(Diff.Dy2*G)]...
                             +diag(Diff.grad*phi)*repmat(diag(fWPP)-Cn*Diff.Lap,2,1)...
                             +diag(repmat(fWP - Cn*Diff.Lap*phi - G,2,1))*Diff.grad;       
            A_mom          = [tauM,...
                              A_mom_phi,...
                              - diag(phiM2)*Diff.grad - [diag(Diff.Dy1*phi);diag(Diff.Dy2*phi)]];
            v_mom         = tauM*uv - phiM2.*(Diff.grad*G) + ...
                             repmat(fWP - Cn*Diff.Lap*phi - G,2,1).*(Diff.grad*phi);
                                                              

            
         end        
         function [v_SeppAdd,A_SeppAdd] = GetSeppecherConditions(this,uv,phi,G,a,deltaX,theta)

            M              = this.IDC.M;
            Ind            = this.IDC.Ind;
            Diff           = this.IDC.Diff;
            IntNormalUp    = this.IDC.borderTop.IntNormal;
            Cn             = this.optsPhys.Cn;
            phi_m          = this.optsPhys.phi_m;    
            phiM2          = repmat(phi+phi_m,2,1);       
            y2Max          = this.optsNum.PhysArea.y2Max;
            lbC            = Ind.left & Ind.bottom;
            rbC            = Ind.right & Ind.bottom;
            F              = false(M,1);
            T              = true(M,1);   
 
            [fWP,fW,fWPP]  = DoublewellPotential(phi,Cn);

            % Three extra conditions [uv;phi;G]
            % (EX 1) int((phi+rho_m)*u_y|_y2Max,y1=-infty..infty) = 2*y2Max        
            A_a            = zeros(1,4*M);        

            A_a([T;T;F;F])   = IntNormalUp*diag(phiM2);
            A_a([F;F;T;F])   = IntNormalUp*[diag(uv(1:end/2));diag(uv(1+end/2:end))];
            A_a([F;F;rbC;F]) = A_a([F;F;rbC;F]) + y2Max;
            A_a([F;F;lbC;F]) = A_a([F;F;lbC;F]) - y2Max;
            A_a              = [0,0,0,A_a];        

            v_a              = IntNormalUp*(phiM2.*uv) + (phi(rbC)-phi(lbC))*y2Max;


            % (EX 2) phi(y2Max/tan(theta) + deltaX,y2Max) = 0
            InterpMatchPos       = this.IDC.SubShapePtsCart(...
                                    struct('y1_kv',deltaX + y2Max/tan(theta),...
                                           'y2_kv',y2Max));
            A_deltaX             = zeros(1,4*M);
            A_deltaX([F;F;T;F])  = InterpMatchPos;
            A_deltaX_deltaX      = InterpMatchPos*(Diff.Dy1*phi);
            A_deltaX_theta       = -y2Max*(1/sin(theta))^2*InterpMatchPos*(Diff.Dy1*phi);
            A_deltaX             = [0,A_deltaX_deltaX,A_deltaX_theta,...
                                    A_deltaX];    
            v_deltaX             = InterpMatchPos*phi;

            % (EX 3) mu(y1=-infty) = 0
             A_theta              = zeros(1,4*M);
             A_theta([F;F;lbC;F]) = -G(lbC)+fWP(lbC);
             A_theta([F;F;F;lbC]) = -(phi(lbC)+phi_m);        
             A_theta              = [0,0,0,A_theta];
             v_theta              = -(phi(lbC)+phi_m)*G(lbC) + fW(lbC);

            %        A_theta = zeros(1,4*M);
            %        A_theta = [0,0,1,A_theta];    
            %        v_theta = theta - pi/2;

             v_SeppAdd = [v_a;v_deltaX;v_theta];
             A_SeppAdd = [A_a;A_deltaX;A_theta];
         end
        
         %Old
         [A,b]       = FullStressTensorIJ(this,phi,i,j)                    
    end
end