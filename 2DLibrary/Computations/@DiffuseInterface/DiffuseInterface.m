classdef DiffuseInterface < Computation
    
   properties (Access = public)                                          
       IntSubArea
              
       phi = [],uv  = [],mu=[]      
       theta=[],a=[],deltaX=[]
       errors = []
       
       IsoInterface=struct('y2',[],'h',[],'theta',[],'kappa',[],...
                                 'mu',[],'p',[],'u_n',[],'u_t',[]);
       StagnationPoint=[]
       %surfaceTension = 4/3;
       
       RightCapillary
       filename
   end
   
   methods (Abstract = true, Access = public)
       [A,b] = FullStressTensorIJ(this,phi,i,j)       
   end
   
   methods (Access = public)          
        function this = DiffuseInterface(config)  
            if(nargin == 0)
                 config = [];
            end
            this@Computation(config);            
            if(~isfield(this.optsNum.PhysArea,'y1Max'))
                this.optsNum.PhysArea.y1Max = inf;
            end      
            if(~isfield(this.optsNum.PhysArea,'NBorder'))
                this.optsNum.PhysArea.NBorder = [];
            end
            
            if(~isfield(this.optsNum.PlotArea,'y2Min'))
                this.optsNum.PlotArea.y2Min = this.optsNum.PhysArea.y2Min;
            end
            if(~isfield(this.optsNum.PlotArea,'y2Max'))
                this.optsNum.PlotArea.y2Max = this.optsNum.PhysArea.y2Max;
            end
        end       
        function Preprocess(this)          
            
            if(~isfield(this.optsNum.PhysArea,'shape'))
                this.optsNum.PhysArea.shape = 'InfCapillaryQuad';
            end
            Preprocess@Computation(this);                        
            
            Sel = {'Dy1' ;'DDy1' ; 'Dy2'; 'DDy2';...
                    'DDDy2';'DDDy1';...
                    'Dy1Dy2'; 'DDy1Dy2'; 'Dy1DDy2';...
                    'Lap' ;'grad' ;'div'; ...
                    'gradLap' ;'gradDiv'; 'LapVec';'LapDiv';'Lap2'};
            this.IDC.ComputeDifferentiationMatrix(Sel);
                                    
            this.IDC.SetUpBorders(this.optsNum.PhysArea.NBorder);            
            
            this.IDC.Ind.fluidInterface = [];
            if(isfield(this.optsPhys,'fluidInterface'))
                this.IDC.Ind.fluidInterface = this.IDC.Ind.(this.optsPhys.fluidInterface);
            end

            bxArea          = struct('y1Min',this.optsNum.PhysArea.IntInterval(1),...
                                     'y1Max',this.optsNum.PhysArea.IntInterval(2),...
                                     'y2Min',this.optsNum.PhysArea.y2Min,...
                                     'y2Max',this.optsNum.PhysArea.y2Max,...
                                     'N',[100,100]);
            
            BX              = Box(bxArea);
            IntBx           = BX.ComputeIntegrationVector();
            this.IntSubArea = IntBx*this.IDC.SubShapePts(BX.GetCartPts());
                     
            
                         Phys_Area = struct('L1',this.optsNum.PhysArea.L1,...
                                'N',[50;40],...
                                'y2Min',this.optsNum.PhysArea.y2Min,...
                                'y2Max',this.optsNum.PhysArea.y2Max);
    
            this.RightCapillary  = HalfInfCapillary(Phys_Area);
            this.RightCapillary.SetUpBorders(100);    
        end       
        
        function sepp  = IsSeppecher(this)
            sepp = (length(this.optsPhys.thetaEq) == 1);
        end        
        function vec = GetInitialCondition(this,theta) 
            if(nargin == 1)
                theta = pi/2;
            end
           
            if(~isempty(this.uv) && ...                
                ~isempty(this.mu) && ...
                ~isempty(this.phi))                
                cprintf('r*','Existing data is ignored in initial Condition.\n')               
            end                       
            
            this.phi = InitialGuessRho(this,theta);                             
            this.uv  = InitialGuessUV(this,theta);
            this.mu  = zeros(this.IDC.M,1);                          
            vec      = [this.uv;this.phi;this.mu];                         
                        
            if(IsSeppecher(this))
                if(isempty(this.theta))
                    vec = [0;0;theta;vec];
                else
                    vec = [0;0;0;theta;vec];
                    %vec = [this.a;this.deltaX;this.theta;vec];
                end
            end
             
         end
        function phi = InitialGuessRho(this,theta)
            if(nargin == 1)
                theta = pi/2;
            end
            PtsCart    = this.IDC.GetCartPts();
            Cn         = this.optsPhys.Cn;
            phi        = tanh((PtsCart.y1_kv-PtsCart.y2_kv/tan(theta))/Cn);            
        end                         
        function uv  = InitialGuessUV(this,theta)
            if(nargin == 1)
                theta = pi/2;
            end
            UWall            = this.optsPhys.UWall;    
            PtsCart          = this.IDC.GetCartPts();            
            y2_kv            = PtsCart.y2_kv;
            y2Min            = this.optsNum.PhysArea.y2Min;
            y2Max            = this.optsNum.PhysArea.y2Max;
            M                = this.IDC.M;
            
            if(IsSeppecher(this))
                deltaX = 0;                 
                uv     = GetSeppecherSolutionCart_Blurred([PtsCart.y1_kv - deltaX,...
                                                PtsCart.y2_kv],1,0,0,theta);                                
            else                
                uv     = [UWall(1) + ...
                            (UWall(2)-UWall(1))*(y2_kv-y2Min)/(y2Max-y2Min);...
                            zeros(M,1)];                
            end
        end
                                          
        function [v_mu_TB,A_mu_TB] = GetPhiBC(this,phi,theta)
            
            %[uv;phi;G]
            Ind            = this.IDC.Ind;
            Diff           = this.IDC.Diff;
            M              = this.IDC.M;            
            F              = false(M,1);   
            T              = true(M,1);   
            EYM            = eye(M);
            
            % (BC4.a) nu*grad(phi) = 0            
            A_mu                               = zeros(M,4*M);
            A_mu(Ind.bottom|Ind.top,[F;F;T;F]) = Diff.Dy2(Ind.bottom|Ind.top,:);    
            v_mu(Ind.bottom|Ind.top)           = Diff.Dy2(Ind.bottom|Ind.top,:)*phi;            
            

            %A_momThree            = [A_mom_a',A_mom_deltaX',A_mom_theta'];            

            if(IsSeppecher(this))
                a_direction               = [cos(theta)*EYM,sin(theta)*EYM];            
                a_direction_theta         = [-sin(theta)*EYM,cos(theta)*EYM];            
            
                A_mu                                      = [zeros(M,3),A_mu];
                A_mu(Ind.top,[false;false;false;F;F;T;F]) = a_direction(Ind.top,:)*Diff.grad;
                A_mu(Ind.top,[false;false;true;F;F;F;F])  = a_direction_theta(Ind.top,:)*(Diff.grad*phi);
                
                v_mu(Ind.top)                             = a_direction(Ind.top,:)*(Diff.grad*phi);                                                   
            end
            
            A_mu_TB = A_mu(Ind.top | Ind.bottom,:);
            v_mu_TB = v_mu(Ind.top | Ind.bottom);
            
        end     
        function [v_mom,A_mom] = GetVelBC(this,uv,a,deltaX,theta)
            
            UWall            = this.optsPhys.UWall;    
            PtsCart          = this.IDC.GetCartPts();            
            y2_kv            = PtsCart.y2_kv;
            y2Min            = this.optsNum.PhysArea.y2Min;
            y2Max            = this.optsNum.PhysArea.y2Max;
            M                = this.IDC.M;
            EYMM             = eye(2*M);            
            Ind              = this.IDC.Ind;
            IBB              = repmat(Ind.bound,2,1); 
            F                = false(M,1);   
            T                = true(M,1);
            
            if(IsSeppecher(this))
                [v_mom,A_mom] = GetSeppecher_Vel(this,uv,a,deltaX,theta);
            else                
                uvBound    = [UWall(1) + ...
                             (UWall(2)-UWall(1))*(y2_kv-y2Min)/(y2Max-y2Min);...
                              zeros(M,1)];                
                                                
                A_mom                = zeros(sum(IBB),4*M);
                A_mom(:,[T;T;F;F])   = EYMM(IBB,:);
                v_mom                = uv(IBB) - uvBound(IBB);
            end            
        end                                         
        function [v_mom_IBB,A_mom_IBB] = GetSeppecher_Vel(this,uv,a,deltaX,theta)
            
            Ind            = this.IDC.Ind;
            Diff           = this.IDC.Diff;
            M              = this.IDC.M;
            PtsCart        = this.IDC.GetCartPts();
            y2Max          = this.IDC.y2Max;
            ITT            = repmat(Ind.top,2,1);
            IBB            = repmat(Ind.bound,2,1); 
        	Dy12           = blkdiag(Diff.Dy1,Diff.Dy1);
            d_theta        = 0.01;
            EYMM           = eye(2*M);
            F              = false(M,1);   
            T              = true(M,1);   
            EYM            = eye(M);
            Cak            = this.optsPhys.Cak;
            
            %Get beta from 
            % theta = beta + Ca*f(beta)
            Ca          = 3/4*Cak;
            beta        = fsolve(@func,theta);
            [~,fb]      = f_stokes(beta);
            dthetadbeta = 1 + Ca*fb;
            
            %Get Delta Z
            deltaZ         = deltaX + y2Max*(1/tan(theta) - 1/tan(beta));
            
            u_flow         = GetSeppecherSolutionCart([PtsCart.y1_kv - deltaZ,PtsCart.y2_kv],1,0,0,beta);                  
            
            %u_flow         = GetSeppecherSolutionCart([PtsCart.y1_kv - deltaX,...%_Blurred
%                                                PtsCart.y2_kv],1,0,0,theta);                  
                                            
            y1_Interface                        = PtsCart.y1_kv - (deltaX + y2Max/tan(theta));
            a_corr                              = (1 + a./(1+y1_Interface.^2));
            a_corr_a                            = 1./(1+y1_Interface.^2);        
            a_corr_deltaX                       = 2*a*y1_Interface./(1+y1_Interface.^2).^2;
            a_corr_deltaX(Ind.left | Ind.right) = 0;
            a_corr_theta                        = -a_corr_deltaX*y2Max/(sin(theta)).^2;            
    
            uvBound          = u_flow .*repmat(a_corr,2,1);                                        

            uvBound_a        = u_flow.*repmat(a_corr_a,2,1);               
            uvBound_deltaX   = -(Dy12*u_flow).*repmat(a_corr,2,1)+u_flow .*repmat(a_corr_deltaX,2,1);        

            u_flow_PTheta    = GetSeppecherSolutionCart([PtsCart.y1_kv - deltaZ,PtsCart.y2_kv],1,0,0,beta+d_theta/dthetadbeta);
            u_flow_d         = (u_flow_PTheta - u_flow)/d_theta;
            uvBound_theta    = u_flow_d.*repmat(a_corr,2,1) + u_flow.*repmat(a_corr_theta,2,1);

            A_mom_ITT              = zeros(sum(ITT),4*M);            
            A_mom_ITT(:,[T;T;F;F]) = EYMM(ITT,:);            

            A_mom_a_ITT            = -uvBound_a(ITT);
            A_mom_deltaX_ITT       = -uvBound_deltaX(ITT); 
            A_mom_theta_ITT        = -uvBound_theta(ITT);        

            A_mom(ITT,:)           = [A_mom_a_ITT,...
                                     A_mom_deltaX_ITT,...
                                     A_mom_theta_ITT,...
                                     A_mom_ITT];
                                 
            v_mom(ITT)              = uv(ITT) - uvBound(ITT);
            
            u_Wall = [ones(M,1);zeros(M,1)];
        
            A_mom(IBB & ~ITT,:)                = 0;
            A_mom(IBB & ~ITT,[false;false;false;IBB & ~ITT;F;F]) = eye(sum(IBB & ~ITT));
            v_mom(IBB & ~ITT)                  = uv(IBB & ~ITT) - u_Wall(IBB & ~ITT);
            
            A_mom_IBB = A_mom(IBB,:);
            v_mom_IBB = v_mom(IBB);
            
            function t = func(b)
                t = b + Ca*f_stokes(b);
            end
                
        end
        function [v_mu,A_mu] = ChemicalPotential(this,phi,G)
            %[uv;phi;G]     
            Cn             = this.optsPhys.Cn;
            M              = this.IDC.M;
            Z              = zeros(M);
            Diff           = this.IDC.Diff;            
            
            [fWP,fW,fWPP]  = DoublewellPotential(phi,Cn);
            
            A_mu           = [Z,Z,diag(fWPP)-Cn*Diff.Lap,-eye(M)];
            v_mu           = fWP - Cn*Diff.Lap*phi - G;                           
        end  
        
        %Auxiliary Cahn-Hilliard
        [A,b] = Div_FullStressTensor(this,phi)
                
        %PostProcessing
        function PostProcessInterface(this)
            y2           = this.IDC.Pts.y2;
            Diff         = barychebdiff(y2,2);            
            phi          = this.phi; 

            fsolveOpts   = optimset('Display','off');
            interface    = zeros(size(y2));

            y1Bottom = this.IDC.GetCartPts.y1_kv(this.IDC.Ind.bottom);
            [~,j]    = min(abs(phi(this.IDC.Ind.bottom)));
            y1I      = y1Bottom(j);
                        
            for i = 1:length(y2)
                pt.y2_kv     = y2(i);        
                [interface(i),~,flag] = fsolve(@phiX1,y1I,fsolveOpts);        
                if(flag < 1)
                    cprintf('*r',['PostProcessInterface: Interface not found at y2 = ',num2str(y2(i)),'\n']);
                    return;
                end
                y1I          = interface(i);        
            end

            this.IsoInterface.h = interface;    
            
            %Compute values on Interface 
            pts.y1_kv = this.IsoInterface.h;
            pts.y2_kv = y2;            
            IP        = this.IDC.SubShapePtsCart(pts);                                    
            u1        = this.uv(1:end/2);
            u2        = this.uv(1+end/2:end);            
            
            this.IsoInterface.IP      = IP;
            this.IsoInterface.y2      = y2;
            y10                       = this.IsoInterface.h(1);
            this.IsoInterface.r       = sqrt(y2.^2 + (pts.y1_kv-y10));
            this.IsoInterface.mu      = IP*this.mu;                        
            this.IsoInterface.mu_ddy2 = IP*(this.IDC.Diff.DDy2*this.mu);     
                        
            th                      = pi/2 - atan(Diff.Dx*interface);
            this.IsoInterface.theta = th;
            
            this.IsoInterface.hatL  = FitSliplength(this);
            
            this.IsoInterface.kappa = (Diff.DDx*interface)./((1+(Diff.Dx*interface).^2).^(3/2));   
            this.IsoInterface.u_n   = sin(th).*(IP*u1)  - cos(th).*(IP*u2);
            this.IsoInterface.u_t   = cos(th).*(IP*u1) + sin(th).*(IP*u2);                         
            
            %Compute jump across interface
            Cn     = this.optsPhys.Cn;
            DeltaZ = 1.5*Cn;
            alpha  = th-pi/2;
            pts_right.y1_kv = pts.y1_kv + DeltaZ*cos(alpha);
            pts_right.y2_kv = pts.y2_kv + DeltaZ*sin(alpha);
            
            pts_left.y1_kv  = pts.y1_kv - DeltaZ*cos(alpha);
            pts_left.y2_kv  = pts.y2_kv - DeltaZ*sin(alpha);
            
            IP_Jump = this.IDC.SubShapePtsCart(pts_left) - this.IDC.SubShapePtsCart(pts_right);
            
            
            J_u1 = IP_Jump*this.uv(1:end/2);
            J_u2 = IP_Jump*this.uv(1+end/2:end);
            this.IsoInterface.Jump_u_1  = J_u1;
            this.IsoInterface.Jump_u_2  = J_u2;
            this.IsoInterface.Jump_u_n  = sin(th).*J_u1 - cos(th).*J_u2;
            this.IsoInterface.Jump_u_t  = cos(th).*J_u1 + sin(th).*J_u2;
                        
            J_dmudy1                           = IP_Jump*(this.IDC.Diff.Dy1*this.mu);
            J_dmudy2                           = IP_Jump*(this.IDC.Diff.Dy2*this.mu);
            this.IsoInterface.Jump_DmuDy1      = J_dmudy1;
            this.IsoInterface.Jump_DmuDy2      = J_dmudy2;
            this.IsoInterface.Jump_Gradmu_dnu  = sin(th).*J_dmudy1 - cos(th).*J_dmudy2;
            this.IsoInterface.Jump_Gradmu_dt   = cos(th).*J_dmudy1 + sin(th).*J_dmudy2;                         
            
            this.IsoInterface.Jump_mu = IP_Jump*this.mu;
            this.IsoInterface.Jump_p  = IP_Jump*this.p;
            
            %[A11,b11] = FullStressTensorIJ(this,this.phi,1,1);
            %tau11 = A11*[this.p;this.uv] + b11;
            %[A12,b12] = FullStressTensorIJ(this,this.phi,1,2);
            %tau12 = A12*[this.p;this.uv] + b12;
            %[A22,b22] = FullStressTensorIJ(this,this.phi,2,2);
            %tau22 = A22*[this.p;this.uv] + b22;
            
            tau11 = this.IDC.Diff.Dy1*u1;
            tau12 = 1/2*(this.IDC.Diff.Dy1*u2+this.IDC.Diff.Dy2*u1);
            tau22 = this.IDC.Diff.Dy2*u2;
            
            Jump_tau11 = IP_Jump*tau11;
            Jump_tau12 = IP_Jump*tau12;
            Jump_tau22 = IP_Jump*tau22;
                        
            this.IsoInterface.Jump_tau11 = Jump_tau11;
            this.IsoInterface.Jump_tau12 = Jump_tau12;
            this.IsoInterface.Jump_tau22 = Jump_tau22;
            
            S = sin(alpha); C = cos(alpha);
            this.IsoInterface.Jump_tau_nn = (C.^2).*Jump_tau11 + 2*C.*S.*Jump_tau12 + (S.^2).*Jump_tau22;
            this.IsoInterface.Jump_tau_nt = -C.*S.*Jump_tau11 + (C.^2 - S.^2).*Jump_tau12 + S.*C.*Jump_tau22;
            this.IsoInterface.Jump_tau_tt = (S.^2).*Jump_tau11 - 2*C.*S.*Jump_tau12 + (C.^2).*Jump_tau22;            
            
            function z = phiX1(y1)
                pt.y1_kv = y1;
                IP_h     = this.IDC.SubShapePtsCart(pt);
                z        = IP_h*phi;
            end   
        end
        function hatL = FitSliplength(this)
           
            thetaEq = this.optsPhys.thetaEq;
            Ca      = 3/4*this.optsPhys.Cak;
            Cn      = this.optsPhys.Cn;
                        
            L       = this.optsNum.PhysArea.y2Max;            
            
            if(isfield(this.optsPhys,'l_diff'))
                l_diff  = this.optsPhys.l_diff;
            else
                l_diff = 1;
            end
            fitPos  = min(max(10*Cn,12*l_diff),L*4/5);%L - min(L/2,2*l_diff);            
            x2M2    = this.IDC.CompSpace2(fitPos);
            IP05    = barychebevalMatrix(this.IDC.Pts.x2,x2M2);
            %Diff    = barychebdiff(this.IDC.Pts.y2,1);
                        
            theta_05L = IP05*this.IsoInterface.theta;
            %dthetaL = Diff.Dx(end,:)*theta;
            
            opts = optimoptions('fsolve','TolFun',1e-8,'TolX',1e-6);
            z = fsolve(@f,1,opts);
            
            hatL = z(1);            
            
            disp(['hat L = ',num2str(hatL)]);
            
            function y = f(x)
                hL = x(1);                
                
                GM   = (Ca*log(fitPos/hL)+GHR_lambdaEta(thetaEq,1));                                                
                y    = GM - GHR_lambdaEta(theta_05L,1);
%                 GM1  = GHR_Inv(Ca*log(L/lam)+GHR_lambdaEta(thetaEq,1),1);
%                 GM2  = GHR_Inv(Ca*log(L/2/lam)+GHR_lambdaEta(thetaEq,1),1);
%                 
%                 y(1) = (GM1 + Ca*C) - thetaL;
%                 y(2) = (GM2 + Ca*C) - theta_05L;
                %y(2) = f_stokes(GM1,1)*Ca/L - dthetaL;                
            end
        end
        function spt = FindStagnationPoint(this,iguess1,iguess2)
            
            uv = this.uv;
            
            if(nargin < 2)
                iguess1 = [-2,2];
            end
            if(nargin < 3)
                iguess2 = [-2,this.IDC.y2Max-2];
            end
            
            fsolveOpts = optimset('Display','off');
            [spt1,~,flag1] = fsolve(@ValueAtXY,iguess1,fsolveOpts);                        
            if(flag1 >= 1)
                spt.y1_kv = spt1(1);
                spt.y2_kv = spt1(2);
                disp(['Stagnation point found at (',num2str(spt1(1)),',',num2str(spt1(2)),')']);
            else
                spt.y1_kv = [];
                spt.y2_kv = [];
            end
                        
            [spt2,~,flag2] = fsolve(@ValueAtXY,iguess2,fsolveOpts);
            if(flag2 >= 1)
                spt.y1_kv = [spt.y1_kv;spt2(1)];
                spt.y2_kv = [spt.y2_kv;spt2(2)];
                disp(['Stagnation point found at (',num2str(spt2(1)),',',num2str(spt2(2)),')']);
            end
            
            if((flag1 < 1) && (flag2 < 1))
                disp('No stagnation point found');
                spt = [];
            else
                this.StagnationPoint = spt;                
            end
                                    
            function z = ValueAtXY(xy)
                pt.y1_kv = xy(1);
                pt.y2_kv = xy(2);
                IP  = this.IDC.SubShapePtsCart(pt);
                z   = [IP*uv(1:end/2);IP*uv(1+end/2:end)];
            end
        end                        
        
        %Plotting & Analysis functions       
        function AnalyzeScalarQuantity(this,f,interval)
            if(nargin == 2)
                interval = [-10,10];
            end
            y2Max = this.optsNum.PhysArea.y2Max;
            noCuts = 5;
            optsC  = {'b','g','m','k','r'};
            
            figure('Position',[0 0 1000 800]);
            subplot(2,1,1);
            leg = {};
            for i= 0:1:(noCuts-1)
                y2 = y2Max*i/(noCuts-1);
                this.IDC.plotLine(interval,[y2 y2],f,struct('color',optsC{i+1}));  hold on;
                
                IP       = this.IDC.SubShapePtsCart(struct('y1_kv',-inf,'y2_kv',y2));
                plot(interval,[1,1]*(IP*f),[optsC{i+1},'--']); hold on;
                IP       = this.IDC.SubShapePtsCart(struct('y1_kv',inf,'y2_kv',y2));
                plot(interval,[1,1]*(IP*f),[optsC{i+1},'-.']);  hold on;
                
                leg{end+1} = ['y2 = ',num2str(y2)];
                leg{end+1} = '';
                leg{end+1} = ['y2 = ',num2str(y2),'y1=-inf'];
                leg{end+1} = ['y2 = ',num2str(y2),'y1=inf'];
            end
            legend(leg,'Location','eastoutside');
            
            subplot(2,1,2);
            this.IDC.plotLine([-inf,-inf],[0 y2Max],f,struct('color','r')); hold on;
            this.IDC.plotLine([inf,inf],[0 y2Max],f,struct('color','b')); hold on;
            
            this.IDC.plotLine([-4 -4],[0 y2Max],f);  hold on;
            this.IDC.plotLine([0 0],[0 y2Max],f);  hold on;
            this.IDC.plotLine([4 4],[0 y2Max],f);  hold on;
            
            legend({'-inf','inf'},'Location','eastoutside');
        end
        function PlotResultsMu(this)                         
            figure('Position',[0 0 800 600],'color','white');
            this.IDC.plot(this.mu,'contour');             
            PlotU(this); hold on;             
            SaveCurrentFigure(this,'ChemPot');            
        end
        function PlotResultsPhi(this,opts)
            if(nargin < 2)
                opts = {};
            end
            if(IsOption(opts,'newFigure'))
                figure('Position',[0 0 800 600],'color','white');
            end
            
       %     PlotU(this); hold on;             
            %this.IDC.plot(this.phi,'contour');                 
                    
            optsPlot.nContours = -0.9;
            optsPlot.linecolor = 'b';
            optsPlot.linestyle = '--';
            this.IDC.plot(this.phi,'contour',optsPlot);  hold on;  

            optsPlot.nContours = 1e-10;
            optsPlot.linecolor = [0 0.75 0];
            this.IDC.plot(this.phi,'contour',optsPlot);  hold on;  

            optsPlot.nContours = 0.9;
            optsPlot.linecolor = 'r';
            this.IDC.plot(this.phi,'contour',optsPlot);  hold on;  
                        
            hold on;
            
            %if(~isempty(this.IsoInterface.h))
%                plot(this.IsoInterface.h,this.IDC.Pts.y2,...
               %                                     'k','linewidth',3);
%            end
      %      if(IsSeppecher(this))
%                PlotSeppecherSolution(this);
%            end
            if(IsOption(opts,'save'))
                SaveCurrentFigure(this,'Density');            
            end
        end        
        function PlotSeppecherSolution(this)                                                                 
           % PlotU(this);            hold on;                               
            y2Max   = this.optsNum.PhysArea.y2Max;
            PtsCart = this.IDC.GetCartPts;
            plot([this.deltaX (this.deltaX+y2Max/tan(this.theta))],[0 y2Max],'k--','linewidth',2.5);
            
            %u_flow         = GetSeppecherSolutionCart_Blurred([PtsCart.y1_kv - deltaX,...
%                                                PtsCart.y2_kv],1,0,0,theta);                  
                                    
            uSepp = GetSeppecherSolutionCart_Blurred([PtsCart.y1_kv - this.deltaX,...
                                                PtsCart.y2_kv],1,0,0,this.theta);
            PlotU(this,uSepp,struct('color','m','linewidth',2.5));
        end         
        function AddStreamlines(this,uv)
            if(nargin == 1)
                uv = this.uv;
            end
            for i = 1:3
                [y10,y20] = ginput(1);   
                this.IDC.plotStreamlines(uv,y10,y20); %IDC.plotFlux(u_flow)(mu);
            end
        end
        function [y1M,y2M,fl_y1,fl_y2,startMask1,startMask2] = PlotU(this,uv,y1Pts,y2Pts,opts) 
            
            n2 = 20;
            n1 = 20;
            if((nargin<2) || isempty(uv))
                if(isempty(this.uv))
                    return;
                else
                    uv = this.uv;
                end
            end
            y2Min = this.optsNum.PlotArea.y2Min;
            y2Max = this.optsNum.PlotArea.y2Max;
            
            y1Min = this.optsNum.PlotArea.y1Min;
            y1Max = this.optsNum.PlotArea.y1Max;
                        
            y2L = y2Min + (y2Max-y2Min)*(0:n2-2)'/(n2-1);            
            y1L = y1Min + (y1Max-y1Min)*(0:n1-1)'/(n1-1);            
            
            startPtsy1    = [y1Max*ones(size(y2L))-0.1;...
                             y1Min*ones(size(y2L))+0.1;...
                             y1L];
            startPtsy2    = [y2L;y2L;y2Max*ones(size(y1L))];
            
            if(nargin >= 4)
               startPtsy1 = [startPtsy1;y1Pts]; 
               startPtsy2 = [startPtsy2;y2Pts];
            end
            
            if((nargin >= 5))
                [y1M,y2M,fl_y1,fl_y2,startMask1,startMask2] = this.IDC.plotStreamlines(uv,startPtsy1,startPtsy2,opts); %IDC.plotFlux(u_flow)(mu);
            else
                [y1M,y2M,fl_y1,fl_y2,startMask1,startMask2] = this.IDC.plotStreamlines(uv,startPtsy1,startPtsy2); %IDC.plotFlux(u_flow)(mu);
            end            
                                   
        end      
        function PlotStagnationPoint(this)
            sp = this.StagnationPoint;
            if(~isempty(sp))          
                hold on;
                plot(sp.y1_kv,sp.y2_kv,'or','MarkerFaceColor','r','MarkerSize',10); 
                hold on;
            end
        end
        function PlotInterfaceAnalysis(this)            
            thetaEq = this.optsPhys.thetaEq;
            Ca      = 3/4*this.optsPhys.Cak;            
            
            y2      = this.IDC.Pts.y2;%IsoInterface.y2;
            L       = this.optsNum.PhysArea.y2Max;
                        
            thetaY2 = this.IsoInterface.theta;
                                                
            figure('color','white','Position',[0 0 800 800]);                        
            if(IsSeppecher(this))                                                
                y2P     = y2(2:end);
                hatL    = this.IsoInterface.hatL;
                theta_L = GHR_Inv(Ca*log(y2P/hatL)+GHR_lambdaEta(thetaEq,1),1);

                plot(y2,thetaY2*180/pi,'ok','linewidth',2); hold on;
                plot(y2P,theta_L*180/pi,'m','linewidth',2); hold on;                
                ylabel('$\theta$','Interpreter','Latex','fontsize',20);
            else
                A = -0.34; k = 0.2; 
	            analytic_BriantYeomans = A* (-(1+exp(-k*L)) + (exp(-k*y2) + exp(-k*(L-y2))));            

                plot(y2,cos(thetaY2),'k','linewidth',2); hold on;
                plot(y2,analytic_BriantYeomans,'k--','linewidth',2);
                ylabel('$\cos(\theta)$','Interpreter','Latex','fontsize',20);
            end
            
            xlabel('$y_2$','Interpreter','Latex','fontsize',20);            
            SaveCurrentFigure(this,'InterfaceSlope');           
            
            legendstr = {};
            figure('color','white','Position',[0 0 800 800]);                        
            plot(y2,this.IsoInterface.kappa/2,'k-o','linewidth',1.5); legendstr(end+1) = {'kappa/2'}; hold on;
            plot(y2,this.IsoInterface.mu,'m-o','linewidth',1.5);    legendstr(end+1) = {'mu'};                        
            plot(y2,this.IsoInterface.p,'r-o','linewidth',1.5);     legendstr(end+1) = {'p'};
            plot(y2,this.IsoInterface.u_n,'b-o','linewidth',1.5);   legendstr(end+1) = {'u_n'};            
            plot(y2,this.IsoInterface.u_t,'g-o','linewidth',1.5);   legendstr(end+1) = {'u_t'};
                        
            ylim([min(this.IsoInterface.u_n),max(this.IsoInterface.p)]*1.2);
            set(gca,'linewidth',1.5);
            xlabel('$y_2$','Interpreter','Latex','fontsize',20);            
            legend(legendstr,'Location','northEast');%,'Orientation','horizontal');
            SaveCurrentFigure(this,'InterfaceValues');           
        end
        function PlotResults(this)                                            
            PlotResultsPhi(this);                       
            PlotResultsMu(this);
            PlotErrorIterations(this);                             
            PlotInterfaceAnalysis(this);
        end        
        function PlotErrorIterations(this)           
            disp('*** Results ***');
            if(~isempty(this.theta))
                disp(['theta = ',num2str(this.theta*180/pi),' [deg]']);
            end
            if(isfield(this.errors,'aIter'))
                disp(['a = ',num2str(this.errors.aIter(end))]);
            end            

            legStr = {};
            figure('color','white');
            if(isfield(this.errors,'errorIterations') && ~isempty(this.errors.errorIterations))
                disp(['Max error of equations excluding boundaries: ',num2str(this.errors.errorIterations(end))]);
                semilogy(this.errors.errorIterations,'ro','MarkerFaceColor','r'); hold on;
                legStr{end+1} = 'Error';
            end
            if(isfield(this.errors,'errorAverage'))
                semilogy(this.errors.errorAverage,'ko','MarkerFaceColor','k');
                legStr{end+1} = 'Average error';
            end
            if(length(legStr) > 1)
                legend(legStr);            
            end
            xlabel('Iteration');
            ylabel('Error');
                        
            SaveCurrentFigure(this,'ErrorIterations');
        end
                
        function SaveCurrentFigure(this,filename,foldername)
            
            [~,fn]   = fileparts(this.filename);
            filename = [fn,'_',filename];            
            if(nargin>2)
                filename = [foldername filesep filename];
            end
            SaveCurrentFigure@Computation(this,filename);            
        end
        %Old                
        function theta  = FindInterfaceAngle(this,phi)

            y2M = 0; y2P = 10;
          %  fsolveOpts   = optimset('Display','off');                

            pt.y2_kv  =  y2M;
            y1CartStart = fsolve(@phiX1,-1);%,fsolveOpts);

            pt.y2_kv  = y2P;
            y1CartEnd = fsolve(@phiX1,2);%,fsolveOpts);                

            alpha = atan((y1CartStart-y1CartEnd)/(y2P-  y2M));
            theta = alpha + pi/2;
            
            disp(['Contact angle: ',num2str(theta*180/pi),'[deg].']);

            function z = phiX1(y1)
                pt.y1_kv = y1;
                IP       = this.IDC.SubShapePtsCart(pt);
                z        = IP*phi;
            end    
        end          
   end
    
end