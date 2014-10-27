function IterationStepFullProblem_Seppecher(this,noIterations)
    %Continuitiy: 0 = div(uv)
    %Momentum:    0 = -Grad(p) + (G+s*phi)*grad(phi)/Cak + Lap(uv)
    %Phasefield   0 = m*Lap(G + s*phi) - u*grad(phi)
    %ChemPot      0 = f_w' - s*phi - Cn*Lap(phi) - G
    %
    % 
    % (BC1) uv = uv_BC    
    % (BC2) nu*grad(phi) = 0
    % (BC3) nu*grad(G) = 0
    % (BC4) p = 0  at y1 = +/- infinity
    
    % A*[uv;phi;G;p] = b corresponds to momentum and continuity Eq. for given phasefield phi                           
    if(nargin == 1)
        noIterations = 20;
    end
    
    Cn             = this.optsPhys.Cn;
    Cak            = this.optsPhys.Cak;    
    zeta           = this.optsPhys.zeta;
    phi_m          = this.optsPhys.phi_m;    
    nParticles     = this.optsPhys.nParticles;
    
    PtsCart        = this.IC.GetCartPts();
    Diff           = this.IC.Diff;
    M              = this.IC.M;  
    Ind            = this.IC.Ind;
    IntSubArea     = this.IntSubArea;    
    y2Max          = this.optsNum.PhysArea.y2Max;
    
    IntPathUpLow   = this.IC.borderTop.IntSc - this.IC.borderBottom.IntSc; 
    IntNormalUp    = this.IC.borderTop.IntNormal;
    tauM           = Cak*(Diff.LapVec + (zeta + 1/3)*Diff.gradDiv);
    
	Dy12           = blkdiag(Diff.Dy1,Diff.Dy1);
	d_theta        = 0.01;
    
    lbC = Ind.left & Ind.bottom;
    rbC = Ind.right & Ind.bottom;

    Z    = zeros(M);
    IBB  = repmat(Ind.bound,2,1);  
    ITT  = repmat(Ind.top,2,1);  
    EYM  = eye(M);  EYMM  = eye(2*M);

    F       = false(M,1);   T       = true(M,1);    
    
    opts = struct('optsNum',this.optsNum,...
                  'optsPhys',this.optsPhys,...
                  'Comments',this.configName);
              
    if(isempty(this.mu))
        in = struct('initialGuess',[0;0;pi/2;GetInitialCondition(this)]);        
    else
        in = struct('initialGuess',[this.a;...
                                    this.deltaX;...
                                    this.theta;...
                                    this.uv;...
                                    this.phi;...
                                    this.mu]);
    end
    
    [res,~,Parameters] = DataStorage([],@SolveSingleFluid,...
                    opts,in);
	   
    this.uv       = res.uv;
    this.mu       = res.mu;
    this.phi      = res.phi;
    this.a        = res.a;
    this.deltaX   = res.deltaX;
    this.theta    = res.theta;
    
    this.filename = Parameters.Filename;
    this.errors.errorIterations = res.errorHistory;
    
    %     this.errors.errorIterations = res.errorIterations;
    
    
    function res = SolveSingleFluid(conf,in)
        [vec,errHistory] = NewtonMethod(in.initialGuess,@f,1e-7,noIterations,0.8);    
    
        %[uv;phi;mu] 
        res.a      = vec(1);
        res.deltaX = vec(2);
        res.theta  = vec(3);
        vec        = vec(4:end);
        res.uv     = vec([T;T;F;F]);
        res.phi    = vec([F;F;T;F]);
        res.mu     = vec([F;F;F;T]);
        
        res.errorHistory = errHistory;
    end    
    function [v_full,A_full] = f(z)
        
        %[uv;phi;G;p]         
        a      = z(1);
        deltaX = z(2);
        theta  = z(3);
        disp(['[a,deltaX,theta] = ',num2str(a),' , ',num2str(deltaX),' , ',num2str(theta*180/pi),'.']);
        
        z      = z(4:end);
        uv  = z([T;T;F;F]);
        phi = z([F;F;T;F]);
        G   = z([F;F;F;T]);        
        
        phiM2          = repmat(phi+phi_m,2,1);                   
            
        [v_cont,A_cont]  = Continuity(this,uv,phi,G);               
        [v_mom,A_mom]    = Momentum(this,uv,phi,G);      
        [v_mu,A_mu]      = ChemicalPotential(this,uv,phi,G);
        
       %% Boundary conditions [uv;phi;G]
        % (BC4.b) a*grad(phi) = 0        
        a_direction               = [cos(theta)*EYM,sin(theta)*EYM];            
        a_direction_theta         = [-sin(theta)*EYM,cos(theta)*EYM];            
        A_cont(Ind.top,:)           = 0;
        A_cont(Ind.top,[F;F;T;F])   = a_direction(Ind.top,:)*Diff.grad;
        
        A_contThree                             = zeros(M,3);
        A_contThree(Ind.top,[false,false,true]) = a_direction_theta(Ind.top,:)*(Diff.grad*phi);
        v_cont(Ind.top)                         = a_direction(Ind.top,:)*(Diff.grad*phi);                            
                       
        % (BC3.a) uv = uv_BC    
        %[uvBound,a]            = GetBoundaryCondition(this);%,theta,phi);           
        u_flow      = GetSeppecherSolutionCart_Blurred([PtsCart.y1_kv - deltaX,...
                                                PtsCart.y2_kv],1,0,0,theta);                  
                                            
        y1_Interface  = PtsCart.y1_kv - (deltaX + y2Max/tan(theta));
        a_corr        = (1 + a./(1+y1_Interface.^2));
        a_corr_a      = 1./(1+y1_Interface.^2);        
        a_corr_deltaX = 2*a*y1_Interface./(1+y1_Interface.^2).^2;
        a_corr_deltaX(Ind.left | Ind.right) = 0;
        a_corr_theta  = -a_corr_deltaX*y2Max/(sin(theta)).^2;
        a_corr_phi    = zeros(M,M);

%         a_corr                = (1 + a*(phi(lbC)-phi).^2.*(phi(rbC)-phi).^2);
%         a_corr_phi            = diag(-2*a*(phi(lbC)-phi).*(phi(rbC)-phi).^2 ...
%                                      -2*a*(phi(lbC)-phi).^2.*(phi(rbC)-phi));
%         a_corr_phi(:,lbC)     =  a_corr_phi(:,lbC) + 2*a*(phi(lbC)-phi).*(phi(rbC)-phi).^2;
%         a_corr_phi(:,rbC)     =  a_corr_phi(:,rbC) + 2*a*(phi(rbC)-phi).*(phi(lbC)-phi).^2;
%         a_corr_a              = ((phi(lbC)-phi).^2).*(phi(rbC)-phi).^2;
%         a_corr_theta          = zeros(M,1);
%         a_corr_deltaX         = zeros(M,1);

        uvBound        = u_flow .*repmat(a_corr,2,1);                
        uvBound_phi    = diag(u_flow)*repmat(a_corr_phi,2,1);                                 
        uvBound_a      = u_flow.*repmat(a_corr_a,2,1);               
        uvBound_deltaX = -(Dy12*u_flow).*repmat(a_corr,2,1)+u_flow .*repmat(a_corr_deltaX,2,1);        
        
        u_flow_PTheta    = GetSeppecherSolutionCart([PtsCart.y1_kv - deltaX,PtsCart.y2_kv],1,0,0,theta+d_theta);
        u_flow_d         = (u_flow_PTheta - u_flow)/d_theta;
        uvBound_theta    = u_flow_d.*repmat(a_corr,2,1) + u_flow.*repmat(a_corr_theta,2,1);
        
        A_mom(ITT,:)           = 0;
        A_mom(ITT,[T;T;F;F])   = EYMM(ITT,:);
        A_mom(ITT,[F;F;T;F])   = -uvBound_phi(ITT,:);
                
        A_mom_a(ITT)           = -uvBound_a(ITT);
        A_mom_deltaX(ITT)      = -uvBound_deltaX(ITT); 
        A_mom_theta(ITT)       = -uvBound_theta(ITT);        
                
        A_momThree             = [A_mom_a',A_mom_deltaX',A_mom_theta'];
        
        v_mom(ITT)             = uv(ITT) - uvBound(ITT);
        
        % (BC3.b)
        u_Wall = [ones(M,1);zeros(M,1)];
        
        A_mom(IBB & ~ITT,:)                = 0;
        A_mom(IBB & ~ITT,[IBB & ~ITT;F;F]) = eye(sum(IBB & ~ITT));
        v_mom(IBB & ~ITT)                  = uv(IBB & ~ITT) - u_Wall(IBB & ~ITT);       
%
%         A_mu(Ind.top,[T;T;F;F]) = Cak*(zeta+4/3)*Diff.LapDiv(Ind.top,:);
%         A_mu(Ind.top,[F;F;T;F]) = - Diff.Lap(Ind.top,:)*diag(G) ...
%                                   + Diff.div(Ind.top,:)*(diag(repmat(G,2,1))*Diff.grad);
% 
%         A_mu(Ind.top,[F;F;F;T]) =  - Diff.Lap(Ind.top,:)*(diag(phi+phi_m))...
%                                    + Diff.Dy1(Ind.top,:)*diag(Diff.Dy1*phi) ...
%                                    + Diff.Dy2(Ind.top,:)*diag(Diff.Dy2*phi) ;                                   
% 
%         A_muThree               = zeros(M,3);
% 
%         v_mu(Ind.top)   = Cak*(zeta+4/3)*Diff.LapDiv(Ind.top,:)*uv...
%                           - Diff.Lap(Ind.top,:)*(G.*(phi+phi_m)) ...
%                           + Diff.div(Ind.top,:)*(repmat(G,2,1).*(Diff.grad*phi));
%
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
        InterpMatchPos       = this.IC.SubShapePtsCart(...
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
%          A_theta              = zeros(1,4*M);
%          A_theta([F;F;lbC;F]) = -G(lbC)+fWP(lbC);
%          A_theta([F;F;F;lbC]) = -(phi(lbC)+phi_m);        
%          A_theta              = [0,0,0,A_theta];
%          v_theta              = -(phi(lbC)+phi_m)*G(lbC) + fW(lbC);

        A_theta = zeros(1,4*M);
        A_theta = [0,0,1,A_theta];    
        v_theta = theta - pi/2;
        
        A = [A_mom;A_mu;A_cont];
        v = [v_mom;v_mu;v_cont];    
        
        A_three = [A_momThree;zeros(M,3);A_contThree];
        
        A_full  = [A_a;...
                   A_deltaX;...
                   A_theta;...
                   [A_three,A]];
      
        v_full  = [v_a;...
                   v_deltaX;...
                   v_theta;...
                   v];                
        
        DisplayError(v_full);
    end    
    function DisplayError(errorFull)
        
        PrintErrorPos(errorFull(1),'consistent mass influx');
        PrintErrorPos(errorFull(2),'zero density at interface');        
        PrintErrorPos(errorFull(3),'chemical potential at -inf');
        
        error = errorFull(4:end);
        PrintErrorPos(error([F;F;F;T]),'continuity equation',this.IC.Pts);
        PrintErrorPos(error([T;F;F;F]),'y1-momentum equation',this.IC.Pts);
        PrintErrorPos(error([F;T;F;F]),'y2-momentum equation',this.IC.Pts);                                    
        PrintErrorPos(error([F;F;T;F]),'mu equation',this.IC.Pts);
    end

end