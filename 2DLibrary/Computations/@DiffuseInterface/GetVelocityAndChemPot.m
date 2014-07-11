function [mu,uv,A,b] = GetVelocityAndChemPot(this,rho,D_B,theta)
%% Equations solved:
%
% Continuity: div(rho*uv) = 0
% Momentum:        div(m) = 0
%   - (rho+rho_m)*grad(mu) +
%          ... +grad(rho)*(W'(rho) - Cn*Lap(rho) - mu)
%              + Cak*((zeta + eta/3)*grad(div(uv)) + eta*Lap(uv) ) = 0
%
% with boundary conditions
% 
% (BC1) uv = [Uwall;1]   for Ind.bound & ~Ind.top
% (BC2) uv = U_Seppecher(D_B,theta) for Ind.top
% (BC3) mu = div(div(m)) 
%          = (zeta + 4/3*eta)*Lap(div(uv))  
%            - Lap( (rho+rho_m)*grad(mu)) 
%            + div((W'(rho) - Cn*Lap(rho))*grad(rho))  for Ind.top | Ind.bottom
% (BC4) mu = 0 for Ind.left
% (BC5) int( m_{1,2} ,y_1=-inf..inf,y2=y2Max) 
%     - int( m_{1,2} ,y_1=-inf..inf,y2=0) = [mu*(rho+rho_m)*y2Max]_y1=inf 
%                                         - [mu*(rho+rho_m)*y2Max]_y1=-inf

    M       = this.IC.M;
    Int     = this.IC.Int;
    Ind     = this.IC.Ind;
    PtsCart = this.IC.GetCartPts();
    D_A     = this.optsPhys.D_A;
    UWall   = this.optsPhys.UWall;
    Cn      = this.optsPhys.Cn;
    Cak     = this.optsPhys.Cak;
	eta     = this.optsPhys.eta;
    zeta    = this.optsPhys.zeta;
    rho_m   = this.optsPhys.rho_m;
    Diff    = this.IC.Diff;
    y2Max   = this.optsNum.PhysArea.y2Max;
    
    IBB     = [Ind.bound;Ind.bound];  
    F       = false(M,1);    FF      = false(2*M,1);
    T       = true(M,1);     TT      = true(2*M,1);
    A       = eye(3*M);  
    b       = zeros(3*M,1);        
        
    %% Continuity and Momentum Equation
    [Af,bf]                  = ContMom_DiffuseInterfaceSingleFluid(this,rho);
    A([~Ind.bound;~IBB],:)   = Af([~Ind.bound;~IBB],:);   
    b([~Ind.bound;~IBB])     = bf([~Ind.bound;~IBB]);

    %% BC1 and BC2
    uwall                                    = [UWall*ones(M,1) ; zeros(M,1);];
    uvBound(repmat(Ind.bound &~Ind.top,2,1)) =  uwall(repmat(Ind.bound &~Ind.top,2,1));                       
    uvBound([Ind.top;Ind.top]) = GetSeppecherSolutionCart(...
                                        [PtsCart.y1_kv(Ind.top),...
                                         PtsCart.y2_kv(Ind.top)],UWall,D_A,D_B,theta);
    b([F;IBB])                 = uvBound(IBB);           

    %% BC3
    ys             = DoublewellPotential(rho,Cn) - Cn*Diff.Lap*rho;    
    Cmu            = -diag(rho + rho_m)*Diff.Lap - diag(Diff.Lap*rho) ...
                     - 2*diag(Diff.Dy1*rho)*Diff.Dy1 ...
                     - 2*diag(Diff.Dy2*rho)*Diff.Dy2;
    Cuv            = Cak*(zeta + 4/3*eta)*Diff.LapDiv; 
    C              = [Cmu,Cuv];
    bBound         = - repmat(ys,2,1).*(Diff.grad*rho); 

    A([Ind.top|Ind.bottom;FF],:)  = C(Ind.top|Ind.bottom,:);                        
    b([Ind.top|Ind.bottom;FF])    = Diff.div(Ind.top|Ind.bottom,:)*bBound;
%    A([Ind.bottom;FF],:)  = C(Ind.bottom,:);                        
%    b([Ind.bottom;FF])    = Diff.div(Ind.bottom,:)*bBound;
    
    
    % BC3' [cos(theta),sin(theta)]^T  grad(mu) = 0
%     E               = eye(M);
% 	ETop            = E(Ind.top,:);
% 	topDirection    = [cos(theta)*ETop,sin(theta)*ETop];
% 	A([Ind.top;FF],[T;FF]) = topDirection*Diff.grad; 
    
    %% BC4
    A([Ind.left;FF],:)             = 0;        
    A([Ind.left;FF],[Ind.left;FF])  = eye(sum(Ind.left)); 
    A([Ind.left;FF],[Ind.right;FF]) = eye(sum(Ind.right)); 
    %A([Ind.left;FF],[Ind.left;FF]) = eye(sum(Ind.left)); 
    %b([Ind.left;FF])               = DoublewellPotential(rho(Ind.left),Cn) - Cn*Diff.DDy2(Ind.left,:)*rho;

    %% BC5
    OR                 = ones(sum(Ind.right),1);        
    IntPathUpLow       = this.IC.borderTop.IntNormal + this.IC.borderBottom.IntNormal; %?? => to check!!
    [Tt11,Tb11]        = FullStressTensorIJ(this,rho,1,1);
    [Tt12,Tb12]        = FullStressTensorIJ(this,rho,1,2);            
    TtHor              = [Tt11;Tt12];        
    TbHor              = [Tb11;Tb12];        

    A([Ind.right;FF],:)                = (OR*IntPathUpLow)*TtHor;
    A([Ind.right;FF],[Ind.right;FF])   = A(Ind.right,[Ind.right;FF]) ...
                                         - diag(rho(Ind.right)+rho_m)*y2Max;
    A([Ind.right;FF],[Ind.left;FF])   = A(Ind.right,[Ind.left;FF])  ...
                                        + diag(rho(Ind.left)+rho_m)*y2Max;
    b([Ind.right;FF])                 = -IntPathUpLow*TbHor;
        
    %% Solve Equation
    x               = A\b;

    disp(['Error: ',num2str(max(abs(A*x-b)))]);

    mu              = x(1:M);
    uv              = x(1+M:end);

    [At,bt] = Div_FullStressTensor(this,rho);
    disp(['Error of divergence of stress tensor: ',num2str(max(abs(At*[mu;uv] + bt)))]);       

end