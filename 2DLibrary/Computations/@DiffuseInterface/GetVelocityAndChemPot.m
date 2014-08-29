function [mu,uv,A,b,a] = GetVelocityAndChemPot(this,rho,theta)
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
%            - Lap( (rho+rho_m)*mu) 
%            + div((W'(rho) - Cn*Lap(rho))*grad(rho))  for Ind.top | Ind.bottom
% (BC4) mu = 0 for Ind.left
% (BC5) int( m_{1,2} ,y_1=-inf..inf,y2=y2Max) 
%     - int( m_{1,2} ,y_1=-inf..inf,y2=0) = [mu*(rho+rho_m)*y2Max]_y1=inf 
%                                         - [mu*(rho+rho_m)*y2Max]_y1=-inf

    M       = this.IC.M;
    Int     = this.IC.Int;
    Ind     = this.IC.Ind;
    PtsCart = this.IC.GetCartPts();
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
    
    OR                 = ones(sum(Ind.right),1);        
    IntPathUpLow       = this.IC.borderTop.IntSc - this.IC.borderBottom.IntSc; 
    [Tt11,Tb11]        = FullStressTensorIJ(this,rho,1,1);
    [Tt12,Tb12]        = FullStressTensorIJ(this,rho,1,2);            
    [Tt22,Tb22]        = FullStressTensorIJ(this,rho,2,2);
    [~,WFull]          = DoublewellPotential(rho,Cn);
    
        
    %% Continuity and Momentum Equation
    [Af,bf]                  = ContMom_DiffuseInterfaceSingleFluid(this,rho);
    A([~Ind.bound;~IBB],:)   = Af([~Ind.bound;~IBB],:);   
    b([~Ind.bound;~IBB])     = bf([~Ind.bound;~IBB]);

    %% BC1 and BC2
    uwall                                    = [UWall*ones(M,1) ; zeros(M,1);];
    uvBound(repmat(Ind.bound &~Ind.top,2,1)) =  uwall(repmat(Ind.bound &~Ind.top,2,1));                       
    %uvBound([Ind.top;Ind.top]) = GetSeppecherSolutionCart(...
%                                        [PtsCart.y1_kv(Ind.top),...
%                                         PtsCart.y2_kv(Ind.top)],UWall,D_A,D_B,theta);                                     
    [uvBound([Ind.top;Ind.top]),a] = CorrectVelocityProfile(this,theta,rho);
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
        
%     nV = (Diff.Dy1(Ind.top|Ind.bottom,:)*rho).^2 + (Diff.Dy2(Ind.top|Ind.bottom,:)*rho).^2;
%     A([Ind.top|Ind.bottom;FF],[T;FF])  = (diag((Diff.Dy1(Ind.top|Ind.bottom,:)*rho)./nV)*Diff.Dy2(Ind.top|Ind.bottom,:)...
%                                         - diag((Diff.Dy2(Ind.top|Ind.bottom,:)*rho)./nV)*Diff.Dy1(Ind.top|Ind.bottom,:));
%     b([Ind.top|Ind.bottom;FF])    = 0;
        
    %% BC4
    A([Ind.left;FF],[T;FF]) = Diff.Dy2(Ind.left,:);
    
    A([Ind.left&Ind.top;FF],:)                      = 0;        
    A([Ind.left&Ind.top;FF],[Ind.left&Ind.top;FF])  = eye(sum(Ind.left&Ind.top)); 
    %A([Ind.left&Ind.top;FF],[Ind.right&Ind.top;FF]) = eye(sum(Ind.right&Ind.top)); 
    
    %b(Ind.left&Ind.top) = 0;
    
    %A([Ind.left&Ind.top;FF],:) = 0;
    %A([Ind.left&Ind.top;FF],:) = IntPathUpLow*Tt22;
    %b([Ind.left&Ind.top;FF])   = -IntPathUpLow*Tb22;    
    
    %% BC5    
    [~,WL] = DoublewellPotential(rho(Ind.left&Ind.top),Cn);
    [~,WR] = DoublewellPotential(rho(Ind.right&Ind.top),Cn);
    
    indBC5 = Ind.right & Ind.top;
    
    A([Ind.right;FF],[T;FF]) = Diff.Dy2(Ind.right,:);
        
    
     A([indBC5;FF],[T;FF]) = 0;
     A([indBC5;FF],:)                = IntPathUpLow*Tt12;
     A([indBC5;FF],[T;FF])           = A(indBC5,[T;FF]) ...
                                         - this.IC.borderRight.IntSc*diag(rho+rho_m);
     A([indBC5;FF],[T;FF])           = A(indBC5,[T;FF])  ...
                                         + this.IC.borderLeft.IntSc*diag(rho+rho_m);
     b([indBC5;FF])                  = -IntPathUpLow*Tb12 + (WL - WR)*y2Max;
        
%     
%    [~,WL] = DoublewellPotential(rho(Ind.left),Cn);
%    [~,WR] = DoublewellPotential(rho(Ind.right),Cn);    
%
%     A([Ind.right;FF],:)                = (OR*IntPathUpLow)*Tt12;
%     A([Ind.right;FF],[Ind.right;FF])   = A(Ind.right,[Ind.right;FF]) ...
%                                          - diag(rho(Ind.right)+rho_m)*y2Max;
%     A([Ind.right;FF],[Ind.left;FF])   = A(Ind.right,[Ind.left;FF])  ...
%                                         + diag(rho(Ind.left)+rho_m)*y2Max;
%     b([Ind.right;FF])                 = -IntPathUpLow*Tb12 + (WL - WR)*y2Max;
    
    %% BC6                
    %A([Ind.left;FF],:)                = (OR*IntPathUpLow)*Tt22;
    %b([Ind.left;FF])                  = -IntPathUpLow*Tb22;    
    
        
    %% Solve Equation
    x               = A\b;

    disp(['Error: ',num2str(max(abs(A*x-b)))]);

    mu              = x(1:M);
    uv              = x(1+M:end);

    [At,bt] = Div_FullStressTensor(this,rho);
    [maxError,j] = max(abs(At*[mu;uv] + bt));
    j = mod(j,M);
    disp(['Error of divergence of stress tensor: ',num2str(maxError),...
                ' at y_1 = ',num2str(this.IC.GetCartPts.y1_kv(j)),...
                ' at y_2 = ',num2str(this.IC.GetCartPts.y2_kv(j))]);
    
    %Accuracy of force balance normal to the substrate:
    errorNormal = IntPathUpLow*Tt22*[mu;uv] + IntPathUpLow*Tb22;
    disp(['Error of force balance normal to substate: ',num2str(errorNormal)])    

    
    %Accuracy of force balance normal to the substrate:    
    h         = WFull - (rho+rho_m).*mu;
    
    errorParallel = IntPathUpLow*Tt12*[mu;uv] + IntPathUpLow*Tb12 + (this.IC.borderRight.IntSc - this.IC.borderLeft.IntSc)*h;
    disp(['Error of force balance parallel to substate: ',num2str(errorParallel)])    
    
end