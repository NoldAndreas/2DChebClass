function [p,uv,A,b,a] = GetVelocityAndChemPot(this,phi,theta)
%% Equations solved:
%
% Continuity: div(uv) = 0
% Momentum:    div(m) = 0
%         0 = -grad(p) + Lap*u + 1/Cak *mu*grad(phi)
%           = Lap*u + 1/Cak * div*PI
%
% PI = (-p*Cak + f_DW + Cn/2*|Grad(phi)|^2) I - Cn (Grad(phi)) X (Grad(phi))
%
% with boundary conditions
% 
% (BC1) uv = [Uwall;1]   for Ind.bound & ~Ind.top
% (BC2) uv = U_Seppecher(D_B,theta) for Ind.top
% (BC3) mu = div(div(m))
%          = -Lap(p) + 1/Cak*div(mu*grad(phi))
% (BC4) p = 0 for Ind.left
%
%
%
% (BC5) int( m_{1,2} ,y_1=-inf..inf,y2=y2Max) 
%     - int( m_{1,2} ,y_1=-inf..inf,y2=0) = [mu*(phi+phi_m)*y2Max]_y1=inf 
%                                         - [mu*(phi+phi_m)*y2Max]_y1=-inf

    M       = this.IC.M;    
    Ind     = this.IC.Ind;       
    Cn      = this.optsPhys.Cn;
    Cak     = this.optsPhys.Cak;
    Diff    = this.IC.Diff;
    y2Max   = this.optsNum.PhysArea.y2Max;
    
    IBB     = [Ind.bound;Ind.bound];  
    F       = false(M,1);    FF      = false(2*M,1);
    T       = true(M,1);     TT      = true(2*M,1);
    A       = eye(3*M);  
    b       = zeros(3*M,1);        
        
    IntPathUpLow       = this.IC.borderTop.IntSc - this.IC.borderBottom.IntSc; 
    [Tt11,Tb11]        = FullStressTensorIJ(this,phi,1,1);
    [Tt12,Tb12]        = FullStressTensorIJ(this,phi,1,2);            
    [Tt22,Tb22]        = FullStressTensorIJ(this,phi,2,2);
    [~,WFull]          = DoublewellPotential(phi,Cn);
    
        
    %% Continuity and Momentum Equation
    [Af,bf]                  = ContinuumMomentumEqs(this,phi);
    A([~Ind.left;~IBB],:)   = Af([~Ind.left;~IBB],:);   
    b([~Ind.left;~IBB])     = bf([~Ind.left;~IBB]);

    %% BC1 and BC2       
	[uvBound,a] = GetBoundaryCondition(this,theta,phi);    
	b([F;IBB])  = uvBound(IBB);

    %% BC3
%     ys             = DoublewellPotential(phi,Cn) - Cn*Diff.Lap*phi;    
%     Cmu            = -diag(phi + phi_m)*Diff.Lap - diag(Diff.Lap*phi) ...
%                      - 2*diag(Diff.Dy1*phi)*Diff.Dy1 ...
%                      - 2*diag(Diff.Dy2*phi)*Diff.Dy2;
%     Cuv            = zeros(M,2*M);
%     C              = [Cmu,Cuv];
%     bBound         = - repmat(ys,2,1).*(Diff.grad*phi); 
% 
%     A([Ind.top|Ind.bottom;FF],:)  = C(Ind.top|Ind.bottom,:);                        
%     b([Ind.top|Ind.bottom;FF])    = Diff.div(Ind.top|Ind.bottom,:)*bBound;
               
    %% BC4
    A([Ind.left;FF],[T;FF]) = Diff.Dy2(Ind.left,:);
    
    A([Ind.left&Ind.top;FF],:)                      = 0;        
    A([Ind.left&Ind.top;FF],[Ind.left&Ind.top;FF])  = eye(sum(Ind.left&Ind.top));    
    
    %% BC5    
    [~,WL] = DoublewellPotential(phi(Ind.left&Ind.top),Cn);
    [~,WR] = DoublewellPotential(phi(Ind.right&Ind.top),Cn);
    
     indBC5 = Ind.right & Ind.top;
    
     A([Ind.right;FF],[T;FF]) = Diff.Dy2(Ind.right,:);
        
    
     A([indBC5;FF],[T;FF]) = 0;
     A([indBC5;FF],:)      = IntPathUpLow*Tt12;
     A([indBC5;FF],[T;FF]) = A(indBC5,[T;FF]) ...
                                 - this.IC.borderRight.IntSc*Cak;
     A([indBC5;FF],[T;FF]) = A(indBC5,[T;FF])  ...
                                + this.IC.borderLeft.IntSc*Cak;
     b([indBC5;FF])        = -IntPathUpLow*Tb12 + (WL - WR)*y2Max;
             
     %reduce by ignoring velocities with y1>y1Max
     markRed = (abs(this.IC.GetCartPts.y1_kv) < this.optsNum.PhysArea.y1Max);
     markRedFull  = [T;markRed;markRed];
     ARed = A(markRedFull,markRedFull);
     bRed = b(markRedFull) - A(markRedFull,~markRedFull)*...
                                uvBound(~[markRed;markRed]);
     
    % Solve Equation
    %x               = A\b;
    x               = [ones(M,1);uvBound];
    x(markRedFull)  = ARed\bRed;

    disp(['Error: ',num2str(max(abs(A*x-b)))]);

    p               = x(1:M);
    uv              = x(1+M:end);

    this.p = p; this.uv = uv;

    [At,bt] = Div_FullStressTensor(this,phi);
    [maxError,j] = max(abs(At*[p;uv] + bt));
    j = mod(j,M);
    disp(['Error of divergence of stress tensor: ',num2str(maxError),...
                ' at y_1 = ',num2str(this.IC.GetCartPts.y1_kv(j)),...
                ' at y_2 = ',num2str(this.IC.GetCartPts.y2_kv(j))]);
    
    %Accuracy of force balance normal to the substrate:
    errorNormal = IntPathUpLow*Tt22*[p;uv] + IntPathUpLow*Tb22;
    disp(['Error of force balance normal to substate: ',num2str(errorNormal)])    

    
    %Accuracy of force balance normal to the substrate:    
    h         =  WFull - Cak*p;
    
    errorParallel = IntPathUpLow*Tt12*[p;uv] + IntPathUpLow*Tb12 + (this.IC.borderRight.IntSc - this.IC.borderLeft.IntSc)*h;
    disp(['Error of force balance parallel to substate: ',num2str(errorParallel)])    
    
end