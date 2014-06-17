function [mu,uv,A,b] = GetVelocityAndChemPot(this,rho,D_B,theta)

    M       = this.IC.M;
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
    F       = false(M,1);
    FF      = false(2*M,1);
    A       = eye(3*M);  
    b       = zeros(3*M,1);
    
    y10     = 0; %this has to be adapted for next iteration
    
	uwall                                    = [UWall*ones(M,1) ; zeros(M,1);];
    uvBound(repmat(Ind.bound &~Ind.top,2,1)) =  uwall(repmat(Ind.bound &~Ind.top,2,1));               
    uvBound([Ind.top;Ind.top]) = ...
                GetSeppecherSolutionCart([PtsCart.y1_kv(Ind.top)- y10,PtsCart.y2_kv(Ind.top)],UWall,D_A,D_B,theta);

    [Af,bf]                  = ContMom_DiffuseInterfaceSingleFluid(this,rho);
    A([~Ind.bound;~IBB],:)   = Af([~Ind.bound;~IBB],:);   
    b([~Ind.bound;~IBB])     = bf([~Ind.bound;~IBB]);

    b([F;IBB])               = uvBound(IBB);        

    %******************************************************
    %BC for chemical Potential
    %
    %1. Top and Bottom Boundary:
    %Boundary condition for chemical potential      

    %******** Take momentum eq normal to boundary ********
    %A([Ind.top;FF],:)    = Ind.normalTop*Af([F;TT],:);
    %b([Ind.top;FF])      = Ind.normalTop*bf([F;TT],:);

    %A([Ind.bottom;FF],:) = Ind.normalBottom*Af([F;TT],:);
    %b([Ind.bottom;FF])   = Ind.normalBottom*bf([F;TT],:);        

    %*****************************************************
    %*** Take divergence of momentum eq ***
    % ******** check    
    ys             = DoublewellPotential(rho,Cn) - Cn*Diff.Lap*rho;

    Cmu            = -Diff.Lap*diag(rho + rho_m);                   
    Cuv            = -Cak*(zeta + 4/3*eta)*(Diff.Lap*Diff.div);
    C              = [Cmu,Cuv];
    bBound         = - repmat(ys,2,1).*(Diff.grad*rho); 

    A([Ind.top|Ind.bottom;FF],:)  = C(Ind.top|Ind.bottom,:);                        
    b([Ind.top|Ind.bottom;FF])    = Diff.div(Ind.top|Ind.bottom,:)*bBound;
    %*****************************************************

    %2. Left boundary:
    %Set Chemical potential at -infinty to zero        
    A([Ind.left;FF],:)             = 0;        
    A([Ind.left;FF],[Ind.left;FF]) = eye(sum(Ind.left)); 
    b([Ind.left;FF])               = 0;

    %3. Right boundary:
    OR                 = ones(sum(Ind.right),1);        
    IntPathUpLow       = this.IC.borderTop.IntNormal + this.IC.borderBottom.IntNormal; %?? => to check!!
    [Tt11,Tb11]        = FullStressTensorIJ(this,rho,1,1);
    [Tt12,Tb12]        = FullStressTensorIJ(this,rho,1,2);            
    TtHor              = [Tt11;Tt12];        
    TbHor              = [Tb11;Tb12];        

    A(Ind.right,:)                = (OR*IntPathUpLow)*TtHor;
    A(Ind.right,[Ind.right;FF])   = A(Ind.right,[Ind.right;FF]) - ...
                                    diag(rho(Ind.right)+rho_m)*y2Max;
    b(Ind.right)                  = -IntPathUpLow*TbHor;
    
    %******************************************************
    x               = A\b;

    disp(['Error: ',num2str(max(abs(A*x-b)))]);

    mu              = x(1:M);
    uv              = x(1+M:end);

    [At,bt] = Div_FullStressTensor(this,rho);
    disp(['Error of divergence of stress tensor: ',num2str(max(abs(At*[mu;uv] + bt)))]);       

  %  figure;
  %  doPlots_SC_Path(InterpPathUpper,TtHor*[mu;uv]+TbHor);
  %  xlim([-40 40]);

end